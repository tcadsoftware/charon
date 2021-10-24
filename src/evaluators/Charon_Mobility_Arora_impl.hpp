
#ifndef CHARON_MOBILITY_ARORA_IMPL_HPP
#define CHARON_MOBILITY_ARORA_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Names.hpp"

/*
Arora mobility model -- reference by Arora, Hauser, and Roulston,
"Electron and Hole Mobilities in Silicon as a Function of
Concentration and Temperature," IEEE Transactions on Electron
Devices, Vol. ED-29, pp.292-295, 1982.

This mobility model is a low-field model, contains effects from phonon scattering,
ionized impurity, and carrier-carrier scatterings. It does not differ majority from
minority mobility and does not include screening. The electron mobility expression
is good for Phosphorous-doped Silicon, 250-500 K and 1e13-1e20 cm^-3 doping.
The hole expresson is good for Boron-doped Silicon, 200-400 K and 1e15-1e20 cm^-3.
Both expressions have about 13% error w.r.t. reported experimental values.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Arora<EvalT, Traits>::
Mobility_Arora(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Retrieve edge data layout and cellType if isEdgedl = true
  isEdgedl = p.get<bool>("Is Edge Data Layout");
  RCP<DataLayout> output_dl = scalar;
  if (isEdgedl)
  {
    RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    output_dl = cellTopoInfo->edge_scalar;
    num_edges = output_dl->dimension(1);
    cellType = cellTopoInfo->getCellTopology();
  }

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize Arora mobility parameters, must be called here since
  // it uses the above variables in the code
  initMobilityParams(matName, mobParamList);

  // Evaluated field
  if (carrType == "Electron")
    mobility = MDField<ScalarT,Cell,Point>(n.field.elec_mobility, output_dl);
  else if (carrType == "Hole")
    mobility = MDField<ScalarT,Cell,Point>(n.field.hole_mobility, output_dl);
  this->addEvaluatedField(mobility);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // Depedent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  std::string name = "Arora_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Arora<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Reference temperature [K]
  ScalarT refTemp = 300.0;

  // FieldContainer to temporarily hold mobility values
  Kokkos::DynRankView<ScalarT,PHX::Device> tmpMob = Kokkos::createDynRankView(mobility.get_static_view(),"tmpMob",workset.num_cells, num_points);

  // Compute the mobility at IP or BASIS
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Obtain temperature [K]
      ScalarT lattTemp = latt_temp(cell,point)*T0;

      // Compute temperature-dependent parameters
      ScalarT normTemp = lattTemp / refTemp;
      ScalarT minMob = mumin * pow(normTemp, ex1);
      ScalarT diffMob = mumax * pow(normTemp, ex2);
      ScalarT refConc = nref * pow(normTemp, ex3);
      ScalarT alpha = exc * pow(normTemp, ex4);

      // Obtain doping concentration
      const ScalarT& Na = acceptor(cell,point); // scaled conc.
      const ScalarT& Nd = donor(cell,point);
      ScalarT Ntot = (Na + Nd) * C0;        // unscale the conc.

      ScalarT mobValue = minMob + diffMob / (1.0 + pow(Ntot/refConc,alpha) );
      tmpMob(cell,point) = mobValue / Mu0;
    }
  }

  // Want mobility available at the center of primary edges
  if (isEdgedl)
  {
    // compute mobility at the center of primary edge
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        // edge low-field mobility
        mobility(cell,edge) = (tmpMob(cell,node0) + tmpMob(cell,node1))/2.0;
      }
    }
  }  // end of if block

  // Want mobility available at IP or BASIS
  else
  {
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
     for (int point = 0; point < num_points; ++point)
       mobility(cell,point) = tmpMob(cell,point);
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_Arora<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for Arora electron mobility model
  if (carrType == "Electron")
  {
    // Retrieve parameters from charon::Material_Properties by default
    mumax = matProperty.getPropertyValue(matName, "Arora Electron mumax");
    mumin = matProperty.getPropertyValue(matName, "Arora Electron mumin");
    nref = matProperty.getPropertyValue(matName, "Arora Electron nref");
    exc = matProperty.getPropertyValue(matName, "Arora Electron nref_exp");
    ex1 = matProperty.getPropertyValue(matName, "Arora Electron exp1");
    ex2 = matProperty.getPropertyValue(matName, "Arora Electron exp2");
    ex3 = matProperty.getPropertyValue(matName, "Arora Electron exp3");
    ex4 = matProperty.getPropertyValue(matName, "Arora Electron exp4");
  }

  // Set up parameters for Arora hole mobility model
  else if (carrType == "Hole")
  {
    // Retrieve parameters from charon::Material_Properties by default
    mumax = matProperty.getPropertyValue(matName, "Arora Hole mumax");
    mumin = matProperty.getPropertyValue(matName, "Arora Hole mumin");
    nref = matProperty.getPropertyValue(matName, "Arora Hole nref");
    exc = matProperty.getPropertyValue(matName, "Arora Hole nref_exp");
    ex1 = matProperty.getPropertyValue(matName, "Arora Hole exp1");
    ex2 = matProperty.getPropertyValue(matName, "Arora Hole exp2");
    ex3 = matProperty.getPropertyValue(matName, "Arora Hole exp3");
    ex4 = matProperty.getPropertyValue(matName, "Arora Hole exp4");
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite parameters when specified by users
  if (mobParamList.isParameter("mumax"))
    mumax = mobParamList.get<double>("mumax");
  if (mobParamList.isParameter("mumin"))
    mumin = mobParamList.get<double>("mumin");
  if (mobParamList.isParameter("nref"))
    nref = mobParamList.get<double>("nref");
  if (mobParamList.isParameter("nref_exp"))
    exc = mobParamList.get<double>("nref_exp");
  if (mobParamList.isParameter("exp1"))
    ex1 = mobParamList.get<double>("exp1");
  if (mobParamList.isParameter("exp2"))
    ex2 = mobParamList.get<double>("exp2");
  if (mobParamList.isParameter("exp3"))
    ex3 = mobParamList.get<double>("exp3");
  if (mobParamList.isParameter("exp4"))
    ex4 = mobParamList.get<double>("exp4");

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Arora<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Material Name", "?");
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Arora", "Arora low-field mobility model");
  p->sublist("Mobility ParameterList").set<double>("mumax", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mumin", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("nref", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("nref_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("exp1", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("exp2", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("exp3", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("exp4", 0., "[1]");

  p->set("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
