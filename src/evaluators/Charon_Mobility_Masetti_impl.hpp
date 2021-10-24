
#ifndef CHARON_MOBILITY_MASETTI_IMPL_HPP
#define CHARON_MOBILITY_MASETTI_IMPL_HPP

#include <cmath>
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Names.hpp"


/*
Masetti mobility model -- reference by Guido Masetti,Maurizio Severi, and Sandro Solmi
"Modeling of Carrier Mobility Against Carrier Concentration in Arsenic-,
Phosphorus-, and Boron-Doped Silicon," IEEE Transactions on Electron
Devices, Vol. ED-30, pp.764-769, 1983.

This mobility model is a low-field model, contains effects from phonon scattering,
and ionized impurity. It has different parameters for Arsenic- and Phosphorus-doped
n-type Silicon, but it does not differ majority from minority mobility and does not
include screening. The expressions are good for T = 300 K and a very wide range
of impurity concentration, i.e., 1e13 to 5e21 cm^-3 for Phosphorus, 1e13 to
4e21 cm^-3 for Arsenic, and 1e14 to 1.2e21 cm^-3 for Boron.

To include temperature dependence, replace mumax in Eq.(1) of the reference with
mumax*(T/300K)^(gamma).
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Masetti<EvalT, Traits>::
Mobility_Masetti(
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

  // Initialize Masetti mobility parameters, must be called here since
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

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  std::string name = "Masetti_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Masetti<EvalT, Traits>::
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
      ScalarT lattMob = mumax * pow(normTemp, gamma);

      // Obtain doping concentration
      const ScalarT& Na = acceptor(cell,point); // scaled conc.
      const ScalarT& Nd = donor(cell,point);
      ScalarT Ntot = (Na + Nd) * C0;         // unscale the conc.

      ScalarT term1 = mumin1*exp(-pc/Ntot);
      ScalarT term2 = (lattMob - mumin2) / (1.0+pow(Ntot/cr,alpha) );
      ScalarT term3 = mu1 / (1.0+pow(cs/Ntot,beta) );

      tmpMob(cell,point) = (term1 + term2 - term3) / Mu0;
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
void Mobility_Masetti<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for Masetti electron mobility model
  if (carrType == "Electron")
  {
    string species = "Arsenic";  // by default
    if (mobParamList.isParameter("Dopant Species"))
      species = mobParamList.get<string>("Dopant Species");

    // Retrieve parameters from charon::Material_Properties by default
    if (species == "Arsenic")
    {
      mumax = matProperty.getPropertyValue(matName, "Masetti Electron(As) mumax");
      mumin1 = matProperty.getPropertyValue(matName, "Masetti Electron(As) mumin1");
      mumin2 = matProperty.getPropertyValue(matName, "Masetti Electron(As) mumin2");
      mu1 = matProperty.getPropertyValue(matName, "Masetti Electron(As) mu1");
      gamma = matProperty.getPropertyValue(matName, "Masetti Electron(As) gamma");
      pc = matProperty.getPropertyValue(matName, "Masetti Electron(As) pc");
      cr = matProperty.getPropertyValue(matName, "Masetti Electron(As) cr");
      cs = matProperty.getPropertyValue(matName, "Masetti Electron(As) cs");
      alpha = matProperty.getPropertyValue(matName, "Masetti Electron(As) alpha");
      beta = matProperty.getPropertyValue(matName, "Masetti Electron(As) beta");
    }
    else if (species == "Phosphorous")
    {
      mumax = matProperty.getPropertyValue(matName, "Masetti Electron(P) mumax");
      mumin1 = matProperty.getPropertyValue(matName, "Masetti Electron(P) mumin1");
      mumin2 = matProperty.getPropertyValue(matName, "Masetti Electron(P) mumin2");
      mu1 = matProperty.getPropertyValue(matName, "Masetti Electron(P) mu1");
      gamma = matProperty.getPropertyValue(matName, "Masetti Electron(P) gamma");
      pc = matProperty.getPropertyValue(matName, "Masetti Electron(P) pc");
      cr = matProperty.getPropertyValue(matName, "Masetti Electron(P) cr");
      cs = matProperty.getPropertyValue(matName, "Masetti Electron(P) cs");
      alpha = matProperty.getPropertyValue(matName, "Masetti Electron(P) alpha");
      beta = matProperty.getPropertyValue(matName, "Masetti Electron(P) beta");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid Dopant Species ! Must be either Arsenic or Phosphorous !");
  }

  // Set up parameters for hole arora mobility model
  else if (carrType == "Hole")
  {
    // Retrieve parameters from charon::Material_Properties by default
    mumax = matProperty.getPropertyValue(matName, "Masetti Hole mumax");
    mumin1 = matProperty.getPropertyValue(matName, "Masetti Hole mumin1");
    mumin2 = matProperty.getPropertyValue(matName, "Masetti Hole mumin2");
    mu1 = matProperty.getPropertyValue(matName, "Masetti Hole mu1");
    gamma = matProperty.getPropertyValue(matName, "Masetti Hole gamma");
    pc = matProperty.getPropertyValue(matName, "Masetti Hole pc");
    cr = matProperty.getPropertyValue(matName, "Masetti Hole cr");
    cs = matProperty.getPropertyValue(matName, "Masetti Hole cs");
    alpha = matProperty.getPropertyValue(matName, "Masetti Hole alpha");
    beta = matProperty.getPropertyValue(matName, "Masetti Hole beta");
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite parameters when specified by users
  if (mobParamList.isParameter("mumax"))
    mumax = mobParamList.get<double>("mumax");
  if (mobParamList.isParameter("mumin1"))
    mumin1 = mobParamList.get<double>("mumin1");
  if (mobParamList.isParameter("mumin2"))
    mumin2 = mobParamList.get<double>("mumin2");
  if (mobParamList.isParameter("mu1"))
    mu1 = mobParamList.get<double>("mu1");
  if (mobParamList.isParameter("gamma"))
    gamma = mobParamList.get<double>("gamma");
  if (mobParamList.isParameter("pc"))
    pc = mobParamList.get<double>("pc");
  if (mobParamList.isParameter("cr"))
    cr = mobParamList.get<double>("cr");
  if (mobParamList.isParameter("cs"))
    cs = mobParamList.get<double>("cs");
  if (mobParamList.isParameter("alpha"))
    alpha = mobParamList.get<double>("alpha");
  if (mobParamList.isParameter("beta"))
    beta = mobParamList.get<double>("beta");

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Masetti<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Masetti", "Masetti low-field mobility model");
  p->sublist("Mobility ParameterList").set<std::string>("Dopant Species", "Arsenic", "Arsenic or Phosphorous as n-type dopant");
  p->sublist("Mobility ParameterList").set<double>("mumax", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mumin1", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mumin2", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu1", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("gamma", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("pc", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cr", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cs", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("alpha", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("beta", 0., "[1]");

  p->set("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
