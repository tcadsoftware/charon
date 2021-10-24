
#ifndef CHARON_MOBILITY_UNIBO_IMPL_HPP
#define CHARON_MOBILITY_UNIBO_IMPL_HPP

#include <cmath>
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Names.hpp"

/*
University of Bologna mobility model -- reference by Susanna Reggiani, et al.,
"Electron and Hole Mobility in Silicon at Large Operating Temperatures -
Part I: Bulk Mobility," IEEE Transactions on Electron Devices, Vol.49, pp.490-499, 2002.

This mobility model is a low-field model. It is based on the Masetti approach with
two major extensions. First, attractive and repulsive scattering are separately
accounted for, therefore, leading to a function of both donor and acceptor
concentrations. This automatically accounts for different mobility values for
majority and minority carriers and ensures continuity at the junctions as long
as the impurity concentrations are continuous functions.

Second, a suitable temperature dependence for most model parameters is introduced
to predict correctly the temperature dependence of carrier mobility in a wide
range of temperatures (300-700K). The model is also good for a wide range of
doping concentrations similar to Massetti model.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_UniBo<EvalT, Traits>::
Mobility_UniBo(
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

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  std::string name = "UniBo_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_UniBo<EvalT, Traits>::
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
      ScalarT mulatt = mumax * pow(normTemp, (gamma+c*normTemp) );
      ScalarT Cr1 = cr1*pow(normTemp,cr1_exp);
      ScalarT Cr2 = cr2*pow(normTemp,cr2_exp);
      ScalarT Cs1 = cs1*pow(normTemp,cs1_exp);
      ScalarT Cs2 = cs2;

      // Obtain doping concentration
      ScalarT Na = acceptor(cell,point)*C0;  // unscaled
      ScalarT Nd = donor(cell,point)*C0;

      ScalarT mu0 = (mu0d*pow(normTemp,mu0d_exp)*Nd + mu0a*pow(normTemp,mu0a_exp)*Na)/(Nd+Na);
      ScalarT mu1 = (mu1d*pow(normTemp,mu1d_exp)*Nd + mu1a*pow(normTemp,mu1a_exp)*Na)/(Nd+Na);

      ScalarT term1 = mu0;
      ScalarT term2 = (mulatt - mu0) / (1.0 + pow(Nd/Cr1,alpha1) + pow(Na/Cr2, alpha2) );
      ScalarT term3 = mu1 / (1.0 + pow((Nd/Cs1 + Na/Cs2),-2.0) );

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
void Mobility_UniBo<EvalT, Traits>::initMobilityParams
  (const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for UniBo electron mobility model
  if (carrType == "Electron")
  {
    string species = "Arsenic";  // by default
    if (mobParamList.isParameter("Dopant Species"))
      species = mobParamList.get<string>("Dopant Species");

    // Retrieve parameters from charon::Material_Properties by default
    if (species == "Arsenic")
    {
      mumax = matProperty.getPropertyValue(matName, "UniBo Electron(As) mumax");
      c = matProperty.getPropertyValue(matName, "UniBo Electron(As) c");
      gamma = matProperty.getPropertyValue(matName, "UniBo Electron(As) gamma");
      mu0d = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu0d");
      mu0d_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu0d_exp");
      mu0a = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu0a");
      mu0a_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu0a_exp");
      mu1d = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu1d");
      mu1d_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu1d_exp");
      mu1a = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu1a");
      mu1a_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) mu1a_exp");
      cr1 = matProperty.getPropertyValue(matName, "UniBo Electron(As) cr1");
      cr1_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) cr1_exp");
      cr2 = matProperty.getPropertyValue(matName, "UniBo Electron(As) cr2");
      cr2_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) cr2_exp");
      cs1 = matProperty.getPropertyValue(matName, "UniBo Electron(As) cs1");
      cs1_exp = matProperty.getPropertyValue(matName, "UniBo Electron(As) cs1_exp");
      cs2 = matProperty.getPropertyValue(matName, "UniBo Electron(As) cs2");
      alpha1 = matProperty.getPropertyValue(matName, "UniBo Electron(As) alpha1");
      alpha2 = matProperty.getPropertyValue(matName, "UniBo Electron(As) alpha2");
    }
    else if (species == "Phosphorous")
    {
      mumax = matProperty.getPropertyValue(matName, "UniBo Electron(P) mumax");
      c = matProperty.getPropertyValue(matName, "UniBo Electron(P) c");
      gamma = matProperty.getPropertyValue(matName, "UniBo Electron(P) gamma");
      mu0d = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu0d");
      mu0d_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu0d_exp");
      mu0a = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu0a");
      mu0a_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu0a_exp");
      mu1d = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu1d");
      mu1d_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu1d_exp");
      mu1a = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu1a");
      mu1a_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) mu1a_exp");
      cr1 = matProperty.getPropertyValue(matName, "UniBo Electron(P) cr1");
      cr1_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) cr1_exp");
      cr2 = matProperty.getPropertyValue(matName, "UniBo Electron(P) cr2");
      cr2_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) cr2_exp");
      cs1 = matProperty.getPropertyValue(matName, "UniBo Electron(P) cs1");
      cs1_exp = matProperty.getPropertyValue(matName, "UniBo Electron(P) cs1_exp");
      cs2 = matProperty.getPropertyValue(matName, "UniBo Electron(P) cs2");
      alpha1 = matProperty.getPropertyValue(matName, "UniBo Electron(P) alpha1");
      alpha2 = matProperty.getPropertyValue(matName, "UniBo Electron(P) alpha2");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid Dopant Species ! Must be either Arsenic or Phosphorous !");
  }

  // Set up parameters for hole arora mobility model
  else if (carrType == "Hole")
  {
    // Retrieve parameters from charon::Material_Properties by default
    mumax = matProperty.getPropertyValue(matName, "UniBo Hole mumax");
    c = matProperty.getPropertyValue(matName, "UniBo Hole c");
    gamma = matProperty.getPropertyValue(matName, "UniBo Hole gamma");
    mu0d = matProperty.getPropertyValue(matName, "UniBo Hole mu0d");
    mu0d_exp = matProperty.getPropertyValue(matName, "UniBo Hole mu0d_exp");
    mu0a = matProperty.getPropertyValue(matName, "UniBo Hole mu0a");
    mu0a_exp = matProperty.getPropertyValue(matName, "UniBo Hole mu0a_exp");
    mu1d = matProperty.getPropertyValue(matName, "UniBo Hole mu1d");
    mu1d_exp = matProperty.getPropertyValue(matName, "UniBo Hole mu1d_exp");
    mu1a = matProperty.getPropertyValue(matName, "UniBo Hole mu1a");
    mu1a_exp = matProperty.getPropertyValue(matName, "UniBo Hole mu1a_exp");
    cr1 = matProperty.getPropertyValue(matName, "UniBo Hole cr1");
    cr1_exp = matProperty.getPropertyValue(matName, "UniBo Hole cr1_exp");
    cr2 = matProperty.getPropertyValue(matName, "UniBo Hole cr2");
    cr2_exp = matProperty.getPropertyValue(matName, "UniBo Hole cr2_exp");
    cs1 = matProperty.getPropertyValue(matName, "UniBo Hole cs1");
    cs1_exp = matProperty.getPropertyValue(matName, "UniBo Hole cs1_exp");
    cs2 = matProperty.getPropertyValue(matName, "UniBo Hole cs2");
    alpha1 = matProperty.getPropertyValue(matName, "UniBo Hole alpha1");
    alpha2 = matProperty.getPropertyValue(matName, "UniBo Hole alpha2");
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite parameters when specified by users
  if (mobParamList.isParameter("mumax"))
    mumax = mobParamList.get<double>("mumax");
  if (mobParamList.isParameter("c"))
    c = mobParamList.get<double>("c");
  if (mobParamList.isParameter("gamma"))
    gamma = mobParamList.get<double>("gamma");
  if (mobParamList.isParameter("mu0d"))
    mu0d = mobParamList.get<double>("mu0d");
  if (mobParamList.isParameter("mu0d_exp"))
    mu0d_exp = mobParamList.get<double>("mu0d_exp");
  if (mobParamList.isParameter("mu0a"))
    mu0a = mobParamList.get<double>("mu0a");
  if (mobParamList.isParameter("mu0a_exp"))
    mu0a_exp = mobParamList.get<double>("mu0a_exp");
  if (mobParamList.isParameter("mu1d"))
    mu1d = mobParamList.get<double>("mu1d");
  if (mobParamList.isParameter("mu1d_exp"))
    mu1d_exp = mobParamList.get<double>("mu1d_exp");
  if (mobParamList.isParameter("mu1a"))
    mu1a = mobParamList.get<double>("mu1a");
  if (mobParamList.isParameter("mu1a_exp"))
    mu1a_exp = mobParamList.get<double>("mu1a_exp");
  if (mobParamList.isParameter("cr1"))
    cr1 = mobParamList.get<double>("cr1");
  if (mobParamList.isParameter("cr1_exp"))
    cr1_exp = mobParamList.get<double>("cr1_exp");
  if (mobParamList.isParameter("cr2"))
    cr2 = mobParamList.get<double>("cr2");
  if (mobParamList.isParameter("cr2_exp"))
    cr2_exp = mobParamList.get<double>("cr2_exp");
  if (mobParamList.isParameter("cs1"))
    cs1 = mobParamList.get<double>("cs1");
  if (mobParamList.isParameter("cs1_exp"))
    cs1_exp = mobParamList.get<double>("cs1_exp");
  if (mobParamList.isParameter("cs2"))
    cs2 = mobParamList.get<double>("cs2");
  if (mobParamList.isParameter("alpha1"))
    alpha1 = mobParamList.get<double>("alpha1");
  if (mobParamList.isParameter("alpha2"))
    alpha2 = mobParamList.get<double>("alpha2");

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_UniBo<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "UniBo", "UniBo low-field mobility model");
  p->sublist("Mobility ParameterList").set<std::string>("Dopant Species", "Arsenic", "Arsenic or Phosphorous as n-type dopant");
  p->sublist("Mobility ParameterList").set<double>("mumax", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("c", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("gamma", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("mu0d", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu0a", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu1d", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu1a", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu0d_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("mu0a_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("mu1d_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("mu1a_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("cr1", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cr2", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cs1", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cs2", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cr1_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("cr2_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("cs1_exp", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("alpha1", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("alpha2", 0., "[1]");

  p->set("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
