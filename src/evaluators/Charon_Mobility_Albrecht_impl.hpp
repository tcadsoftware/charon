
#ifndef CHARON_MOBILITY_ALBRECHT_IMPL_HPP
#define CHARON_MOBILITY_ALBRECHT_IMPL_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"


/*
The low-field electron mobility for GaN comes from
J. Appl. Phys. 83, 4777 (1998) J.D. Albrecht et al.

For the SUPG-FEM formulation, hiField lives at IPs; while for CVFEM-SG and EFFPG-FEM,
hiField lives at centers of primary edges.

An example of using the model is given below:

<ParameterList name="Electron Mobility">
    <Parameter name="Value" type="string" value="Albrecht" />
</ParameterList>

<ParameterList name="Hole Mobility">
    <Parameter name="Value" type="string" value="Albrecht" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Albrecht<EvalT, Traits>::
Mobility_Albrecht(
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

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // BASIS
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis_scalar->dimension(1);
  basis_name = basis->name();

  // Check if edge data layout is specified; if yes, isEdgedl = true
  isEdgedl = p.get<bool>("Is Edge Data Layout");

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize the GaN mobility parameters, must be called here since
  // it uses the above variables in the code
  initMobilityParams(matName, mobParamList);

  // Set data layouts for the fields
  RCP<DataLayout> output_scalar = ip_scalar;  // default for SUPG-FEM
  RCP<DataLayout> input_scalar = ip_scalar;
  RCP<DataLayout> input_vector = ip_vector;
  num_points = num_ips;

  // Retrieve edge data layout and cellType if isEdgedl = true (for CVFEM-SG and EFFPG-FEM)
  if (isEdgedl)
  {
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
    RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
    num_edges = edge_scalar->dimension(1);
    cellType = cellTopoInfo->getCellTopology();

    // Overwrite the data layouts
    output_scalar = edge_scalar;
    input_scalar = basis_scalar;
    input_vector = edge_vector;
    num_points = num_nodes;
  }

  // Evaluated field
  if (carrType == "Electron")
    mobility = MDField<ScalarT,Cell,Point>(n.field.elec_mobility, output_scalar);
  else if (carrType == "Hole")
    mobility = MDField<ScalarT,Cell,Point>(n.field.hole_mobility, output_scalar);
  this->addEvaluatedField(mobility);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,input_scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,input_scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  if (isEdgedl)  // for CVFEM-SG and EFFPG-FEM
  {
    // Always need these fields when isEdgedl = true
    intrin_fermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi, input_scalar);
    bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap, input_scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, input_scalar);

    this->addDependentField(intrin_fermi);
    this->addDependentField(bandgap);
    this->addDependentField(eff_bandgap);

  }  // end of if (isEdgedl)

  std::string name = "Albrecht_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Albrecht<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Albrecht<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold low-field mobility
  Kokkos::DynRankView<ScalarT,PHX::Device> lowFieldMob = Kokkos::createDynRankView(mobility.get_static_view(),"lowFieldMob",workset.num_cells, num_points);

  // Compute the low-field mobility at IP or BASIS depending on isEdgedl
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Obtain temperature [K]
      ScalarT latt = latt_temp(cell,point)*T0;

      // obtain doping and carrier concentrations
      ScalarT Na = acceptor(cell,point);
      ScalarT Nd = donor(cell,point);
      //ScalarT Ntot = (Na + Nd) * C0; // unscaled conc.


      ScalarT lfMob = evaluateLowFieldMobility(Na, Nd, latt);
      lowFieldMob(cell,point) = lfMob ;  // in [cm2/V.s]
    }
  }

 // ----------------------------------------------------------------------------
 // Want mobility at the primary edge center for CVFEM-SG and EFFPG-FEM
 // ----------------------------------------------------------------------------

 if (isEdgedl)
 {
   // compute high-field mobility at the center of primary edge
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
   {
     for (int edge = 0; edge < num_edges; ++edge)
     {
       // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
       int node0 = cellType->getNodeMap(1,edge,0);
       int node1 = cellType->getNodeMap(1,edge,1);

       // primary edge low-field mobility in [cm2/V.s]
       ScalarT edgelfMob = (lowFieldMob(cell,node0) + lowFieldMob(cell,node1)) / 2.0;

       // primary edge mobility (scaled)
       mobility(cell,edge) = edgelfMob / Mu0;
     }
   }

 }  // end of if (isEdgedl)


 // ----------------------------------------------------------------------------
 // ----- Want mobility available at IP for SUPG-FEM ---------------------------
 // ----------------------------------------------------------------------------

 else
 {
   // loop over cells
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
   {
     for (int point = 0; point < num_points; ++point)
     {
       ScalarT lfMob = lowFieldMob(cell,point);

       // mobility at IP (scaled)
       mobility(cell,point) = lfMob / Mu0;
     }
   }

 }  // end of else block

}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_Albrecht<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up default parameters for the Albrecht mobility model
  if (carrType == "Electron")
  {
    lowMob = matProperty.getPropertyValue(matName, "Electron Mobility");  // default
    exa = matProperty.getPropertyValue(matName, "Albrecht AN");
    exb = matProperty.getPropertyValue(matName, "Albrecht BN");
    exc = matProperty.getPropertyValue(matName, "Albrecht CN");
    exN0 = matProperty.getPropertyValue(matName, "Albrecht N0N");
    exT0 = matProperty.getPropertyValue(matName, "Albrecht T0N");
    exT1 = matProperty.getPropertyValue(matName, "Albrecht T1N");
  }
  else if (carrType == "Hole")
  {
    lowMob = matProperty.getPropertyValue(matName, "Hole Mobility");  // default
    exa = matProperty.getPropertyValue(matName, "Albrecht AP");
    exb = matProperty.getPropertyValue(matName, "Albrecht BP");
    exc = matProperty.getPropertyValue(matName, "Albrecht CP");
    exN0 = matProperty.getPropertyValue(matName, "Albrecht N0P");
    exT0 = matProperty.getPropertyValue(matName, "Albrecht T0P");
    exT1 = matProperty.getPropertyValue(matName, "Albrecht T1P");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite the default parameters when given in the input xml
  if (mobParamList.isParameter("expa"))
    exa = mobParamList.get<double>("expa");
  if (mobParamList.isParameter("expb"))
    exb = mobParamList.get<double>("expb");
  if (mobParamList.isParameter("expc"))
    exc = mobParamList.get<double>("expc");
  if (mobParamList.isParameter("expn0"))
    exN0 = mobParamList.get<double>("expN0");
  if (mobParamList.isParameter("expT0"))
    exT0 = mobParamList.get<double>("expT0");
  if (mobParamList.isParameter("expT1"))
    exT1 = mobParamList.get<double>("expT1");

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateLowFieldMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Albrecht<EvalT,Traits>::ScalarT
Mobility_Albrecht<EvalT, Traits>::evaluateLowFieldMobility(const ScalarT& Na, const ScalarT& Nd, const ScalarT& latt)
{

  ScalarT lfMob = 1e-20;
  ScalarT Ntot = (Na + Nd) * C0;
  ScalarT lattratio = 0.0;
  ScalarT doperatio = 0.0;

  lattratio = latt / exT0;
  doperatio = Ntot / exN0;
  lfMob = exa * Ntot/exN0 * std::pow(lattratio,-1.5) *
    std::log( 1.0 + 3.0 * std::pow(lattratio,2.0) * std::pow(doperatio,(-2.0/3.0))) +
    exb * std::pow(lattratio,1.5) + (exc / (std::exp(exT1 / latt) - 1.0));

  return 1.0/lfMob;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Albrecht<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Albrecht", "Albrecht mobility model");

  // for low-field part
  //p->sublist("Mobility ParameterList").set<std::string>("Low Field Mobility Model", "?", "Albrecht");
  p->sublist("Mobility ParameterList").set<double>("expa", 0., "[(V.s)/cm^2]");
  p->sublist("Mobility ParameterList").set<double>("expb", 0., "[(V.s)/cm^2]");
  p->sublist("Mobility ParameterList").set<double>("expc", 0., "[(V.s)/cm^2]");
  p->sublist("Mobility ParameterList").set<double>("expN0", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("expT0", 0., "[K]");
  p->sublist("Mobility ParameterList").set<double>("expT1", 0., "[K]");

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
