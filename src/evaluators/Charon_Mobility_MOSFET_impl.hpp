
#ifndef CHARON_MOBILITY_MOSFET_IMPL_HPP
#define CHARON_MOBILITY_MOSFET_IMPL_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_MOSFET<EvalT, Traits>::
Mobility_MOSFET(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
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

  // Initialize the MOSFET mobility parameters, must be called here since
  // it uses the above variables in the code
  initMobilityParams(matName, mobParamList);

  //create a mobility scaling factor for simple scaling or (eventually) homotopy if it's determined it's needed

  mobilityScaling = rcp(new panzer::ScalarParameterEntry<EvalT>);
  mobilityScaling->setRealValue(1.0);
  if(mobParamList.isParameter("Mobility Scaling"))
    if(mobParamList.isType<double>("Mobility Scaling"))
      {
	mobilityScaling->setRealValue(mobParamList.get<double>("Mobility Scaling"));
      }

  // Set data layouts for the fields
  RCP<DataLayout> output_scalar = ip_scalar;  // default for SUPG-FEM
  RCP<DataLayout> input_scalar = ip_scalar;
  RCP<DataLayout> input_vector = ip_vector;
  num_points = num_ips;

  // Retrieve edge data layout and cellType if isEdgedl = true (for CVFEM-SG and EFFPG-FEM)
  RCP<DataLayout> edge_scalar;
  if (isEdgedl)
  {
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    edge_scalar = cellTopoInfo->edge_scalar;
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
  X0 = scaleParams->scale_params.X0;
  E0 = scaleParams->scale_params.E0;
  T0 = scaleParams->scale_params.T0;

  // Dependent field
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,input_scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,input_scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  //Dependent mobility models
  std::string bulkName;
  if(carrType == "Electron")
    {
      if(bulkMobilityModel == "Klaassen")
	bulkName = n.field.elec_philips_thomas_mobility;
      else
	std::cout<<" FAIL THIS ELECTRON "<<std::endl;
    }
  else if(carrType == "Hole")
    {
      if(bulkMobilityModel == "Klaassen")
	bulkName = n.field.hole_philips_thomas_mobility;
      else
	std::cout<<" FAIL THIS HOLE"<<std::endl;
    }
  //There must always be a bulk mobility model specified.
  if(isEdgedl)
    bulk_mobility = MDField<ScalarT,Cell,Edge>(bulkName, edge_scalar);
  else
    bulk_mobility = MDField<ScalarT,Cell,Point>(bulkName, input_scalar);
  this->addDependentField(bulk_mobility);

  //Dependent mobility models
  std::string perpName = "";
  if(carrType == "Electron")
    {
      if(perpMobilityModel == "Shirahata")
	perpName = n.field.elec_shirahata_mobility;
      else
	{
	  perpMobilityModel = "";  //Sanity step
	}
    }
  else if(carrType == "Hole")
    {
      if(perpMobilityModel == "Shirahata")
	perpName = n.field.hole_shirahata_mobility;
      else
	{
	  perpMobilityModel = "";  //Sanity step
	}
    }
  //There doesn't need to be a perpendicular field model specified.

  if(perpMobilityModel != "")
    {
      if(isEdgedl)
	perp_mobility = MDField<ScalarT,Cell,Edge>(perpName, edge_scalar);
      else
	perp_mobility = MDField<ScalarT,Cell,Point>(perpName, input_scalar);
      this->addDependentField(perp_mobility);
    }


  if (hiFieldOn and !isEdgedl) //This adds an electric field dependence if supg--if SG or EEFPG, field has to be calculated on edge.
    {
      if (carrType == "Electron")
        eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield, input_vector);
      else if (carrType == "Hole")
        eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield, input_vector);
      this->addDependentField(eff_field);
    }

  //Create objects that compute the various parts of the mobility

  std::string name = "MOSFET_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_MOSFET<EvalT, Traits>::
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
Mobility_MOSFET<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;


  //Grab the scaling factor (defaults to 1.0)
  double mobilityScaleFactor = mobilityScaling->getRealValue();

 // ----------------------------------------------------------------------------
 // Want mobility at the primary edge center for CVFEM-SG and EFFPG-FEM
 // ----------------------------------------------------------------------------

 if (isEdgedl)
 {
   // compute low or high-field mobility at the center of primary edge
   // depending on the hiFieldOn flag from the bulk mobility used
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
   {
     for (int edge = 0; edge < num_edges; ++edge)
     {
       // Get the bulk mobility
       ScalarT edgelhfMob =  bulk_mobility(cell,edge) * Mu0;

       // primary edge perp-field mobility in [cm2/V.s]
       //ScalarT edgeperpMob = (perpFieldMob(cell,node0) + perpFieldMob(cell,node1)) / 2.0;
       ScalarT edgeperpMob = 0.0;
       if(perpMobilityModel != "")
	 edgeperpMob = perp_mobility(cell,edge) * Mu0;

       ScalarT tempMob = 1/edgelhfMob;
       if(perpMobilityModel != "" and edgeperpMob > fabs(1e-9))
	 tempMob += 1.0/edgeperpMob;

       // primary edge mobility (scaled);
       mobility(cell,edge) = mobilityScaleFactor / (tempMob) / Mu0;
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
       // compute low or high-field mobility depending on the hiFieldOn 
       // flag from the bulk mobility used
       ScalarT lhfMob =  bulk_mobility(cell,point) * Mu0;
       
       ScalarT tempipMob = 1/lhfMob;
       ScalarT edgeperpMob = 0.0;
       if(perpMobilityModel != "")
	 edgeperpMob = perp_mobility(cell,point) * Mu0;

       if(perpMobilityModel != "" and edgeperpMob > 0.0)
	 tempipMob += 1.0/(edgeperpMob);

       // mobility at IP (scaled)
       mobility(cell,point) = mobilityScaleFactor / (tempipMob) / Mu0;

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
void Mobility_MOSFET<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up default parameters for the MOSFET mobility model
  if (carrType == "Electron")
  {
    sign = 1.0;   // + for n and - for p
    lowMob = matProperty.getPropertyValue(matName, "Electron Mobility");  // default
    vsat300 = matProperty.getPropertyValue(matName, "Electron Saturation Velocity");
  }
  else if (carrType == "Hole")
  {
    sign = -1.0;   // + for n and - for p
    lowMob = matProperty.getPropertyValue(matName, "Hole Mobility");  // default
    vsat300 = matProperty.getPropertyValue(matName, "Hole Saturation Velocity");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  //The bulk mobility model
  bulkMobilityModel = mobParamList.get<std::string>("Bulk Mobility");

  //The perp mobility model
  if(mobParamList.isParameter("Perpendicular Field Model"))
    perpMobilityModel = mobParamList.get<std::string>("Perpendicular Field Model");

  // Overwrite the default parameters when given in the input xml
  if (mobParamList.isParameter("Saturation Velocity"))
    vsat300 = mobParamList.get<double>("Saturation Velocity");


}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_MOSFET<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList").set<double>("Mobility Scaling",1.0);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "MOSFET", "MOSFET mobility model");


  p->sublist("Mobility ParameterList").set<std::string>("Bulk Mobility","","Klaassen");
  p->sublist("Mobility ParameterList").set<std::string>("Perpendicular Field Model","","");

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
