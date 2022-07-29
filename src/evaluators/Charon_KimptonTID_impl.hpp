
#ifndef CHARON_INS_KIMPTONTID_IMPL_HPP
#define CHARON_INS_KIMPTONTID_IMPL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <typeinfo>

// Charon
#include "Charon_KimptonTID.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

// Panzer
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

// Intrepid
#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"

#include "Teuchos_DefaultComm.hpp"


namespace {
  // compute triangle area given nodes coordinates
  double compute_triangle_area(
		 std::vector<double>& coordX,
		 std::vector<double>& coordY,
		 std::vector<double>& coordZ) {
    double a = std::sqrt( 
		   (coordX[1]-coordX[0])*(coordX[1]-coordX[0]) + 
                   (coordY[1]-coordY[0])*(coordY[1]-coordY[0]) +	
		   (coordZ[1]-coordZ[0])*(coordZ[1]-coordZ[0])
		 );
    double b = std::sqrt( 
		   (coordX[2]-coordX[0])*(coordX[2]-coordX[0]) + 
                   (coordY[2]-coordY[0])*(coordY[2]-coordY[0]) +	
		   (coordZ[2]-coordZ[0])*(coordZ[2]-coordZ[0])
		 );
    double c = std::sqrt( 
		   (coordX[2]-coordX[1])*(coordX[2]-coordX[1]) + 
                   (coordY[2]-coordY[1])*(coordY[2]-coordY[1]) +	
		   (coordZ[2]-coordZ[1])*(coordZ[2]-coordZ[1])
		 );
    double s = 0.5*(a+b+c);
    return std::sqrt(s*(s-a)*(s-b)*(s-c)); // um^2
  }


  // compute quadrilateral area given nodes coordinates
  double compute_quadrilateral_area(
		 std::vector<double>& coordX,
		 std::vector<double>& coordY,
		 std::vector<double>& coordZ) {	
    // assumption is that the face nodes area ordered
    // this is not checked yet
    double a = std::sqrt( 
		   (coordX[1]-coordX[0])*(coordX[1]-coordX[0]) + 
                   (coordY[1]-coordY[0])*(coordY[1]-coordY[0]) +	
		   (coordZ[1]-coordZ[0])*(coordZ[1]-coordZ[0])
		 );
    double b = std::sqrt( 
		   (coordX[2]-coordX[1])*(coordX[2]-coordX[1]) + 
                   (coordY[2]-coordY[1])*(coordY[2]-coordY[1]) +	
		   (coordZ[2]-coordZ[1])*(coordZ[2]-coordZ[1])
		 );
    double c = std::sqrt( 
		   (coordX[3]-coordX[2])*(coordX[3]-coordX[2]) + 
                   (coordY[3]-coordY[2])*(coordY[3]-coordY[2]) +	
		   (coordZ[3]-coordZ[2])*(coordZ[3]-coordZ[2])
		 );
    double d = std::sqrt( 
		   (coordX[0]-coordX[3])*(coordX[0]-coordX[3]) + 
                   (coordY[0]-coordY[3])*(coordY[0]-coordY[3]) +	
		   (coordZ[0]-coordZ[3])*(coordZ[0]-coordZ[3])
		 );
    double s = 0.5*(a+b+c+d);
    double p = std::sqrt( 
		   (coordX[0]-coordX[2])*(coordX[0]-coordX[2]) + 
                   (coordY[0]-coordY[2])*(coordY[0]-coordY[2]) +	
		   (coordZ[0]-coordZ[2])*(coordZ[0]-coordZ[2])
		 );
    double q = std::sqrt( 
		   (coordX[1]-coordX[3])*(coordX[1]-coordX[3]) + 
                   (coordY[1]-coordY[3])*(coordY[1]-coordY[3]) +	
		   (coordZ[1]-coordZ[3])*(coordZ[1]-coordZ[3])
		 );
    return 
      std::sqrt((s-a)*(s-b)*(s-c)*(s-d)-0.25*(a*c+b*d+p*q)*(a*c+b*d-p*q)); // um^2;
  }


  // compute normal to line pointing to direction
  void compute_normal2D(const double* p1,
			const double* p2,
			const double* dir_point,
			std::vector<double>& norm) {	
    norm.push_back(0.0); norm.push_back(0.0); 
    // compute direction vector
    const double vec[2] = {dir_point[0] - p1[0], 
			   dir_point[1] - p1[1]};
    // compute normals to the line defined by point1 and point2  
    const double n1[2] = {p2[1] - p1[1], -(p2[0] - p1[0])};
    const double n2[2] = {-(p2[1] - p1[1]), p2[0] - p1[0]};
    // select normal pointing towards dir_point
    if((n1[0]*vec[0] + n1[1]*vec[1]) > 0) {
      norm[0] = n1[0]; norm[1] = n1[1];
    } else {
      norm[0] = n2[0]; norm[1] = n2[1];
    }
    // normalize normal to line
    double len = std::sqrt(norm[0]*norm[0]+norm[1]*norm[1]);
    norm[0] = norm[0]/len; norm[1] = norm[1]/len;
  }


  // compute normal to triangle pointing to direction
  void compute_normal3D(const double* p1,
		   const double* p2,
		   const double* p3,
		   const double* dir_point,
		   std::vector<double>& norm) {
    norm.push_back(0.0); norm.push_back(0.0); norm.push_back(0.0);
    // compute direction vector
    const double vec[3] = {dir_point[0] - p1[0], 
			   dir_point[1] - p1[1],
                           dir_point[2] - p2[1]};
    
    // compute normals to the triangle defined by p1, p2, p3  
    const double a[3] = {p2[0] - p1[0],
			 p2[1] - p1[1],
			 p2[2] - p1[2]};
    const double b[3] = {p3[0] - p1[0],
			 p3[1] - p1[1],
			 p3[2] - p1[2]};
    // compute cross products of a and b (a x b), a = P1P2, b = P1P3
    const double n1[3] = {a[1]*b[2] - a[2]*b[1],
			  a[2]*b[0] - a[0]*b[2],
			  a[0]*b[1] - a[1]*b[0]};
    const double n2[3] = {-a[1]*b[2] + a[2]*b[1],
			  -a[2]*b[0] + a[0]*b[2],
			  -a[0]*b[1] + a[1]*b[0]};

    // select normal pointing towards dir_point
    if((n1[0]*vec[0] + n1[1]*vec[1] + n1[2]*vec[2]) > 0) {
      norm[0] = n1[0]; norm[1] = n1[1]; norm[2] = n1[2];
    } else {
      norm[0] = n2[0]; norm[1] = n2[1]; norm[2] = n2[2];
    }
    // normalize normal to triangle
    double len = std::sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
    norm[0] = norm[0]/len; norm[1] = norm[1]/len; norm[2] = norm[2]/len;
  }
			  
}



namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
KimptonTID<EvalT, Traits>::
KimptonTID(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  comm = p.get<RCP<Teuchos::Comm<int> const> >("Comm");
 
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);
  
  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();
  num_basis = basis_scalar->dimension(1);

  isSGCVFEM = p.get<bool>("Is CVFEM");

  hcurl_basis_name = "";
  if (isSGCVFEM) {
    RCP<BasisIRLayout> hcurl_basis = p.get<RCP<BasisIRLayout> >("HCurlBasis");
    hcurl_basis_name = hcurl_basis->name();
  }

  num_edges = 0;
  if (isSGCVFEM) {
    // Edge data layout
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
    num_edges = edge_scalar->dimension(1);

    // Get reference edge length
    Intrepid2::Basis_HGRAD_LINE_C1_FEM<PHX::Device> lineBasis;
    Kokkos::DynRankView<double,PHX::Device> dofCoords("dofCoords",2,1);
    lineBasis.getDofCoords(dofCoords);
    refEdgeLen = dofCoords(1,0)-dofCoords(0,0);
  } 
  // Get the primary cell topology
  cellType = basis->getCellTopologyInfo()->getCellTopology();

  // geometry
  dev_mesh = p.get<RCP<const panzer_stk::STK_Interface>>("Mesh");
  blockID = p.get<string>("Block ID");

  const string& materialName = p.get<string>("Material Name");
  mass_dens = matProperty.getPropertyValue(materialName, "Mass Density");
  
  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  E0 = scaleParams->scale_params.E0;

  dose = 0.0;
  if(p.isParameter("Dose")) 
    dose = p.get<double>("Dose"); // [rad]
  if(dose <= 0) {
    string msg = "A strictly positive radiation dose in TID model must be specified.\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  dose = dose*1e-5; // [J/g]

  DEF = 1.0;
  if(p.isParameter("Effective Dose Enhancement Factor"))
    DEF = p.get<double>("Effective Dose Enhancement Factor");
  if(DEF <= 0) {
    string msg = "A strictly positive Effective Dose Enhancement Factor in TID model must be specified.\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  Eform = 0.0;
  if(p.isParameter("Electron-Hole Pair Formation Energy"))
    Eform = p.get<double>("Electron-Hole Pair Formation Energy"); // [eV]
  if(Eform <= 0) {
    string msg = "A strictly positive Electron-Hole Pair Formation Energy in TID model must be specified.\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  Eform = Eform*cpc.q; // [J]

  pow_dep = 0.0;
  if(p.isParameter("Electric Field Power Dependency"))
    pow_dep = p.get<double>("Electric Field Power Dependency"); 
  if(pow_dep <= 0) {
    string msg = "A strictly positive Electric Field Power Dependency in TID model must be specified.\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }
  
  // interface traps
  WithInterfaceTraps = false;
  if(p.isParameter("WithInterfaceTraps"))
    WithInterfaceTraps = p.get<bool>("WithInterfaceTraps");
  if(WithInterfaceTraps) {
    // sideset ID
    if(p.isParameter("Sideset ID")) 
      sidesetID = p.get<string>("Sideset ID");
    if(sidesetID == "") {
      std::string msg = "A Sideset ID must be specified in the TID model.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    // interface trap density
    Nti = -1.0; // [cm^-2]
    bool trapDensitySweep = false;
    if (p.isParameter("Interface Trap Total Density Sweep is On"))
      trapDensitySweep = p.get<bool>("Interface Trap Total Density Sweep is On");

    if (not trapDensitySweep)
      {
	if(p.isParameter("Interface Trap Density"))
	  Nti = p.get<double>("Interface Trap Density");
	if(Nti < 0.0) {
	  std::string msg = "A positive Interface Trap Density in TID model must be specified.\n";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
	}
      }
    // interface initial filling factor
    fill_facti = -1.0;
    if(p.isParameter("Interface Initial Filling Factor")) 
      fill_facti = p.get<double>("Interface Initial Filling Factor");
    if(fill_facti < 0.0) {
      std::string msg = "A positive Interface Initial Filling Factor in TID model must be specified.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    // interface trap capture cross section
    sigma_i = 0.0; // [cm^2]
    bool crossSectionSweep = false;
    if (p.isParameter("Interface Trap Capture Cross Section Sweep is On"))
      crossSectionSweep = p.get<bool>("Interface Trap Capture Cross Section Sweep is On");

    if (not crossSectionSweep)
      {
	if(p.isParameter("Interface Trap Capture Cross Section"))
	  sigma_i = p.get<double>("Interface Trap Capture Cross Section");
	if(sigma_i <= 0.0) {
	  std::string msg = 
	    "A strictly positive Interface Trap Capture Cross Section in TID model must be specified.\n";
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
	}
      }
  }

  // volume traps
  WithVolumeTraps = false;
  ins_vol = 0.0;
  if(p.isParameter("WithVolumeTraps"))
    WithVolumeTraps = p.get<bool>("WithVolumeTraps");
  if(WithVolumeTraps) {
    // volume trap density
    Ntv = -1.0; // [cm^-3]
    if(p.isParameter("Volume Trap Density"))
      Ntv = p.get<double>("Volume Trap Density");
    if(Ntv < 0.0) {
      std::string msg = "A positive Volume Trap Density in TID model must be specified.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    // volume initial filling factor
    fill_factv = -1.0;
    if(p.isParameter("Volume Initial Filling Factor")) 
      fill_factv = p.get<double>("Volume Initial Filling Factor");
    if(fill_factv < 0.0) {
      std::string msg = "A positive Volume Initial Filling Factor in TID model must be specified.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    // volume trap capture cross section
    sigma_v = 0.0; // [cm^2]
    if(p.isParameter("Volume Trap Capture Cross Section"))
      sigma_v = p.get<double>("Volume Trap Capture Cross Section");
    if(sigma_v <= 0.0) {
      std::string msg = 
	"A strictly positive Volume Trap Capture Cross Section in TID model must be specified.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
    // trap capture volume (at E = 1MV/cm)
    trap_vol = 4./3. * cpc.pi * std::pow(sigma_v/cpc.pi,1.5); // [cm^-3]

    // volume trap critical capture cross section
    sigma_v_crit = 0.0; // [cm^2]
    if(p.isParameter("Volume Trap Critical Capture Cross Section"))
      sigma_v_crit = p.get<double>("Volume Trap Critical Capture Cross Section");
    if(sigma_v_crit <= 0.0) {
      std::string msg = 
	"A strictly positive Volume Trap Critical Capture Cross Section in TID model must be specified.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    } 
    trap_crit_vol = 4./3. * cpc.pi * std::pow(sigma_v_crit/cpc.pi,1.5);
  }

  withVaryingV = false;
  V_freeze = 0.0;
  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  if(p.isParameter("ParamLib")) {
    if (p.isParameter("Interface Trap Capture Cross Section Sweep is On"))
      {
	if (p.get<bool>("Interface Trap Capture Cross Section Sweep is On"))
	  {
	    user_value = panzer::createAndRegisterScalarParameter<EvalT>(std::string("Interface Trap Capture Cross Section Sweep"),
									 *p.get<RCP<panzer::ParamLib> >("ParamLib"));
	    interfaceTrapSweep = p.get<bool>("Interface Trap Capture Cross Section Sweep is On");
	    initialInterfaceTrapCrossSection = p.get<double>("Interface Trap Initial Capture Cross Section");
	    finalInterfaceTrapCrossSection = p.get<double>("Interface Trap Final Capture Cross Section");
	    user_value->setRealValue(0.0);
	  }
      }


    if (p.isParameter("Interface Trap Total Density Sweep is On"))
      {
	if (p.get<bool>("Interface Trap Total Density Sweep is On"))
	  {
	    user_value = panzer::createAndRegisterScalarParameter<EvalT>(std::string("Interface Trap Total Density Sweep"),
									 *p.get<RCP<panzer::ParamLib> >("ParamLib"));
	    interfaceDensitySweep = p.get<bool>("Interface Trap Total Density Sweep is On");
	    initialInterfaceTotalDensity = p.get<double>("Interface Trap Initial Total Density");
	    finalInterfaceTotalDensity = p.get<double>("Interface Trap Final Total Density");
	    user_value->setRealValue(0.0);
	  }
      }


    if (p.isParameter("Voltage Sweep"))
      {
	if (p.get<bool>("Voltage Sweep"))
	  {
	    withVaryingV = true;
	    user_value = panzer::createAndRegisterScalarParameter<EvalT>(std::string("Varying Voltage"),
                                          *p.get<RCP<panzer::ParamLib> >("ParamLib"));
	    V_freeze = p.get<double>("Freeze Voltage");
	  }
      }
  }

  // factor pre-computation
  gen_const = (mass_dens * dose * DEF) / Eform; // [cm-3]

  // input fields
  if (!isSGCVFEM) {
    grad_pot = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi,ip_vector);
    this->addDependentField(grad_pot);
  } else {
    potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,basis_scalar);
    this->addDependentField(potential);
  }

  // Evaluated fields
  ins_genpair_dens = MDField<ScalarT,Cell,Point>(n.field.ins_genpair_density,ip_scalar);
  ins_trappedcharge = MDField<ScalarT,Cell,Point>(n.field.ins_htrappedcharge,ip_scalar);
  this->addEvaluatedField(ins_genpair_dens);
  this->addEvaluatedField(ins_trappedcharge);
 
  std::string name = "Ins_TrappedCharge";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
KimptonTID<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

  if (isSGCVFEM) 
    hcurl_basis_index = panzer::getBasisIndex(hcurl_basis_name,(*sd.worksets_)[0]);

  // compute geo/top info 
  comp_geo_info(dev_mesh);

  /*
  if (WithVolumeTraps) {
    // compute insulator volume
    double local_vol = 0.0;
    for (const auto& wkst : (*sd.worksets_)) {
      double wrk_vol = 0.0;
      for (panzer::index_t cell = 0; cell < wkst.num_cells; ++cell) {
	// compute workset volume
	const panzer::index_t num_ip = 
	  (wkst.int_rules[int_rule_index])->ip_coordinates.extent_int(1);
	for (panzer::index_t ip = 0; ip < num_ip; ++ip) {
	  wrk_vol += 
	    (wkst.int_rules[int_rule_index])->weighted_measure(cell, ip); 
	}
      } 
      if (num_dims == 2) wrk_vol *= 1e4; // [um^3]
      //std::cout << "wkst.num_cells = " << wkst.num_cells << std::endl;
      //std::cout << "wrk_volume = " << wrk_vol << std::endl;
      local_vol += wrk_vol;
    }
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), 
			 &local_vol, &ins_vol);    
  }
  */

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
KimptonTID<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  using std::vector;

  const vector<stk::mesh::Entity>& localElements = 
    *dev_mesh->getElementsOrderedByLID();
  const vector<std::size_t>& localCellIds = workset.cell_local_ids;
  const panzer::IntegrationValues2<double> &irv = *workset.int_rules[int_rule_index];
  
  ScalarT scale_const = 1/C0;
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();

  ScalarT cnt_voltage = 0.0;
  if (withVaryingV) {
    cnt_voltage = user_value->getValue();
  }

  // for a LOCA sweep freeze the charge for gate voltages 
  // strictly greater than V_freeze
  if (withVaryingV and (cnt_voltage > V_freeze)) {
    return;
  }

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    std::size_t cellLocalId = localCellIds[cell];
    stk::mesh::Entity elem = localElements[cellLocalId];
    bool isOnInterface = (area_map[elem] > 0.0) ? true : false;

    // pre-compute cell quantities
    ScalarT Ntib = 0.0; // converted interface trap density
    vector<double> inter_normal; // interface normal for interface elems
    if(isOnInterface and WithInterfaceTraps) {
      // convert trap surface density to volume density
      // Ntotal = Nti*area = Ntb*vol => Ntb = (area/vol)*Nti
      // compute cell volume
      ScalarT cell_vol = 0.0; // [um^3]
      const size_t num_ip = irv.ip_coordinates.extent_int(1);
      for (size_t ip = 0; ip < num_ip; ++ip) {
	cell_vol += irv.weighted_measure(cell, ip); 
      }
      if (num_dims == 2) cell_vol *= 1e4; // [um^3]
      // compute interface face area
      ScalarT cell_inter_area = area_map[elem]; // [um^2]
      ScalarT vol_fac = cell_inter_area/cell_vol; // [um^-1] 
      vol_fac *= 1e4; // [cm^-1]
      if (interfaceDensitySweep)
	{
	  double fraction = user_value->getRealValue();
	  Nti = initialInterfaceTotalDensity*(1-fraction) + finalInterfaceTotalDensity*fraction;
	}
      Ntib = Nti;
      Ntib *= vol_fac; // [cm^-3]
      // interface normal
      for (int dim = 0; dim < num_dims; ++dim)
	inter_normal.push_back((normal_map[elem])[dim]);
    }
    
    vector<ScalarT> E_vec_centroid;
    if (isSGCVFEM) {
      // compute electric field at centroids
      computeCentroidField(workset,cell,E_vec_centroid);
    }

    // trapped charge at IPs
    for (int ip = 0; ip < num_ips; ++ip)
    {    
      ins_genpair_dens(cell,ip) = 0.0;
      ins_trappedcharge(cell,ip) = 0.0;
      
      // compute concentration of generated electrons-hole pairs
      //  gen_const = (mass_dens * dose * DEF) / Eform
      //  ins_genpair_dens(cell,ip) = (mass_dens * dose * DEF) / Eform
      ins_genpair_dens(cell,ip) = gen_const; // [cm^-3]
      //ins_genpair_dens(cell,ip) *= scale_const;

      // compute electric field magnitude 
      vector<ScalarT> E_vec;
      E_vec.resize(num_dims);
      if (!isSGCVFEM) {
	for (int dim = 0; dim < num_dims; ++dim) 
	  E_vec[dim] = -grad_pot(cell,ip,dim);
      } else {
	// SGCVFEM
	for (int dim = 0; dim < num_dims; ++dim)
	  E_vec[dim] = E_vec_centroid[ip*num_dims + dim];
      }
      ScalarT Emag = 0.0;
      for (int dim = 0; dim < num_dims; ++dim) 
	Emag += E_vec[dim]*E_vec[dim];
      if (Emag > 0.0) Emag = std::sqrt(Emag);
      Emag = Emag*E0/1e6; // [MV/cm]
      
      // compute fractional hole yield
      ScalarT fox = 0.0;
      if (Emag >= 1e-3) 
	fox = std::pow((1.35/Emag) + 1.0, -0.9);
      else
	fox = std::pow((1.35/0.001) + 1.0, -0.9);
       
      // compute effective concentration of generated electrons-hole pairs
      ins_genpair_dens(cell,ip) *= fox; // [cm^-3]

      // interface trap contribution
      if (isOnInterface and WithInterfaceTraps) { 
	// check if field is positive; we project the IP field on the 
	// local normal to the semic/ins interface face
	ScalarT field_dir = 0.0;
	for (int dim = 0; dim < num_dims; ++dim) 
	  field_dir += E_vec[dim] * inter_normal[dim];
	bool isFieldPositive = field_dir < 0.0 ? true : false;

	if (isFieldPositive) {
	  // interface trap TID is defined 
	  // only for positive fields
	  if (interfaceTrapSweep)
	    {
	      double fraction = user_value->getRealValue();
	      sigma_i = initialInterfaceTrapCrossSection*(1-fraction) + finalInterfaceTrapCrossSection*fraction;
	    }
	  ScalarT sigma_eff = sigma_i * std::pow(Emag,-pow_dep);

	  // correct for electron-hole pairs being generated 
	  // within the capturing volume of volume traps 
	  ScalarT gen_val = gen_const;
	  if(WithVolumeTraps) {
	    // probability that an electron-hole pair is generated
            // within the capturing volume of a volume trap within 
            // this particular element(cell) volume located on interface 
            // is elem_vol*Ntv*eff_trap_vol/elem_vol = Ntv*eff_trap_vol
	    ScalarT eff_trap_vol = trap_vol * std::pow(Emag,-1.5*pow_dep); // [cm^-3]
	    if (eff_trap_vol > 0.5/Ntv) eff_trap_vol = 0.5/Ntv;
	    ScalarT Pgen_tvol = Ntv * eff_trap_vol;
	    // probability that an electron-hole pair is generated
            // outside the capturing volume of volume traps
	    ScalarT Pgen = 1.0 - Pgen_tvol;
	    // compute correction of the generated electron-hole pair 
            // density due to volume trapping 
	    gen_val *= Pgen;	     
	  }
   
          // transform the capture cross section to a equivalent capturing
          // volume as the interface traps have been converted to volume
          // traps
	  ScalarT vol_capt_eff = 4./3. * cpc.pi * std::pow(sigma_eff/cpc.pi,1.5);
	  // density of new holes trapped
	  ScalarT gen_pairs = gen_val * fox;
	  ScalarT Nttp = 
           (Ntib - fill_facti*Ntib) * vol_capt_eff * gen_pairs;

	  // detrapping; primary source of energy for detrapping is 
          // the energy released by the recombination of electron-hole
          // pairs during irradiation
	  ScalarT gen_pairs_recomb = gen_val * (1.0 - fox);
	  ScalarT Ntde = fill_facti*Ntib * vol_capt_eff * gen_pairs_recomb;

          // net density of interface trapped holes
	  ins_trappedcharge(cell,ip) = 
	    (fill_facti*Ntib +  Nttp - Ntde) * scale_const;
	} // positive field
      } // interface trap contribution


      // volume trap contribution
      if (WithVolumeTraps) {
	// compute field-dependent trap capture volume
        // trap_vol is trap capture volume at Emag = 1MV/cm
	ScalarT eff_trap_vol = trap_vol * std::pow(Emag,-1.5*pow_dep); // [cm^-3]
	// for Emag -> 0, eff_trap_vol would become infinite; 
        // physically the total trap capture volume V_capt = ins_vol*Ntv*eff_trap_vol
	// is limited to ins_val/2 since at zero field the holes have a equal
        // chance of diffusing to gate or interface
	if (eff_trap_vol > 0.5/Ntv) eff_trap_vol = 0.5/Ntv;

	// volume hole capture probability; the assumption is that
        // electron-hole pair is generated within the capture volume of the trap 
	ScalarT Pcapt = trap_crit_vol / (trap_crit_vol + eff_trap_vol); 

        // trapping/detrapping
	ScalarT Nttp = Pcapt * ins_genpair_dens(cell,ip) *
	  (Ntv - fill_factv*Ntv) * eff_trap_vol;
	ScalarT Ntde = (1-Pcapt) * ins_genpair_dens(cell,ip) *
	  fill_factv*Ntv * eff_trap_vol;

	// net density of volume trapped holes
	ins_trappedcharge(cell,ip) += 
	  (fill_factv*Ntv +  Nttp - Ntde) * scale_const;
      }

    } // IPs
  } // cells

}


template<typename EvalT, typename Traits>
void
KimptonTID<EvalT, Traits>::comp_geo_info(const Teuchos::RCP<const panzer_stk::STK_Interface> mesh)
{
  using std::vector;
  using Teuchos::RCP;

  std::vector<stk::mesh::Entity> ins_elems;

  RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();
  vector<stk::mesh::Entity> side_entities;
  mesh->getAllSides(sidesetID,blockID,side_entities);
  vector<stk::mesh::Entity> inter_entities;
  for(std::size_t et = 0; et < side_entities.size(); et++) { // 2D or 3D
    if( (mesh->getDimension() == 2 and bulkData->entity_rank(side_entities[et]) == 
         stk::topology::EDGE_RANK) or (dev_mesh->getDimension() == 3 and
	 bulkData->entity_rank(side_entities[et]) == stk::topology::FACE_RANK) ) {
      inter_entities.push_back(side_entities[et]);
    }
  }

  mesh->getMyElements(blockID,ins_elems);
  vector<double> areas; 
  for(std::size_t el = 0; el < ins_elems.size(); el++) areas.push_back(0.0);
  for(std::size_t el = 0; el < ins_elems.size(); el++) {
    if(mesh->getDimension() == 2) {
      stk::mesh::Entity const *edge_i = bulkData->begin_edges(ins_elems[el]);
      stk::mesh::Entity const *edge_e = bulkData->end_edges(ins_elems[el]);
      for(; edge_i != edge_e; ++edge_i) {
	// find if any element edge is in the sideset
	vector<stk::mesh::Entity>::iterator it = 
	  std::find(inter_entities.begin(),inter_entities.end(),*edge_i);
	if(it != inter_entities.end()) {
	  //compute edge length / area
	  vector<double> coordX;
	  vector<double> coordY;
	  stk::mesh::Entity const *node_i = bulkData->begin_nodes(*edge_i);
	  stk::mesh::Entity const *node_e = bulkData->end_nodes(*edge_i);
	  vector<stk::mesh::Entity> edge_nodes; 
	  for(; node_i != node_e; ++node_i) {
	    const double * coords = mesh->getNodeCoordinates(*node_i);
	    coordX.push_back(coords[0]);
	    coordY.push_back(coords[1]);
	    edge_nodes.push_back(*node_i);
	  }
	  double l = std::sqrt( (coordX[1]-coordX[0])*(coordX[1]-coordX[0]) + 
				(coordY[1]-coordY[0])*(coordY[1]-coordY[0]));  
	  areas[el] = l*1e4; // um^2 

	  // compute normal to the edge pointing inside the ins element
	  stk::mesh::Entity const *node1_i = bulkData->begin_nodes(ins_elems[el]);
	  stk::mesh::Entity const *node1_e = bulkData->end_nodes(ins_elems[el]);
	  // find extra node in element inside the ins element
	  stk::mesh::Entity extra_node;
	  for(; node1_i != node1_e; ++node1_i) {
	    vector<stk::mesh::Entity>::iterator it1 = 
	      std::find(edge_nodes.begin(),edge_nodes.end(),*node1_i);
	    if(it1 == edge_nodes.end()) {
	      extra_node = *node1_i;
	      break;
	    } 
	  }
	  vector<double> norm; 
	  compute_normal2D(mesh->getNodeCoordinates(edge_nodes[0]), 
			   mesh->getNodeCoordinates(edge_nodes[1]), 
			   mesh->getNodeCoordinates(extra_node),
			   norm);
	  
	  normal_map.insert({ins_elems[el],norm});  
	  break; 
	} // element edge in sideset	
      } // element edge
    } else if(mesh->getDimension() == 3) {
      stk::mesh::Entity const *face_i = bulkData->begin_faces(ins_elems[el]);
      stk::mesh::Entity const *face_e = bulkData->end_faces(ins_elems[el]);
      for(; face_i != face_e; ++face_i) {
	// find if any element face in sideset
	vector<stk::mesh::Entity>::iterator it = 
	  std::find(inter_entities.begin(),inter_entities.end(),*face_i);
	if(it != inter_entities.end()) {
	  // compute face area
	  vector<double> coordX;
	  vector<double> coordY;
	  vector<double> coordZ;
	  stk::mesh::Entity const *node_i = bulkData->begin_nodes(*face_i);
	  stk::mesh::Entity const *node_e = bulkData->end_nodes(*face_i);
	  vector<stk::mesh::Entity> face_nodes;
	  for(; node_i != node_e; ++node_i) {
	    const double * coords = mesh->getNodeCoordinates(*node_i);
	    coordX.push_back(coords[0]);
	    coordY.push_back(coords[1]);
	    coordZ.push_back(coords[2]);
	    face_nodes.push_back(*node_i);
	  }
	  if(coordX.size() == 3) { // triangular face
	    areas[el] = compute_triangle_area(coordX,coordY,coordZ); // um^2
	  } else { // quadrilateral face
            // assumption is that the face nodes area ordered
            // this is not checked yet
	    areas[el] = compute_quadrilateral_area(coordX,coordY,coordZ); // um^2
	  }

	  // compute normal to the face pointing inside the ins element
	  stk::mesh::Entity const *node1_i = bulkData->begin_nodes(ins_elems[el]);
	  stk::mesh::Entity const *node1_e = bulkData->end_nodes(ins_elems[el]);
	  // find extra node in element inside the ins element
	  stk::mesh::Entity extra_node;
	  for(; node1_i != node1_e; ++node1_i) {
	    vector<stk::mesh::Entity>::iterator it1 = 
	      std::find(face_nodes.begin(),face_nodes.end(),*node1_i);
	    if(it1 == face_nodes.end()) {
	      extra_node = *node1_i;
	      break;
	    } 
	  }
	  vector<double> norm; 
	  compute_normal3D(mesh->getNodeCoordinates(face_nodes[0]), 
			   mesh->getNodeCoordinates(face_nodes[1]), 
			   mesh->getNodeCoordinates(face_nodes[2]), 
			   mesh->getNodeCoordinates(extra_node),
			   norm);  
	  normal_map.insert({ins_elems[el],norm});  
	  break;
	} // element face in sideset	
      } // element face
    } // 3D
  }

  // create an element to area map for fast processing
  for(std::size_t k = 0; k < ins_elems.size(); k++) 
    area_map.insert({ins_elems[k],areas[k]});  
}



template<typename EvalT, typename Traits>
void
KimptonTID<EvalT, Traits>::computeCentroidField(
   typename Traits::EvalData& workset,
   panzer::index_t cell, 
   std::vector<ScalarT>& E)
{
  E.resize(num_ips*num_dims);
  // zero out vector at subcontrol volume centroids
  for(int node = 0; node < num_ips; ++node) {
    for(int dim = 0; dim < num_dims; ++dim)
      E[node*num_dims + dim] = 0.0;
  }
  // loop over primary edges
  for(int iedge = 0; iedge < num_edges; ++iedge) {
    // get local node ids: first index 1 for edge
    // (0 for vertex, 2 for face, 3 for volume)
    const int node0 = cellType->getNodeMap(1,iedge,0);
    const int node1 = cellType->getNodeMap(1,iedge,1);
    
    // get local node coordinates
    double x0 = workset.cell_vertex_coordinates(cell,node0,0);
    double x1 = workset.cell_vertex_coordinates(cell,node1,0);
    double y0 = 0.0, y1 = 0.0;
    double z0 = 0.0, z1 = 0.0;
    if(num_dims > 1)  { // 2D or 3D
      y0 = workset.cell_vertex_coordinates(cell,node0,1);
      y1 = workset.cell_vertex_coordinates(cell,node1,1);
    }
    if (num_dims > 2) { // 3D
      z0 = workset.cell_vertex_coordinates(cell,node0,2);
      z1 = workset.cell_vertex_coordinates(cell,node1,2);
    }
    // compute the primary cell edge length
    double edgeLen =
      std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));

    const ScalarT& pot0 = potential(cell,node0);
    const ScalarT& pot1 = potential(cell,node1);
    ScalarT grad_pot_edge = (pot1 - pot0)/edgeLen;
    ScalarT val_edge = -grad_pot_edge;

    // evaluate driving force at the subcv centroids
    // note: number of subcv centroids is equal to the number of primary
    // nodes for quad, tri, hex, and tet mesh elements.
    //
    // In this stabilized formulation edge values are mapped to the interior
    // of the element using HCurl basis functions. In order for the values
    // to scale properly with the definitions of the HCurl and HGrad basis
    // functions in Intrepid2 as of 11/2020 we must divide by the reference
    // edge length. 
    for(int ip = 0; ip < num_ips; ++ip) {
      for(int dim = 0; dim < num_dims; ++dim) {
	E[ip*num_dims + dim] += val_edge
	  * (workset.bases[hcurl_basis_index])->basis_vector(cell,iedge,ip,dim)*edgeLen/refEdgeLen;
      }
    }
  } // end of loop over primary edges

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
KimptonTID<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::BasisIRLayout> hcurl_basis;
  p->set("HCurlBasis", basis);

  p->set<std::string>("Material Name","?");

  Teuchos::RCP<const panzer_stk::STK_Interface> mesh;
  p->set("Mesh", mesh);

  p->set<std::string>("Block ID","?");

  p->set<double>("Dose",0.0);
  
  p->set<bool>("Is CVFEM", false);

  p->set<double>("Effective Dose Enhancement Factor",1.0);
  p->set<double>("Electron-Hole Pair Formation Energy",0.0);
  p->set<double>("Electric Field Power Dependency",0.0);
  
  p->set<bool>("WithInterfaceTraps", false);
  p->set<std::string>("Sideset ID","?");
  p->set<double>("Interface Trap Density",0.0);
  p->set<double>("Interface Initial Filling Factor",-1.0);
  p->set<double>("Interface Trap Capture Cross Section",0.0);

  p->set<bool>("WithVolumeTraps", false);
  p->set<double>("Volume Trap Density",0.0);
  p->set<double>("Volume Initial Filling Factor",-1.0);
  p->set<double>("Volume Trap Capture Cross Section",0.0);
  p->set<double>("Volume Trap Critical Capture Cross Section",0.0);

  p->set<bool>("Interface Trap Capture Cross Section Sweep is On",false);
  p->set<double>("Interface Trap Initial Capture Cross Section",0.0);
  p->set<double>("Interface Trap Final Capture Cross Section",0.0);
  p->set<std::string>("Interface Trap Capture Cross Section Sweep","Parameter");

  p->set<bool>("Interface Trap Total Density Sweep is On",false);
  p->set<double>("Interface Trap Initial Total Density",0.0);
  p->set<double>("Interface Trap Final Total Density",0.0);
  p->set<std::string>("Interface Trap Total Density Sweep","Parameter");

  p->set<Teuchos::RCP<Teuchos::Comm<int> const> >("Comm",
       Teuchos::DefaultComm<int>::getComm());
 
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));
  p->set<double>("Freeze Voltage",0.0);
  
  p->set<bool>("Voltage Sweep",false);

  return p;
}


}
#endif //CHARON_KIMPTONTID_IMPL_HPP
