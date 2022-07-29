

#ifndef CHARON_BC_STRATEGY_NEUMANN_DYNAMICTRAPS_IMPL_HPP
#define CHARON_BC_STRATEGY_NEUMANN_DYNAMICTRAPS_IMPL_HPP
 
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_NeumannBC_DynamicTraps.hpp"
#include "Charon_Integrator_HJFluxDotNorm.hpp"
#include <vector>


const int MAX_NUM_TRAPS = 50;


///////////////////////////////////////////////////////////////////////////////
//
//  BCStrategy_Neumann_DynamicTraps()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
BCStrategy_Neumann_DynamicTraps(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Neumann Dynamic Traps" );
  std::cout << "Warning: Dynamic Traps do NOT work for a heterojunction!" << std::endl; 
}


///////////////////////////////////////////////////////////////////////////////
//
//  setup()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::string;
  using std::vector;
  using std::pair;

  // get the physics block parameter list
  auto pbParamList = side_pb.getParameterList();
  
  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  
  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // get the object of charon::Names
  auto names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  // check if Equation Set Name = ALL_DOFS
  string eqsn = this->m_bc.equationSetName();
  if (eqsn != "ALL_DOFS")
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Equation Set Name must be ALL_DOFS"
        << ". But you entered \"" << eqsn << "\" ! \n");

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // validate the Data parameter list
  RCP<ParameterList> valid_params = this->getValidParameters();
  dataPList->validateParameters(*valid_params);

  initDynamicTrapsParams(dataPList);
 
  // For now, assume that only one IR is used for everything.
  const auto& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);  // Doesn't support mixed ir yet!

  // obtain the integration order
  const int integration_order = ir.begin()->second->order();

  // get all dofs 
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();
  
  string dof_name = ""; 
  string residual_name = "";
  string flux_name = ""; 
  
  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
       dofs.begin(); dof_it != dofs.end(); ++dof_it)
  {
    dof_name = dof_it->first; 
    if (dof_name == names->dof.phi) // for ELECTRIC_POTENTIAL 
    {
      residual_name = "Residual_" + dof_name; 
      flux_name = "DynTraps_Charge"; 
      fluxDynTrapsCharge = flux_name; 
      this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
    }
    
    // for ELECTRON_DENSITY and HOLE_DENSITY
    if ((dof_name == names->dof.edensity) || (dof_name == names->dof.hdensity)) 
    {
      residual_name = "Residual_" + dof_name; 
      if (dof_name == names->dof.edensity) {
	flux_name = "DynTraps_eRecombination";
	eFluxDynTrapsRecomb = flux_name; 
      } else {
	flux_name = "DynTraps_hRecombination";
	hFluxDynTrapsRecomb = flux_name; 
      } 
      this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
    }
  }  // end of loop over dofs
  
  
}


///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                           const Teuchos::ParameterList& models,
                           const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector; 

  vector<string> vec_dof_name;
  vector<string> vec_flux_name;
  vector<string> vec_residual_name; 

  // get the element block id and physics id for the given physics block
  string ebID_pb = pb.elementBlockID();
  string pbID = pb.physicsBlockID();

  // get the element block ID for the given boundary condition
  string ebID_bc = this->m_bc.elementBlockID();

  if (ebID_pb != ebID_bc)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: " << pbID << " corresponds to "
      << ebID_pb << ", while the BC corresponds to " << ebID_bc << "! \n");

  const auto data = this->getResidualContributionData();

  // Iterate over each residual contribution
  for (auto eq = data.begin(); eq != data.end(); ++eq) 
  {
    const string& residual_name = std::get<0>(*eq);
    const string& dof_name = std::get<1>(*eq);
    const string& flux_name = std::get<2>(*eq);
    vec_residual_name.push_back(residual_name);
    vec_dof_name.push_back(dof_name);
    vec_flux_name.push_back(flux_name); 
  }

  // All DOFs use the same basis and ir, so can use data[0]
  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(std::get<1>(data[0]));

  // get sgcvfem boundary integration rule and basis
  std::string cvfem_type = "boundary";
  auto cvfem_bc_ir = Teuchos::rcp(new panzer::IntegrationRule(pb.cellData(),cvfem_type));
  auto Hgradbasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),ir->topology));
  auto hgrad_bc_cvfem = Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_bc_ir));

  // build and register all closure models
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);
  pb.buildAndRegisterDOFProjectionsToIPEvaluators(fm,Teuchos::null,user_data);

  // get the physics block parameter list
  auto pbParamList = pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";
  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  
  // get the equation set type
  string eqnType = eqSetPList.get<string>("Type"); 
  std::size_t foundSG = eqnType.find("SGCVFEM");
  
  // determine if the SG discretization is specified
  bool isSG = false;  
  if (foundSG != string::npos)  
    isSG = true;  // find SGCVFEM

  // get scaling parameters
  auto scaleParams = 
    user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // compute the scaled flux at the IPs of an interface
  {
    ParameterList p("Dynamic Traps Recombination");
    p.set<RCP<const charon::Names>>("Names", names);
 
    Teuchos::RCP<std::map<string,string>> block_mat =
      user_data.get<Teuchos::RCP<std::map<string,string>>>("block2mat");
    const string matName = (*block_mat)[ebID_pb];
    p.set<string>("Material Name", matName);

    p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
  
    if (isSG)  // SG
    {
      p.set("IR", cvfem_bc_ir);
      p.set("Basis", hgrad_bc_cvfem);
    }  
    else  // non-SG, i.e., FEM or EFFPG
    {
      p.set("IR", ir);
      p.set("Basis", basis);
    }  
    
    p.set<string>("Flux Dynamic Traps Charge", fluxDynTrapsCharge);        // output flux name 
    p.set<string>("Flux Dynamic Traps eRecombination", eFluxDynTrapsRecomb); // output flux name 
    p.set<string>("Flux Dynamic Traps hRecombination", hFluxDynTrapsRecomb); // output flux name
    p.set<RCP<ParameterList>>("Dynamic Traps ParameterList", dynTrapsPList); 

    if (withField) {
      // compute normals to Dynamic Traps interface
      string normal_name = pb.cellData().side() + "_DynaTrapsNorm";
      ParameterList p1(normal_name);
      p1.set("Name", normal_name);
      p1.set("Side ID", pb.cellData().side());
      if (isSG)
	p1.set("IR", cvfem_bc_ir);
      else 
	p1.set("IR", ir);
      p1.set("Normalize", true);
      const RCP< PHX::Evaluator<panzer::Traits> >
	op1 = rcp(new panzer::Normals<EvalT,panzer::Traits>(p1));
      this->template registerEvaluator<EvalT>(fm, op1);

      // compute normal dot the corresponding potential gradient
      string normal_dot_grad_name = "DynaTrapsNormDotGrad";
      ParameterList p2(normal_dot_grad_name);
      p2.set("Result Name", normal_dot_grad_name);
      p2.set("Vector A Name", names->grad_dof.phi);
      p2.set("Vector B Name", normal_name);
      p2.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
      const RCP< PHX::Evaluator<panzer::Traits> >
	op2 = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p2));
      this->template registerEvaluator<EvalT>(fm, op2);

      // pass electric field name to the NeumannBC_DynamicTraps evaluator
      p.set("EdotNorm", normal_dot_grad_name);
    }

    auto op = rcp(new charon::NeumannBC_DynamicTraps<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);

  }

  // loop over the dof names
  for (std::size_t id = 0; id < vec_dof_name.size(); ++id)
  {
    const string& residual_name = vec_residual_name[id];
    const string& flux_name = vec_flux_name[id]; 
    const string& dof_name = vec_dof_name[id];

    double multiplier = 1.0; 
    if (dof_name == names->dof.phi)   multiplier = -1.0;
    
    if (isSG)  // SG
    {
      ParameterList p(residual_name);
      p.set("Residual Name", residual_name);
      p.set("Flux Name", flux_name);
      p.set("Basis", hgrad_bc_cvfem);
      p.set("IR", cvfem_bc_ir);
      p.set("Multiplier", multiplier);
    
      const RCP< PHX::Evaluator<panzer::Traits> >
      op = rcp(new charon::Integrator_HJFluxDotNorm<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
    else  // FEM
    {
      using panzer::EvaluatorStyle;
      using panzer::Integrator_BasisTimesScalar;
      using panzer::Traits;
      using PHX::Evaluator;
    
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
        residual_name, flux_name, *basis, *ir, multiplier));
      fm.template registerEvaluator<EvalT>(op);
    }
  }  
}


///////////////////////////////////////////////////////////////////////////////
//
//  
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
initDynamicTrapsParams(Teuchos::RCP<const Teuchos::ParameterList> plist)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
    
  // dynamic traps parameter list
  if (plist->isSublist("Dynamic Traps"))
  { 
    const ParameterList& trapsPList = plist->sublist("Dynamic Traps"); 
    dynTrapsPList = rcp(new ParameterList(trapsPList));
 
    // check if electric field needs to be passed 
    withField = false;
    for (auto itr = trapsPList.begin(); itr != trapsPList.end(); ++itr) {
      const Teuchos::ParameterEntry& entry = itr->second;
      const ParameterList& trapPList = Teuchos::getValue<ParameterList>(entry);
      double e_xdep = 0.0;
      if (trapPList.isParameter("Electron Electric Field Power Dependency")) {
	e_xdep = trapPList.get<double>("Electron Electric Field Power Dependency");
	if (e_xdep > 0.0) {
	  withField = true; break;
	}
      }
      double h_xdep = 0.0;
      if (trapPList.isParameter("Hole Electric Field Power Dependency")) {
	h_xdep = trapPList.get<double>("Hole Electric Field Power Dependency");
	if (h_xdep > 0.0) {
	  withField = true; break;
	}
      }
    } // electric field check

  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<Teuchos::ParameterList>
charon::BCStrategy_Neumann_DynamicTraps<EvalT>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<ParameterList> p = Teuchos::rcp(new ParameterList);

  ParameterList& trapsPL = p->sublist("Dynamic Traps", false, 
			    "Sublist defining Dynamic Traps");

  for (int i = 0; i < MAX_NUM_TRAPS; i++)
  {  
    std::stringstream ss;
    ss << i;
    string subListName("Trap " + ss.str());
    trapsPL.sublist(subListName, false, 
		   "Sublist defining the parameters for an individual trap");

    // check trap parameters
    trapsPL.sublist(subListName).set<string>("Trap Type", "", "Either Acceptor or Donor");
    trapsPL.sublist(subListName).set<double>("Energy Level", 0.0, 
		  "Trap energy level measured from the band edge in [eV]");
    trapsPL.sublist(subListName).set<double>("Trap Density", 0.0, 
		  "Trap density in [cm^-2]");
    trapsPL.sublist(subListName).set<string>("Energy Distribution", "", 
		  "Trap energy distribution type");
    trapsPL.sublist(subListName).set<string>("Thermal Velocity Calculation", "", 
		  "Either Mean or Root Mean Square");
    trapsPL.sublist(subListName).set<int>("Number of Levels", 20, 
		  "Number of discrete energy levels for continuous distribution");
    trapsPL.sublist(subListName).set<double>("Energy Width", 0.0, 
		  "Energy distribution width in [eV]");
    trapsPL.sublist(subListName).set<double>("Degeneracy Factor", 0.0, 
		  "Degeneracy factor unitless");
    trapsPL.sublist(subListName).set<double>("Electron Cross Section", 0.0, 
		  "Electron capture cross section in [cm^2]");
    trapsPL.sublist(subListName).set<double>("Hole Cross Section", 0.0, 
		  "Hole capture cross section in [cm^2]");
    trapsPL.sublist(subListName).set<double>("Electron Electric Field Power Dependency",0.0,
		  "Electric Field Power Dependency Exponenent for Electron Capture Cross Section");
    trapsPL.sublist(subListName).set<double>("Hole Electric Field Power Dependency",0.0,
		  "Electric Field Power Dependency Exponenent for Hole Capture Cross Section");
    trapsPL.sublist(subListName).set<string>("Electron Field Dependence", "", 
		  "Electron Capture Cross Section field dependency type");
    trapsPL.sublist(subListName).set<string>("Hole Field Dependence", "", 
		  "Hole Capture Cross Section field dependency type");
    trapsPL.sublist(subListName).set<double>("X Min", 0.0, 
                  "X min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("X Max", 0.0, 
                  "X max for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Y Min", 0.0, 
                  "Y min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Y Max", 0.0, 
                  "Y max for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Z Min", 0.0, 
                  "Z min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Z Max", 0.0, 
                  "Z max for the spatial box containing trap");
  }
  
  return p;
}


///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterGatherAndOrientationEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::PhysicsBlock& pb,
                                               const panzer::LinearObjFactory<panzer::Traits> & lof,
                                               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Gather
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_DynamicTraps<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
