

#ifndef CHARON_BC_STRATEGY_NEUMANN_SCHOTTKYCONTACT_IMPL_HPP
#define CHARON_BC_STRATEGY_NEUMANN_SCHOTTKYCONTACT_IMPL_HPP

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
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"

#include "Charon_Names.hpp"
#include "Charon_BCStrategy_Neumann_SchottkyContact.hpp"
#include "Charon_Scaling_Parameters.hpp"

// evaluators
#include "Charon_BC_NeumannSchottkyContact.hpp"
#include "Charon_Integrator_HJFluxDotNorm.hpp"

// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
BCStrategy_Neumann_SchottkyContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Neumann Schottky Contact" );
}

template <typename EvalT>
void charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
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
  RCP<const ParameterList> pbParamList = side_pb.getParameterList();
  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // get the object of charon::Names
  auto names = rcp(new charon::Names(1,prefix,discfields,discsuffix));

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  if (not dataPList->isParameter("Type")) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	<< "Schottky Contact Error: Contact must have electron or hole type!" <<  "\n");
  } else {
    if (dataPList->template get<string>("Type") == "Electron")
      contact_type = -1;
    else if (dataPList->template get<string>("Type") == "Hole")
      contact_type = 1;
    else 
      contact_type = 0;
    if(contact_type == 0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	  << "Schottky Contact Error: Contact must have electron or hole type!" <<  "\n");
  }

  // barrier
  if (not dataPList->isParameter("Work Function")) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	  << "Schottky Contact Error: Contact must specify a work function!" <<  "\n");
  } else {
    Wf = dataPList->template get<double>("Work Function");
    if(Wf <= 0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	<< "Schottky Contact Error: Wf must be positive value!" <<  "\n");
  }

  An = 0.0; Ap = 0.0;
  // Electron Effective Richardson Constant in [A/K^2-cm^2] 
  if (dataPList->isParameter("Electron Richardson Constant") ) 
    An = dataPList->template get<double>("Electron Richardson Constant");
  // Hole Effective Richardson Constant in [A/K^2-cm^2] 
  if (dataPList->isParameter("Hole Richardson Constant") ) 
    Ap = dataPList->template get<double>("Hole Richardson Constant");
  if (An <= 0.0) 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	<< "Schottky Contact Error: 'An' must be defined and have a positive value!" <<  "\n");
  if (Ap <= 0.0) 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
	<< "Schottky Contact Error: 'Ap' must be defined and have a positive value!" <<  "\n");   

  withBL = false;
  alpha = 1.0, beta = 0.0, gamma = 1.0;
  if (dataPList->isSublist("Barrier Lowering")) {
    withBL = true;
    const ParameterList& BL = dataPList->sublist("Barrier Lowering");
    if (BL.isParameter("alpha"))
      alpha = BL.get<double>("alpha");
    if (BL.isParameter("beta"))
      beta = BL.get<double>("beta");
    if (BL.isParameter("gamma"))
      gamma = BL.get<double>("gamma");
  }

  withTunneling = false;
  m_tun = 1.0;
  if (dataPList->isSublist("Tunneling")) {
    withTunneling = true;
    const ParameterList& tun = dataPList->sublist("Tunneling");
    if (tun.isParameter("mass"))
      m_tun = tun.get<double>("mass");
  }

  // check if Equation Set Name = ELECTRON_DENSITY HOLE_DENSITY or HOLE_DENSITY ELECTRON_DENSITY
  string eqsn = this->m_bc.equationSetName();

  if ( (eqsn != "ELECTRON_DENSITY HOLE_DENSITY") && (eqsn != "HOLE_DENSITY ELECTRON_DENSITY") )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Schottky Contact Error: Equation Set Name must be either ELECTRON_DENSITY HOLE_DENSITY or "
        << "HOLE_DENSITY ELECTRON_DENSITY !" << "\n");

  // For now, assume that only one IR is used for everything.
  const auto& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);  // Doesn't support mixed ir yet!

  // obtain the integration order
  const int integration_order = ir.begin()->second->order();

  // get all dofs 
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();
  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
       dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    // for ELECTRON_DENSITY and HOLE_DENSITY
    if ((dof_it->first == names->dof.edensity) || (dof_it->first == names->dof.hdensity)) {
      const string residual_name = "Residual_" + dof_it->first; 
      const string flux_name = (dof_it->first == names->dof.edensity) 
	                          ? "SchottkySurface_eCurrent" 
	                          : "SchottkySurface_hCurrent";
      this->addResidualContribution(residual_name,dof_it->first,flux_name,
				    integration_order,side_pb);
    }
  }  // end of loop over dofs
}



template <typename EvalT>
void charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
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
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? 
    eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? 
    eqSetPList.get<string>("Discontinuous Suffix") : "";

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

  // compute fluxes at Schottky contact
  {
    ParameterList p("NeumannBC Schottky Current");
    if (isSG) { // SG
      p.set("Data Layout", cvfem_bc_ir->dl_scalar);  
    } else { // non-SG, i.e., FEM or EFFPG
      p.set("Data Layout", ir->dl_scalar);
    }  
    p.set("Names", names);
    p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
    p.set("Electron Flux Name", "SchottkySurface_eCurrent");
    p.set("Hole Flux Name", "SchottkySurface_hCurrent");
    p.set("Contact Type", (contact_type == -1) ? "Electron" : "Hole" );
    p.set("An", An);
    p.set("Ap", Ap);
    p.set("Work Function", Wf);

    if (withBL) {
      p.set("BL_alpha", alpha);
      p.set("BL_beta", beta);
      p.set("BL_gamma", gamma);

      // compute normals to the Schottky contact
      string normal_name = pb.cellData().side() + "_SchottkyNorm";
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
      string normal_dot_grad_name = "SchottkyNormDotGrad";
      ParameterList p2(normal_dot_grad_name);
      p2.set("Result Name", normal_dot_grad_name);
      p2.set("Vector A Name", names->grad_dof.phi);
      p2.set("Vector B Name", normal_name);
      p2.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
      const RCP< PHX::Evaluator<panzer::Traits> >
	op2 = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p2));
      this->template registerEvaluator<EvalT>(fm, op2);

      // pass electric field name to the BC_NeumannSchottkyContact evaluator
      p.set("EdotNorm", normal_dot_grad_name);
    }

    if (withTunneling) {
      p.set("tun_m", m_tun);
      if (this->m_bc.params()->isParameter("Voltage")) {
	p.setEntry("Voltage", this->m_bc.params()->getEntry("Voltage"));
      } else if (this->m_bc.params()->isParameter("Varying Voltage")) {
	p.setEntry("Varying Voltage",
		   this->m_bc.params()->getEntry("Varying Voltage"));
      } else {
	p.set<double>("Voltage", 0);
      }
      p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);
      
      if(not withBL) {
	// compute normals to the Schottky contact
	string normal_name = pb.cellData().side() + "_SchottkyNorm";
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
	string normal_dot_grad_name = "SchottkyNormDotGrad";
	ParameterList p2(normal_dot_grad_name);
	p2.set("Result Name", normal_dot_grad_name);
	p2.set("Vector A Name", names->grad_dof.phi);
	p2.set("Vector B Name", normal_name);
	p2.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
	const RCP< PHX::Evaluator<panzer::Traits> >
	  op2 = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p2));
	this->template registerEvaluator<EvalT>(fm, op2);

	// pass electric field name to the BC_NeumannSchottkyContact evaluator
	p.set("EdotNorm", normal_dot_grad_name);
      }
    }

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::BC_NeumannSchottkyContact<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  } 
   
  // create the residuals contributors
  for (std::size_t id = 0; id < vec_dof_name.size(); ++id) {
    const string& residual_name = vec_residual_name[id];
    const string& flux_name = vec_flux_name[id]; 
    double multiplier = 1.0; 
    if (isSG) { // SG  
      ParameterList p(residual_name);
      p.set("Residual Name", residual_name);
      p.set("Flux Name", flux_name);
      p.set("Basis", hgrad_bc_cvfem);
      p.set("IR", cvfem_bc_ir);
      p.set("Multiplier", multiplier);
      const RCP< PHX::Evaluator<panzer::Traits> >
	op = rcp(new charon::Integrator_HJFluxDotNorm<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    } else { // FEM
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



template <typename EvalT>
void charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
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


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_SchottkyContact<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
