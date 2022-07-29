

#ifndef CHARON_BC_GATETUNNELING_IMPL_HPP
#define CHARON_BC_GATETUNNELING_IMPL_HPP
 
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

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_GateTunnelingCurrentDensity.hpp"
#include "Charon_Integrator_HJFluxDotNorm.hpp"
#include <vector>


///////////////////////////////////////////////////////////////////////////////
//
//  BCStrategy_GateTunneling()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
charon::BCStrategy_GateTunneling<EvalT>::
BCStrategy_GateTunneling(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Gate Tunneling" );
}


///////////////////////////////////////////////////////////////////////////////
//
//  setup()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_GateTunneling<EvalT>::
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
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? 
    eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? 
    eqSetPList.get<string>("Discontinuous Suffix") : "";

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

  string gate_cnt = dataPList->isParameter("Gate Sideset ID") 
    ? dataPList->get<string>("Gate Sideset ID") : "";
  if (gate_cnt == "") {
    std::stringstream msg;
    msg << "'Gate Sideset ID' must be specified for '"
	<< this->m_bc.sidesetID() << "' sideset!" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  if (!dataPList->isParameter("Gate Distance")) {
    std::stringstream msg;
    msg << "'Gate Distance' must be specified for '"
	<< this->m_bc.sidesetID() << "' sideset!" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  } 

  bool with_eTunneling = false;
  if (dataPList->isParameter("Electron Tunneling"))
    if (dataPList->get<string>("Electron Tunneling") == "True")
      with_eTunneling = true;

  bool with_hTunneling = false;
  if (dataPList->isParameter("Hole Tunneling"))
    if (dataPList->get<string>("Hole Tunneling") == "True")
      with_hTunneling = true;

  if ((with_eTunneling == false) && (with_hTunneling == false)) {
    std::stringstream msg;
    msg << "Electron Tunneling and Hole Tunneling cannot be both 'False' for '" << this->m_bc.sidesetID() << 
	  "' sideset!" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

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
       dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    dof_name = dof_it->first; 
    // for ELECTRON_DENSITY and HOLE_DENSITY
    if ((dof_name == names->dof.edensity && with_eTunneling) || 
        (dof_name == names->dof.hdensity && with_hTunneling)) {
      residual_name = "Residual_" + dof_name; 
      if (dof_name == names->dof.edensity && with_eTunneling) 
	flux_name = "eGateTunnelingCurrentDensity"; 
      else
	flux_name = "hGateTunnelingCurrentDensity";
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
void charon::BCStrategy_GateTunneling<EvalT>::
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
  using Teuchos::Comm;
  
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

  // get tunneling info
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  /*
  bool with_eTunneling = false;
  if (dataPList->isParameter("Electron Tunneling"))
    if (dataPList->get<string>("Electron Tunneling") == "True")
      with_eTunneling = true;
  bool with_hTunneling = false;
  if (dataPList->isParameter("Hole Tunneling"))
    if (dataPList->get<string>("Hole Tunneling") == "True")
      with_hTunneling = true;
  */
  string gate_cnt = dataPList->get<string>("Gate Sideset ID");
  double gate_distance = dataPList->get<double>("Gate Distance");

  const auto data = this->getResidualContributionData();

  // Iterate over each residual contribution with tunneling activated
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

  // loop over the dof names
  for (std::size_t id = 0; id < vec_dof_name.size(); ++id) {
    // compute the scaled tunneling current at the IPs of the double interface
    {
      ParameterList p("Tunneling Parameters");
      if (isSG) { // SG
	p.set("IR", cvfem_bc_ir);  // for output fields
	p.set("Basis", hgrad_bc_cvfem);  // for input fields
      } else { // non-SG, i.e., FEM or EFFPG
	p.set("IR", ir);  // for output fields
	p.set("Basis", basis);  // for input fields 
      }  
      p.set<RCP<const charon::Names>>("Names", names);
      p.set<string>("Sideset ID", this->m_bc.sidesetID());
      p.set<string>("Gate Sideset ID", gate_cnt);
      p.set<string>("Block ID", ebID_bc);
      p.set<double>("Gate Distance", gate_distance);
      
      p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
      p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);
      p.set<string>("Tunneling Current Density", vec_flux_name[id]);  // output flux name 

      auto op = rcp(new charon::GateTunnelingCurrentDensity<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }

    const string& residual_name = vec_residual_name[id];
    const string& flux_name = vec_flux_name[id]; 
    //const string& dof_name = vec_dof_name[id];
    double multiplier = 1.0; 
    
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
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<Teuchos::ParameterList>
charon::BCStrategy_GateTunneling<EvalT>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<ParameterList> p = Teuchos::rcp(new ParameterList);
  p->set<string>("Gate Sideset ID", "?", "Specifies Gate Sideset ID");
  p->set<double>("Gate Distance", 0.0 , "Gate distance to the tunneling gate in cm");
  p->set<string>("Electron Tunneling", "True", "Enable Electron Tunneling");
  p->set<string>("Hole Tunneling", "True", "Enable Hole Tunneling");  

  return p;
}


///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterGatherAndOrientationEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_GateTunneling<EvalT>::
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
void charon::BCStrategy_GateTunneling<EvalT>::
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
void charon::BCStrategy_GateTunneling<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
