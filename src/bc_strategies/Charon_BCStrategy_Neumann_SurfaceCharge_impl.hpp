

#ifndef CHARON_BC_STRATEGY_NEUMANN_SURFACECHARGE_IMPL_HPP
#define CHARON_BC_STRATEGY_NEUMANN_SURFACECHARGE_IMPL_HPP
 
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
#include "Charon_NeumannBC_SurfaceCharge.hpp"
#include "Charon_Integrator_HJFluxDotNorm.hpp"
#include <vector>


const int MAX_NUM_TRAPS = 50;


///////////////////////////////////////////////////////////////////////////////
//
//  BCStrategy_Neumann_SurfaceCharge()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
BCStrategy_Neumann_SurfaceCharge(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Neumann Surface Charge" );
  std::cout << "Warning: Neumann Surface Charge with Surface Trap or Surface Recombination does NOT work for a heterojunction!" << std::endl; 
}


///////////////////////////////////////////////////////////////////////////////
//
//  setup()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
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
  if ((eqsn != "ALL_DOFS") && (eqsn != names->dof.phi) )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Equation Set Name must be either ALL_DOFS or ELECTRIC_POTENTIAL"
        << ". But you entered \"" << eqsn << "\" ! \n");

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // validate the Data parameter list
  RCP<ParameterList> valid_params = this->getValidParameters();
  dataPList->validateParameters(*valid_params);

  // call a member function to initialize the input parameters
  initialize(dataPList); 
  
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
      // The Poisson eqn residual needs to be modified when either bFixCharge or bSurfTrap is true 
      if (bFixCharge || bVaryingCharge || bSurfTrap || bPolar)  
      {
        residual_name = "Residual_" + dof_name; 
        flux_name = "Surface_Charge"; 
        fluxSurfCharge = flux_name; 
        this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
      }
    }
    
    // for ELECTRON_DENSITY and HOLE_DENSITY
    if ((dof_name == names->dof.edensity) || (dof_name == names->dof.hdensity)) 
    {
      // The carrier eqn residuals need to be modified when either bSurfTrap or bSurfRecomb is true
      if (bSurfTrap || bSurfRecomb)
      {
        residual_name = "Residual_" + dof_name; 
        flux_name = "Surface_Recombination"; 
        fluxSurfRecomb = flux_name; 
        this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
      }
    }
  }  // end of loop over dofs
  
  
}


///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
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
    ParameterList p("NeumannBC Surface Charge");
    if (isSG)  // SG
    {
      p.set("Output Data Layout", cvfem_bc_ir->dl_scalar);  // for output fields
      p.set("Basis", hgrad_bc_cvfem);  // for input fields
    }  
    else  // non-SG, i.e., FEM or EFFPG
    {
      p.set("Output Data Layout", ir->dl_scalar);  // for output fields
      p.set("Basis", basis);  // for input fields 
    }  
    p.set<RCP<const charon::Names>>("Names", names);
    p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
    p.set<string>("Flux Surface Charge", fluxSurfCharge);        // output flux name 
    p.set<string>("Flux Surface Recombination", fluxSurfRecomb); // output flux name 
    p.set<bool>("Include Fixed Charge", bFixCharge);  // check if Fixed Charge is given 
    p.set<bool>("Include Varying Charge", bVaryingCharge);  // check if Varying Charge is given 
    p.set<bool>("Include Surface Trap", bSurfTrap);   // check if Surface Trap is given 
    p.set<bool>("Include Surface Recombination", bSurfRecomb);  // check if Surface Recombination is given
    p.set<bool>("Include Polarization", bPolar); // check if Polarization is given
    if (bFixCharge) p.set<double>("Fixed Charge", fixedCharge);
    if (bVaryingCharge) p.set<string>("Varying Charge", "Parameter");
    if (bSurfTrap) p.set<RCP<ParameterList>>("Surface Trap ParameterList", surfTrapPList);
    if (bSurfRecomb) p.set<RCP<ParameterList>>("Surface Recombination ParameterList", surfRecombPList); 
    if (bPolar) p.set<RCP<ParameterList>>("Polarization ParameterList", polarPList);
    p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);

    auto op = rcp(new charon::NeumannBC_SurfaceCharge<EvalT,panzer::Traits>(p));
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
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
initialize(Teuchos::RCP<const Teuchos::ParameterList> plist)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  
  // fixed surface charge is given in [#/cm^2], can be positive or negative
  fixedCharge = 0.0;  
  bFixCharge = false; 
  if (plist->isParameter("Fixed Charge")) 
  {
    bFixCharge = true; 
    fixedCharge = plist->template get<double>("Fixed Charge");
  }
  bVaryingCharge = false;
  if (plist->isParameter("Varying Charge")) 
  {
    bVaryingCharge = true; 
  }

  bPolar = false;
  //Polarazation given as wurtzite material either side of interface
  //xcomp is metal composition e.g. Al in AlGaN/GaN interface
  if ( plist->isSublist("Polarization") )
  {
    bPolar = true;
    const auto& spl = plist->sublist("Polarization");
    polarPList = rcp(new ParameterList(spl));
  }
  
  // surface trap parameter list
  bSurfTrap = false; 
  if (plist->isSublist("Surface Trap"))
  {
    bSurfTrap = true; 
    const ParameterList& trapPList = plist->sublist("Surface Trap"); 
    surfTrapPList = rcp(new ParameterList(trapPList));
  }

  // surface recombination parameter list
  bSurfRecomb = false; 
  if (plist->isSublist("Surface Recombination"))
  {
    bSurfRecomb = true; 
    const ParameterList& surfRecPList = plist->sublist("Surface Recombination"); 
    surfRecombPList = rcp(new ParameterList(surfRecPList));
  }
  
  // throw out an error if NONE of Fixed Charge, Polarization, Piezo, 
  //Surface Trap, or Surface Recombination is specified
  if (!bFixCharge && !bVaryingCharge && !bSurfTrap && !bSurfRecomb &&!bPolar)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error!\
     One of Fixed Charge, Polarization, Surface Trap,\
      or Surface Recombination should be specified!" << std::endl);
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<Teuchos::ParameterList>
charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<ParameterList> p = Teuchos::rcp(new ParameterList);
  p->set<double>("Fixed Charge", 0.0, "Fixed surface charge in unit of cm^(-2)"); 
  p->set<std::string>("Varying Charge", "Parameter", "sweeping surface charge in unit of cm^(-2)"); 
  auto& polarPL = p->sublist("Polarization", false, "Polarization sublist");
  polarPL.set<string>("Type", "total", "total or piezo");
  polarPL.set<string>("Top", "", "Top material");
  polarPL.set<string>("Bottom", "", "Bottom material");
  polarPL.set<double>("Xcomp", 0.3, "x composition i.e. AlxGaN");
  polarPL.set<double>("Scale", 1.0, "Scale polarization to help convergence");
  
  ParameterList& trapPL = p->sublist("Surface Trap", false, "Sublist defining Surface Trap");

  for (int i = 0; i < MAX_NUM_TRAPS; i++)
  {
    trapPL.set<double>("Electron Effective Mass", 0.0, "Electron effective mass in unit of m0");
    trapPL.set<double>("Hole Effective Mass", 0.0, "Hole effective mass in unit of m0");
      
    std::stringstream ss;
    ss << i;
    string subListName("Trap " + ss.str());
    trapPL.sublist(subListName, false, "Sublist defining the parameters for one type of trap");

    // trap related parameters
    trapPL.sublist(subListName).set<double>("Trap Energy", 0.0, 
	"Trap energy level measured from the mid-band gap in [eV], + for above, - for below");
    trapPL.sublist(subListName).set<double>("Trap Density", 0.0, "Trap density in [cm^-2] or [cm^-2 eV^-1]");
    trapPL.sublist(subListName).set<string>("Trap Type", "", 
       "Either Acceptor (0 if unoccupied and -1 if occupied, electron capture) or Donor (0 if unoccupied and +1 if occupied, hole capture)");
    trapPL.sublist(subListName).set<string>("Energy Distribution", "", "Energy distribution type");
    trapPL.sublist(subListName).set<double>("Energy Width", 0.0, "Distribution energy width [eV]");
    trapPL.sublist(subListName).set<int>("Number of Levels", 20, 
		  "Number of discrete energy levels for continuous a distribution");
    trapPL.sublist(subListName).set<double>("Electron Cross Section", 0.0, "Electron capture cross section in [cm^2]");
    trapPL.sublist(subListName).set<double>("Hole Cross Section", 0.0, "Hole capture cross section in [cm^2]");
  }
  
  ParameterList& surfRecPL = p->sublist("Surface Recombination", false, "Sublist defining Surface Recombination");
  surfRecPL.set<double>("Electron Surface Velocity", 0.0, "Electron surface recombination velocity in unit of cm/s");
  surfRecPL.set<double>("Hole Surface Velocity", 0.0, "Hole surface recombination velocity in unit of cm/s");
  surfRecPL.set<double>("Energy Level", 0.0, "Trap energy level measured from the mid-band gap in [eV], + for above, - for below");

  return p;
}


///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterGatherAndOrientationEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
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
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
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
void charon::BCStrategy_Neumann_SurfaceCharge<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
