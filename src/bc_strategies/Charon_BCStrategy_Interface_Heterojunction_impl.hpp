

#ifndef CHARON_BC_STRATEGY_INTERFACE_HETEROJUNCTION_IMPL_HPP
#define CHARON_BC_STRATEGY_INTERFACE_HETEROJUNCTION_IMPL_HPP

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

// Evaluators
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Product.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_FieldSpy.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Subtract.hpp"
#include "Charon_Heterojunction_CurrentDensity.hpp"
#include "Charon_Heterojunction_LocalTunneling.hpp"
#include "Charon_Heterojunction_SurfaceCharge.hpp"
#include "Charon_Integrator_HJFluxDotNorm.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Interface_Heterojunction<EvalT>::
BCStrategy_Interface_Heterojunction(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Interface_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Interface Heterojunction" );
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_Heterojunction<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using std::string;
  using Teuchos::RCP;

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // are we setting up the left (side 1) or right side (side 2) of the interface?
  const int di = this->getDetailsIndex();

  // primary side's DOF name
  const string dof_name = (di == 0) ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  // other side's DOF name
  other_dof_name = (di == 1) ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  std::size_t foundedens = dof_name.find("ELECTRON_DENSITY");
  std::size_t foundhdens = dof_name.find("HOLE_DENSITY");
  std::size_t foundpot = dof_name.find("ELECTRIC_POTENTIAL");
  if (foundedens != string::npos || foundhdens != string::npos) {
    // Effective Mass is given in unit of m0 (free electron effective mass)
    double effMass = dataPList->template get<double>("Effective Mass");

    // Obtain physical constants
    charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
    double q = cpc.q;         // [C]
    double kb = cpc.kb * q;   // [J/K]
    double m0 = cpc.m0;       // [kg]
    double h = cpc.h;         // [J.s]
    double pi = cpc.pi;

    // Compute Richardson constant in unit of [Ampere/(cm^2.Kelvin^2)],
    // RC = q*meff*m0*kb^2 / (2*pi^2*hbar^3) = 4*pi*q*meff*m0*kb^2 / h^3
    // 1e4 converts from A/(m^2.K^2) to A/(cm^2.K^2)
    double prefactor = 4.0*pi*q*m0*kb*kb / std::pow(h, 3) / 1e4;  
    richConst = prefactor * effMass;

    // Band Offset is given in unit of [eV] and can be either positive or negative
    // dEc = Ec1 - Ec2, while dEv = Ev2 - Ev1.
    bandOffset = dataPList->template get<double>("Band Offset");
    // obtain the carrier density in unit of [cm^-3] above which the Fermi-Dirac expressions are used
    fdDensity = dataPList->template get<double>("Fermi Dirac Density");

    // determine if local tunneling is turned on
    string localTunnelStr = "False";
    if (dataPList->isParameter("Local Tunneling"))
      localTunnelStr = dataPList->template get<string>("Local Tunneling");

    localTunnel = false;  // default
    if (localTunnelStr == "True") localTunnel = true;

    tunnelMass = effMass; // default
    if (localTunnel) tunnelMass = dataPList->template get<double>("Tunneling Effective Mass");
  } else if (foundpot != string::npos) {
    // not used
    richConst = 0.0;
    bandOffset = 0.0;
    fdDensity = 0.0;
    localTunnel = false;  
    tunnelMass = 0.0;
    if (dataPList->isParameter("Fixed Charge")) {
      surfCharge = dataPList->template get<double>("Fixed Charge");
    } else if(dataPList->isParameter("Fixed Charge")) {
      surfCharge = dataPList->template get<double>("Fixed Charge");
    } else {
      std::stringstream msg;
      msg << "'Fixed Charge' must be specified on '" << this->m_bc.elementBlockID() << 
	  "' side of the heterointerface!" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }
  }
  
  // determine which discretization method is used (either FEM-SUPG or CVFEM-SG)
  discMethod = dataPList->template get<string>("Discretization Method");
  
  // get the physics block parameter list
  RCP<const Teuchos::ParameterList> pbParamList = side_pb.getParameterList();

  // get the equation set parameter list
  const Teuchos::ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any suffix parameters
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // unique name to give the electric potential on this side
  indexer_electric_potential_name = "ELECTRIC_POTENTIAL"; // Should we use charon::names here?
  electric_potential_name = indexer_electric_potential_name + discsuffix;

  // unique residual name, need to use equationSetName(), not dof_name
  // otherwise, try to add residual contribution when di = 1
  const string residual_name = "Residual_" + this->m_bc.equationSetName();

  // define a flux name
  string flux_name = "";
  if (foundedens != string::npos || foundhdens != string::npos) {
    flux_name = "HJTE_CurrentDensity";
  }  else {
    flux_name = "HJTE_ChargeFlux";
  }

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);

  const int integration_order = ir.begin()->second->order();

  if (foundedens != string::npos || foundhdens != string::npos) 
    this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);
  else if (foundpot != string::npos)
    //this->addResidualContribution(residual_name, dof_name, "", integration_order, side_pb);
    this->addResidualContribution(residual_name, dof_name,flux_name, integration_order, side_pb);
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_Heterojunction<EvalT>::
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

  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
    Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  const string residual_name = std::get<0>(data[0]);
  const string dof_name = std::get<1>(data[0]);
  const string flux_name = std::get<2>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name);

  // get cvfem boundary integration rule and basis
  std::string cvfem_type = "boundary";
  RCP<panzer::IntegrationRule> cvfem_bc_ir = Teuchos::rcp(new panzer::IntegrationRule(pb.cellData(),cvfem_type));
  Teuchos::RCP<const panzer::PureBasis> Hgradbasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),ir->topology));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_bc_cvfem = Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_bc_ir));

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // make sure we have at least one discontinuous field
  TEUCHOS_ASSERT( discfields != "" );

  // apply the suffix to the electric potential too
  discfields += ","+indexer_electric_potential_name;

  // create a charon names object
  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // build and register all closure models
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // does it project non-DOFs from nodes to IPs at a boundary ?
  pb.buildAndRegisterDOFProjectionsToIPEvaluators(fm,Teuchos::null,user_data);

  std::size_t foundedens = dof_name.find("ELECTRON_DENSITY");
  std::size_t foundhdens = dof_name.find("HOLE_DENSITY");
  std::size_t foundpot = dof_name.find("ELECTRIC_POTENTIAL");

  // are we setting up the left (side 1) or right side (side 2) of the interface?
  const int di = this->getDetailsIndex();

  string eqnSetName = this->m_bc.equationSetName();

  if (foundedens != string::npos || foundhdens != string::npos) {
    // determine if the primary side (di = 0) is on the left or right
    string primarySide = "";

    if (eqnSetName.find("1") != string::npos)   // left side always corresponds to side 1
      primarySide = "Left";
    else if (eqnSetName.find("2") != string::npos)  // right side always corresponds to side 2
      primarySide = "Right";
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: Equation Set Name must contain 1 or 2 !");

    // define string names
    string flux_name_side1 = "HJTE_CurrentDensity_Side1";
    string flux_name_side2 = "HJTE_CurrentDensity_Side2";
    string local_tunnel_side1 = "HJLT_LocalTunnel_Side1";
    string local_tunnel_side2 = "HJLT_LocalTunnel_Side2";
    string normal_side1 = "Normal_Side1";
    string normal_side2 = "Normal_Side2";
    string grad_phi_side1 = "Potential_Gradient_Side1";
    string grad_phi_side2 = "Potential_Gradient_Side2";
    string normal_dot_grad_side1 = "Normal_Dot_GradPhi_Side1";
    string normal_dot_grad_side2 = "Normal_Dot_GradPhi_Side2";
    string local_tunnel_name = "HJLT_LocalTunnel";
    
    string flux_name_side = "";
    string local_tunnel_side = "";
    string normal_name = "";
    string grad_phi_name = "";
    string normal_dot_grad_name = "";

    double multiplier = 1;

    // ===========================================================================
    // set up string names and the proper multiplier according to the primary side
    // ===========================================================================
    if (primarySide == "Left")   // di=0 for the left side (side 1), di=1 for the right side (side 2)
      {
	flux_name_side = (di == 0) ? flux_name_side1 : flux_name_side2;
	local_tunnel_side = (di == 0)? local_tunnel_side1: local_tunnel_side2;
	normal_name = (di == 0) ? normal_side1 : normal_side2;
	grad_phi_name = (di == 0)? grad_phi_side1: grad_phi_side2;
	normal_dot_grad_name = (di == 0)? normal_dot_grad_side1: normal_dot_grad_side2;

	if (foundedens != string::npos)       // for ELECTRON_DENSITY
	  multiplier = -1.0;                  // Jhj \cdot normal1 > 0
	else if (foundhdens != string::npos)  // for HOLE_DENSITY
	  multiplier = 1.0;                   // Jhj \cdot normal1 > 0
      }

    else if (primarySide == "Right")  // di=1 for the left side (side 1), di=0 for the right side (side 2)
      {
	flux_name_side = (di == 1) ? flux_name_side1 : flux_name_side2;
	local_tunnel_side = (di == 1)? local_tunnel_side1: local_tunnel_side2;
	normal_name = (di == 1) ? normal_side1 : normal_side2;
	grad_phi_name = (di == 1)? grad_phi_side1: grad_phi_side2;
	normal_dot_grad_name = (di == 1)? normal_dot_grad_side1: normal_dot_grad_side2;
	
	if (foundedens != string::npos)       // for ELECTRON_DENSITY
	  multiplier = 1.0;                   // Jhj \cdot normal2 < 0
	else if (foundhdens != string::npos)  // for HOLE_DENSITY
	  multiplier = -1.0;                  // Jhj \cdot normal2 < 0
      }


    // ===========================================================================
    // compute the heterojunction thermionic emission current density
    // ===========================================================================
    {
      // compute the heterojunction thermionic emission current density due to each side
      ParameterList p(flux_name_side);
      p.set<string>("DOF Name", dof_name);
      p.set<string>("Other DOF Name", other_dof_name);
      p.set<string>("Flux Name", flux_name_side);
      p.set<double>("Richardson Constant", richConst);
      p.set<double>("Band Offset", bandOffset);
      p.set<double>("Fermi Dirac Density", fdDensity);
      p.set("Names", names);
      p.set<int>("Details Index", di);
      p.set<string>("Primary Side", primarySide);
      p.set<string>("Discretization Method", discMethod);
      p.set("Scaling Parameters", scaleParams);
      
      if (discMethod == femsupg)
	{
	  p.set("Basis", basis);
	  p.set("IR", ir);
	  p.set("Data Layout", ir->dl_scalar);
	}
      else if (discMethod == cvfemsg)
	{
	  p.set("Basis", hgrad_bc_cvfem);
	  p.set("IR", cvfem_bc_ir);
	  p.set("Data Layout", basis->functional);
	}
      else
	TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
		 <<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);

      const RCP< PHX::Evaluator<panzer::Traits> >
	op = rcp(new charon::Heterojunction_CurrentDensity<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }

    // ===========================================================================
    // compute the heterojunction local tunneling factor
    // ===========================================================================
    if (localTunnel)  // Local tunneling is turned on
      {
	// Get either Normal_Side1 or Normal_Side2
	{
	  ParameterList p(normal_name);
	  p.set("Name", normal_name);
	  p.set("Side ID", pb.cellData().side());

	  if (discMethod == femsupg)
	    p.set("IR", ir);
	  else if (discMethod == cvfemsg)
	    p.set("IR", cvfem_bc_ir);
	  else
	    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
		  <<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);

	  p.set("Normalize", true);
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
	  this->template registerEvaluator<EvalT>(fm, op);
	}
	
	// Get either Potential_Gradient_Side1 or Potential_Gradient_Side2
	{
	  ParameterList p(grad_phi_name);
	  p.set("Name", electric_potential_name);
	  p.set("Gradient Name", grad_phi_name);

	  if (discMethod == femsupg)
	    {
	      p.set("Basis", basis);
	      p.set("IR", ir);
	    }
	  else if (discMethod == cvfemsg)
	    {
	      p.set("Basis", hgrad_bc_cvfem);
	      p.set("IR", cvfem_bc_ir);
	    }
	  else
	    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
		<<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);

	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
	  this->template registerEvaluator<EvalT>(fm, op);
	}

	// Compute the side normal dot the corresponding potential gradient
	{
	  ParameterList p(normal_dot_grad_name);
	  p.set("Result Name", normal_dot_grad_name);
	  p.set("Vector A Name", grad_phi_name);
	  p.set("Vector B Name", normal_name);
	  p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
	  this->template registerEvaluator<EvalT>(fm, op);
	}

	// Compute the local tunneling factor
	{
	  ParameterList p(local_tunnel_name);
	  p.set<double>("Band Offset", bandOffset);
	  p.set<double>("Tunneling Effective Mass", tunnelMass);
	  p.set<int>("Details Index", di);
	  p.set<string>("Primary Side", primarySide);
	  p.set("Names", names);
	  p.set<string>("Flux Name", local_tunnel_name);
	  p.set("Scaling Parameters", scaleParams);

	  if (discMethod == femsupg)
	    p.set("IR", ir);
	  else if (discMethod == cvfemsg)
	    p.set("IR", cvfem_bc_ir);
	  else
	    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
	       <<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);

	  if (bandOffset > 0.0)  // use the normal_dot_grad_side1 from side 1, independent of primarySide and di;
	    // the same holds for both ELECTRON_DENSITY and HOLE_DENSITY
	    p.set<string>("Normal Dot Gradient", normal_dot_grad_side1);

	  else  // use the normal_dot_grad_side2 from side 2, independent of primarySide and di;
	    p.set<string>("Normal Dot Gradient", normal_dot_grad_side2);

	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new charon::Heterojunction_LocalTunneling<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}
      }

    // ==============================================================================
    // obtain the net TE current density and add HJ flux contribution to dof residual
    // ==============================================================================

    if (di == 0)  // for ELECTRON_ or HOLE_DENSITY1 when primarySide = Left,
      {             // for ELECTRON_ or HOLE_DENSITY2 when primarySide = Right

	// subtract the current density from both sides to obtain a net thermionic emission current density
	{
	  ParameterList p(flux_name);
	  p.set<string>("Difference Name", flux_name);

	  if (discMethod == femsupg)
	    p.set("Data Layout", ir->dl_scalar);
	  else if (discMethod == cvfemsg)
	    p.set("Data Layout", cvfem_bc_ir->dl_scalar);
	  
	  if (foundedens != string::npos)  // for ELECTRON_DENSITY
	    {
	      p.set<string>("Value A", flux_name_side2);
	      p.set<string>("Value B", flux_name_side1);
	    }
	  else if (foundhdens != string::npos)  // for HOLE_DENSITY
	    {
	      p.set<string>("Value A", flux_name_side1);
	      p.set<string>("Value B", flux_name_side2);
	    }

	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new charon::Subtract<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}

	// add HJ flux contribution to the proper dof residual
	if (discMethod == femsupg)   // for FEM-SUPG discretization
	  {
	    using panzer::EvaluatorStyle;
	    using panzer::Integrator_BasisTimesScalar;
	    using panzer::Traits;
	    using PHX::Evaluator;
	    using std::string;
	    using std::vector;
	    vector<string> fieldMultipliers;
	    if (localTunnel)
	      fieldMultipliers.push_back(local_tunnel_name);
	    const RCP<Evaluator<Traits>> op = rcp(new
		    Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
		    residual_name, flux_name, *basis, *ir, multiplier, fieldMultipliers));
	    fm.template registerEvaluator<EvalT>(op);
	  }

	else if (discMethod == cvfemsg)  // for CVFEM-SG discretization
	  {
	    ParameterList p(residual_name);
	    p.set("Residual Name", residual_name);
	    p.set("Flux Name", flux_name);
	    p.set("Basis", hgrad_bc_cvfem);
	    p.set("IR", cvfem_bc_ir);
	    p.set("Multiplier", multiplier);

	    if (localTunnel)  // local tunneling is turned on
	      {
		RCP<std::vector<std::string> > fms = Teuchos::rcp(new std::vector<std::string>);
		fms->push_back(local_tunnel_name);
		p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);
	      }

	    const RCP< PHX::Evaluator<panzer::Traits> >
	      op = rcp(new charon::Integrator_HJFluxDotNorm<EvalT,panzer::Traits>(p));
	    fm.template registerEvaluator<EvalT>(op);
	  }
      }  // end of if (di == 0)
    // end ELECTRON_DENSITY or HOLE_DENSITY

  } else if (foundpot != string::npos) { // ELECTRIC_POTENTIAL
    // primary side's DOF name
    const string dof_name = (di == 0) ? 
      this->m_bc.equationSetName() : this->m_bc.equationSetName2();
    // other side's DOF name
    const string other_dof_name = (di == 1) ? 
      this->m_bc.equationSetName() : this->m_bc.equationSetName2();

    // for continuous potential at interface
    if ((dof_name == "ELECTRIC_POTENTIAL") and (dof_name == other_dof_name)) {
      if (di == 0) {
	string surface_charge_flux = flux_name;
	// add half charge flux to the residual
	{
	  // compute the flux charge on interface
	  ParameterList p("Heterojunction Surface Charge");
	  if (discMethod == cvfemsg) {
	    p.set("Output Data Layout", cvfem_bc_ir->dl_scalar);  // for output fields
	    p.set("Basis", hgrad_bc_cvfem);  // for input fields
	  }  else { // non-SG, i.e., FEM or EFFPG
	    p.set("Output Data Layout", ir->dl_scalar);  // for output fields
	    p.set("Basis", basis);  // for input fields 
	  }  
	  p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
	  p.set<string>("Flux Surface Charge", surface_charge_flux);
	  p.set<double>("Fixed Charge", surfCharge/2.0);
	  p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);
	  auto op = rcp(new charon::Heterojunction_SurfaceCharge<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}

	// add HJ charge flux contribution to the proper dof residual
	double multiplier = -1.0; 
	if (discMethod == femsupg) {
	  // for FEM-SUPG discretization
	  using panzer::EvaluatorStyle;
	  using panzer::Integrator_BasisTimesScalar;
	  using panzer::Traits;
	  using PHX::Evaluator;
	  using std::string;
	  using std::vector;
	  const RCP<Evaluator<Traits>> op = rcp(
	     new Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
	     residual_name, flux_name, *basis, *ir, multiplier));
	  fm.template registerEvaluator<EvalT>(op);
	} else if (discMethod == cvfemsg) {
	  // for CVFEM-SG discretization 
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
      } // di==0
    } else if ((dof_name == "ELECTRIC_POTENTIAL1" and other_dof_name == "ELECTRIC_POTENTIAL2") or
	       (dof_name == "ELECTRIC_POTENTIAL2" and other_dof_name == "ELECTRIC_POTENTIAL1") ) {
      // discontinuous potential
      /*
      // determine if the primary side (di = 0) is on the left or right
      string primarySide = "";
      if (eqnSetName.find("1") != string::npos)   // left side always corresponds to side 1
	primarySide = "Left";
      else if (eqnSetName.find("2") != string::npos)  // right side always corresponds to side 2
	primarySide = "Right";
      else
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: Equation Set Name must contain 1 or 2 !");
      
      // define string names
      string flux_name_side1 = "HJTE_ChargeFlux_Side1";
      string flux_name_side2 = "HJTE_ChargeFlux_Side2";
      string normal_side1 = "Norm_Side1";
      string normal_side2 = "Norm_Side2";
      string grad_phi_side1 = "Potential_Gradient_Side1";
      string grad_phi_side2 = "Potential_Gradient_Side2";
      string normal_dot_grad_side1 = "Norm_Dot_GradPhi_Side1";
      string normal_dot_grad_side2 = "Norm_Dot_GradPhi_Side2";
      string eps_normal_dot_grad_side1 = "EpsNorm_Dot_GradPhi_Side1";
      string eps_normal_dot_grad_side2 = "EpsNorm_Dot_GradPhi_Side2";
      string surface_charge_flux_side1 = "SurfChargeFlux_Side1";
      string surface_charge_flux_side2 = "SurfChargeFlux_Side2";
      
      string flux_name_side = "";
      string normal_name = "";
      string grad_phi_name = "";
      string normal_dot_grad_name = "";
      string eps_normal_dot_grad_name = "";
      string other_eps_normal_dot_grad_name = "";
      string surface_charge_flux = "";
      
      double multiplier = 1.0;

      // set up string names and the proper multiplier according to the primary side
      if (primarySide == "Left")  {
	// di=0 for the left side (side 1), di=1 for the right side (side 2)
	flux_name_side = (di == 0) ? flux_name_side1 : flux_name_side2;
	normal_name = (di == 0) ? normal_side1 : normal_side2;
	grad_phi_name = (di == 0) ? grad_phi_side1: grad_phi_side2;
	normal_dot_grad_name = (di == 0) ? normal_dot_grad_side1 : normal_dot_grad_side2;
	eps_normal_dot_grad_name = (di == 0) ? eps_normal_dot_grad_side1 : eps_normal_dot_grad_side2;
	other_eps_normal_dot_grad_name = (di == 0) ? eps_normal_dot_grad_side2 : eps_normal_dot_grad_side1;
	surface_charge_flux = (di == 0) ? surface_charge_flux_side1 : surface_charge_flux_side2;
	multiplier = 1.0;                   
      } else if (primarySide == "Right")  {
	// di=1 for the left side (side 1), di=0 for the right side (side 2)
	flux_name_side = (di == 1) ? flux_name_side1 : flux_name_side2;
	normal_name = (di == 1) ? normal_side1 : normal_side2;
	grad_phi_name = (di == 1) ? grad_phi_side1: grad_phi_side2;
	normal_dot_grad_name = (di == 1) ? normal_dot_grad_side1 : normal_dot_grad_side2;
	eps_normal_dot_grad_name = (di == 1) ? eps_normal_dot_grad_side1 : eps_normal_dot_grad_side2;
	other_eps_normal_dot_grad_name = (di == 1) ? eps_normal_dot_grad_side2 : eps_normal_dot_grad_side1;
	surface_charge_flux = (di == 1) ? surface_charge_flux_side1 : surface_charge_flux_side2;
	multiplier = 1.0;    
      }
      
      {
	// get either Normal_Side1 or Normal_Side2
	ParameterList p(normal_name);
	p.set("Name", normal_name);
	p.set("Side ID", pb.cellData().side());
	if (discMethod == femsupg)
	  p.set("IR", ir);
	else if (discMethod == cvfemsg)
	  p.set("IR", cvfem_bc_ir);
	else
	  TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
	       <<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);
	p.set("Normalize", true);
	const RCP< PHX::Evaluator<panzer::Traits> >
	  op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
	this->template registerEvaluator<EvalT>(fm, op);
      }
       
      {
	// get either Potential_Gradient_Side1 or Potential_Gradient_Side2
	ParameterList p(grad_phi_name);
	p.set("Name", electric_potential_name);
	p.set("Gradient Name", grad_phi_name);
	if (discMethod == femsupg) {
	  p.set("Basis", basis);
	  p.set("IR", ir);
	} else if (discMethod == cvfemsg) {
	  p.set("Basis", hgrad_bc_cvfem);
	  p.set("IR", cvfem_bc_ir);
	} else {
	  TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
	    <<"Invalid Discretization Method ! Must be either FEM-SUPG or CVFEM-SG !" << std::endl);
	}
	const RCP< PHX::Evaluator<panzer::Traits> >
	  op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
	this->template registerEvaluator<EvalT>(fm, op);
      }

      if (di == 0) {
	string surface_charge_flux = flux_name;
	// add charge flux to the residual
	// for ELECTRIC_POTENTIAL1 when primarySide = Left
	// for ELECTRIC_POTENTIAL2 when primarySide = Right
	{
	  // compute the flux charge on interface
	  ParameterList p("Heterojunction Surface Charge");
	  if (discMethod == cvfemsg) {
	    p.set("Output Data Layout", cvfem_bc_ir->dl_scalar);  // for output fields
	    p.set("Basis", hgrad_bc_cvfem);  // for input fields
	  }  else { // non-SG, i.e., FEM or EFFPG
	    p.set("Output Data Layout", ir->dl_scalar);  // for output fields
	    p.set("Basis", basis);  // for input fields 
	  }  
	  p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
	  p.set<string>("Flux Surface Charge", surface_charge_flux);
	  p.set<double>("Fixed Charge", surfCharge);
	  p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);
	  auto op = rcp(new charon::Heterojunction_SurfaceCharge<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}

	// add HJ charge flux contribution to the proper dof residual
	if (discMethod == femsupg) {
	  // for FEM-SUPG discretization
	  using panzer::EvaluatorStyle;
	  using panzer::Integrator_BasisTimesScalar;
	  using panzer::Traits;
	  using PHX::Evaluator;
	  using std::string;
	  using std::vector;
	  const RCP<Evaluator<Traits>> op = rcp(
	     new Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
	     residual_name, flux_name, *basis, *ir), 0.5);
	  fm.template registerEvaluator<EvalT>(op);
	} else if (discMethod == cvfemsg) {
	  // for CVFEM-SG discretization 
	  ParameterList p(residual_name);
	  p.set("Residual Name", residual_name);
	  p.set("Flux Name", flux_name);
	  p.set("Basis", hgrad_bc_cvfem);
	  p.set("IR", cvfem_bc_ir);
	  p.set("Multiplier", 0.5);
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new charon::Integrator_HJFluxDotNorm<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	} 

	// add charge flux to the residual
	// for ELECTRIC_POTENTIAL1 when primarySide = Left
	// for ELECTRIC_POTENTIAL2 when primarySide = Right
	{
	  // compute the flux charge on interface
	  ParameterList p("Heterojunction Surface Charge");
	  if (discMethod == cvfemsg) {
	    p.set("Output Data Layout", cvfem_bc_ir->dl_scalar);  // for output fields
	    p.set("Basis", hgrad_bc_cvfem);  // for input fields
	  }  else { // non-SG, i.e., FEM or EFFPG
	    p.set("Output Data Layout", ir->dl_scalar);  // for output fields
	    p.set("Basis", basis);  // for input fields 
	  }  
	  p.set<RCP<charon::Scaling_Parameters>>("Scaling Parameters", scaleParams);
	  p.set<string>("Flux Surface Charge", surface_charge_flux);
	  p.set<double>("Fixed Charge", surfCharge);
	  p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);
	  auto op = rcp(new charon::Heterojunction_SurfaceCharge<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}
      
	{
	  // compute total charge flux 
	  ParameterList p(flux_name);
	  p.set<string>("Sum Name", flux_name);
	  if (discMethod == femsupg)
	    p.set("Data Layout", ir->dl_scalar);
	  else if (discMethod == cvfemsg)
	    p.set("Data Layout", cvfem_bc_ir->dl_scalar);
	  std::vector<std::string> values_names(2);
	  values_names[0] = surface_charge_flux;
	  values_names[1] = other_eps_normal_dot_grad_name;
	  p.set< Teuchos::RCP<std::vector<std::string> > >(
		"Values Names", Teuchos::rcp(new std::vector<std::string>(values_names)));
	  std::vector<double> scalars(2);
	  scalars[0] = 1.0;
	  scalars[1] = -scaleParams->scale_params.Lambda2*multiplier;  
	  
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}
  
	// add HJ charge flux contribution to the proper dof residual
	if (discMethod == femsupg) {
	  // for FEM-SUPG discretization
	  using panzer::EvaluatorStyle;
	  using panzer::Integrator_BasisTimesScalar;
	  using panzer::Traits;
	  using PHX::Evaluator;
	  using std::string;
	  using std::vector;
	  const RCP<Evaluator<Traits>> op = rcp(
	     new Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
	     residual_name, flux_name, *basis, *ir));
	  fm.template registerEvaluator<EvalT>(op);
	} else if (discMethod == cvfemsg) {
	  // for CVFEM-SG discretization 
	  ParameterList p(residual_name);
	  p.set("Residual Name", residual_name);
	  p.set("Flux Name", flux_name);
	  p.set("Basis", hgrad_bc_cvfem);
	  p.set("IR", cvfem_bc_ir);
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new charon::Integrator_HJFluxDotNorm<EvalT,panzer::Traits>(p));
	  fm.template registerEvaluator<EvalT>(op);
	}
      } else if (di == 1) {
	{
	  // compute the side normal dot the corresponding potential gradient
	  ParameterList p(normal_dot_grad_name);
	  p.set("Result Name", normal_dot_grad_name);
	  p.set("Vector A Name", grad_phi_name);
	  p.set("Vector B Name", normal_name);
	  p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
	  this->template registerEvaluator<EvalT>(fm, op);
	}
	
	{ // multiply with eps_rel
	  ParameterList p(eps_normal_dot_grad_name);
	  p.set("Product Name", eps_normal_dot_grad_name);
	  std::vector<std::string> values_names(2);
	  values_names[0] = normal_dot_grad_name;
	  values_names[1] = names->field.rel_perm;
	  p.set< Teuchos::RCP<std::vector<std::string> > >(
	     "Values Names", Teuchos::rcp(new std::vector<std::string>(values_names)));
	  if (discMethod == femsupg)
	    p.set("Data Layout", ir->dl_scalar);
	  else if (discMethod == cvfemsg)
	    p.set("Data Layout", cvfem_bc_ir->dl_scalar);
	  //p.set<double>("Scaling", -scaleParams->scale_params.Lambda2);
	  const RCP< PHX::Evaluator<panzer::Traits> >
	    op = rcp(new panzer::Product<EvalT,panzer::Traits>(p));
	  this->template registerEvaluator<EvalT>(fm, op);
	} 
      } // di == 1 
      */
      ;
    } // discontinuous potential


  } // ELECTRIC_POTENTIAL
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_Heterojunction<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::PhysicsBlock& pb,
                                               const panzer::LinearObjFactory<panzer::Traits> & lof,
                                               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // First do the standard gather of all the DOFs within this physics block
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm, lof, user_data);

  // Get tangent fields
  const std::vector<panzer::StrPureBasisPair> tangent_fields = pb.getTangentFields();

  // We don't support tangent fields yet
  TEUCHOS_ASSERT( 0 == tangent_fields.size() );

  // Second do an extra gather for any continuous fields needed on both sides of the interface

  // Get the vector of residual contributions data
  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
    Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  // Charon assumes that there is only one residual contribution
  TEUCHOS_ASSERT( 1 == data.size() );

  // Get the basis from the tuple
  const RCP<const panzer::PureBasis> basis = std::get<4>(data[0]);

  // Register the extra gather evaluator
  {
    ParameterList p("Extra Gather");
    p.set("Basis", basis);
    RCP<std::vector<string> > dof_names = rcp(new std::vector<string>);
    RCP<std::vector<string> > indexer_names = rcp(new std::vector<string>);
    dof_names->push_back(electric_potential_name);
    indexer_names->push_back(indexer_electric_potential_name);
    p.set("DOF Names", dof_names);
    p.set("Indexer Names", indexer_names);
    p.set("Sensitivities Name", "");
    p.set("First Sensitivities Available", true);

    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);

    this->template registerEvaluator<EvalT>(fm, op);
  }

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_Heterojunction<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_Heterojunction<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
