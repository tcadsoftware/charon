
#ifndef CHARON_EQUATIONSET_SGCVFEM_DRIFTDIFFUSION_IMPL_HPP
#define CHARON_EQUATIONSET_SGCVFEM_DRIFTDIFFUSION_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Panzer_Sum.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PoissonSource.hpp"
#include "Charon_SGCVFEM_EdgeCurrDens.hpp"
#include "Charon_SGCVFEM_SubCVCurrDens.hpp"
#include "Charon_SGCVFEM_CentroidCurrDens.hpp"
#include "Charon_SGCVFEM_CentroidDriveForce.hpp"
#include "Charon_SGCVFEM_PotentialFlux.hpp"
#include "Charon_Integrator_SubCVNodeScalar.hpp"
#include "Charon_Integrator_SubCVFluxDotNorm.hpp"


// CVFEM-SG formulation of the Poisson + Continuity equations

// ***********************************************************************
template <typename EvalT>
charon::EquationSet_SGCVFEM_DriftDiffusion<EvalT>::
EquationSet_SGCVFEM_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                   const int& default_integration_order,
                                   const panzer::CellData& cell_data,
                                   const Teuchos::RCP<panzer::GlobalData>& global_data,
                                   const bool build_transient_support) :
  charon::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support)
{

  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;

    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Electron",
      "True",
      "Determines if users want to solve the Electron Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Hole",
      "True",
      "Determines if users want to solve the Hole Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Fermi Dirac",
      "False",
      "Determine if users want to use the Fermi-Dirac statistics",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "FD Formula",
      "Schroeder",  // default
      "Determine if users want to use Schroeder's or the Diffusion method for FD impl.",
      Teuchos::tuple<std::string>("Schroeder","Diffusion"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determines if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determines if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determines if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determines if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Driving Force",
      "EffectiveField",
      "Determines driving force for different models in the continuity equation(s)",
      Teuchos::tuple<std::string>("EffectiveField","GradQuasiFermi","GradPotential"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Defect Cluster",
      "Off",
      "Determines if users want to include Defect Cluster generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Particle Strike",
      "Off",
      "Determine if users want to include ionization due to a particle strike",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Empirical Defect",
      "Off",
      "Determines if users want to include Empirical Defect generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determines if users want to turn on BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Optical Generation",
      "Off",
      "Determine if users want to turn on the optical generation",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Acceptor Incomplete Ionization",
      "Off",
      "Determine if users want to include Acceptor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Donor Incomplete Ionization",
      "Off",
      "Determine if users want to include Donor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Trap SRH",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Trap Charge",
      "False",
      "Determine if users want to include the trap charge in the Poisson equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");

    // Must solve at least one continuity equation
    //if ((solveElectron == "False") && (solveHole == "False"))
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: must solve at least one"
    //       << " continuity equation !" << "\n");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH");
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");
    const std::string& defect_cluster_recomb = params->sublist("Options").get<std::string>("Defect Cluster");
    const std::string& empirical_defect_recomb = params->sublist("Options").get<std::string>("Empirical Defect");
    const std::string& ionization_particle_strike = params->sublist("Options").get<std::string>("Particle Strike");
    const std::string& trap_srh_recomb = params->sublist("Options").get<std::string>("Trap SRH");
    const std::string& add_trap_charge = params->sublist("Options").get<std::string>("Add Trap Charge");
    const std::string& opt_gen = params->sublist("Options").get<std::string>("Optical Generation"); 

    // Must solve both continuity equations when any recomb./gen. model is on.
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") ||
        (trap_srh_recomb == "On") || (ava_gen == "On") || (opt_gen == "On") || 
        (defect_cluster_recomb == "On") || (empirical_defect_recomb == "On") ||
        (ionization_particle_strike == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
      haveSource = true;
    }
    else
      haveSource = false;  // NO recomb. model

    withAvaGen = (haveSource && ava_gen == "On") ? true : false;
    withTrapSRH = (haveSource && trap_srh_recomb == "On") ? true : false;
    addTrapCharge = (trap_srh_recomb == "On" && add_trap_charge == "True") ? true : false;
    drForce = params->sublist("Options").get<std::string>("Driving Force");

    // error check for Incomplete Ionization models
    const std::string& acc_ioniz =
      params->sublist("Options").get<std::string>("Acceptor Incomplete Ionization");
    const std::string& don_ioniz =
      params->sublist("Options").get<std::string>("Donor Incomplete Ionization");
    if(acc_ioniz == "On" && solveHole == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve hole continuity equation when " <<
          "Acceptor Incomplete Ionization is On!" << "\n");
    if(don_ioniz == "On" && solveElectron == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve electron continuity equation when " <<
          "Donor Incomplete Ionization is On!" << "\n");
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int integration_order = params->get<int>("Integration Order");

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

  // *************************************
  // Assemble DOF names and Residual names
  // *************************************

  m_names = Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix));

  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_names));

  // Always solve the Poisson equation
  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);

  // Solve the Electron equation by default and users can disable it
  if (solveElectron == "True")
  {
    this->addDOF(m_names->dof.edensity,basis_type,basis_order,integration_order,m_names->res.edensity);
    this->addDOFGrad(m_names->dof.edensity,m_names->grad_dof.edensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.edensity,m_names->dxdt.edensity);
  }

  // Solve the Hole equation by default and users can disable it
  if (solveHole == "True")
  {
    this->addDOF(m_names->dof.hdensity,basis_type,basis_order,integration_order,m_names->res.hdensity);
    this->addDOFGrad(m_names->dof.hdensity,m_names->grad_dof.hdensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.hdensity,m_names->dxdt.hdensity);
  }

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_SGCVFEM_DriftDiffusion<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<charon::Names> m_names = (this->isFreqDom ? this->base_names : this->m_names);
  const charon::Names& n = (this->isFreqDom ? *(this->base_names) : *(this->m_names));

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // For now, assume that all dofs use the same basis
  Teuchos::RCP<panzer::IntegrationRule> ir = (this->isFreqDom ? this->base_ir : this->getIntRuleForDOF(m_names->dof.phi));
  Teuchos::RCP<panzer::BasisIRLayout> basis = (this->isFreqDom ? this->base_basis : this->getBasisIRLayoutForDOF(m_names->dof.phi));
  // Get basis and IR layout for CVFEM
  panzer::CellData celldata(basis->numCells(),ir->topology);

  std::string cvfem_type = "volume";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_vol_ir =
    Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));
  cvfem_type = "side";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_side_ir =
    Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));

  Teuchos::RCP<const panzer::PureBasis> Hgradbasis =
    Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),ir->topology));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_vol_cvfem =
    Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_vol_ir));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_side_cvfem =
    Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_side_ir));
  Teuchos::RCP<const panzer::PureBasis> Hcurlbasis =
    Teuchos::rcp(new panzer::PureBasis("HCurl",1,basis->numCells(),ir->topology));
  Teuchos::RCP<panzer::BasisIRLayout> hcurl_side_cvfem =
    Teuchos::rcp(new panzer::BasisIRLayout(Hcurlbasis,*cvfem_side_ir));
  Teuchos::RCP<panzer::BasisIRLayout> hcurl_vol_cvfem =
    Teuchos::rcp(new panzer::BasisIRLayout(Hcurlbasis,*cvfem_vol_ir));

  bool haveSource = (this->isFreqDom ? this->base_haveSource : this->haveSource);
  bool addTrapCharge = (this->isFreqDom ? this->base_addTrapCharge : this->addTrapCharge);
  std::string drForce = (this->isFreqDom ? this->base_drForce : this->drForce);

  std::string solveElectron = (this->isFreqDom ? this->base_solveElectron : this->solveElectron);
  std::string solveHole = (this->isFreqDom ? this->base_solveHole : this->solveHole);

  // ***************************************************************************
  // Construct Poisson Equation (solved with continuity equations)
  // based on the CVFEM-SG formulation
  // ***************************************************************************

  // Laplacian Operator: \int_cv_boundary (-\lambda2 \epsilon_r \grad_phi) \cdot {\bf n} dS
  {
    // Compute the potential flux: lambda2*epsilon_r*grad_phi @ the midpoints
    // of subcv edges(2D)/faces(3D) using nodal basis functions (scalar)
    {
      ParameterList p("CVFEM-SG Potential Flux");
      p.set("Flux Name", n.field.phi_flux);
      p.set("DOF Name", n.dof.phi);
      p.set("Basis", hgrad_side_cvfem);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Scaling Parameters", scaleParams);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::SGCVFEM_PotentialFlux<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the CVFEM-SG potential flux over the subcv boundary
    {
      ParameterList p("Laplacian Residual");
      p.set("Residual Name", n.res.phi+n.op.laplacian);
      p.set("Flux Name", n.field.phi_flux);
      p.set<RCP<const charon::Names> >("Names", m_names);
      p.set("Basis", basis);
      p.set("IR", cvfem_side_ir);
      p.set("Multiplier", -1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Source Operator: Assembles \int_cv g dV
  {
    // Compute the Poisson source term: (p-n+Nd-Na) @ BASIS
    {
      ParameterList p("Poisson Source");
      p.set("Source Name", n.field.psrc);
      p.set("Data Layout", basis->functional);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Solve Electron", solveElectron);
      p.set("Solve Hole", solveHole);
      p.set("Scaling Parameters", scaleParams);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::PoissonSource<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Interpolate the source term at primary nodes to the center point of a
    // subcontrol volume and integrate over subcv
    {
      ParameterList p("Source Residual");
      p.set("Residual Name", n.res.phi+n.op.src);
      p.set("Value Name", n.field.psrc);
      p.set("Basis", hgrad_vol_cvfem);
      p.set("IR", cvfem_vol_ir);
      p.set("Multiplier", -1.0);
      p.set("WithInterpolation", true);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }

    // Add the trap charges located at centroids of subcv to the Poisson source residual
    // The trap charges are computed in charon::RecombRate_TrapSRH
    if (addTrapCharge)
    {
      ParameterList p("Trap Source Residual");
      p.set("Residual Name", n.res.phi+"_TRAP_SOURCE_OP");
      p.set("Value Name", n.field.trap_srh_charge);
      p.set("Basis", hgrad_vol_cvfem);
      p.set("IR", cvfem_vol_ir);
      p.set("Multiplier", -1.0);
      p.set("WithInterpolation", false);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    ParameterList p;
    p.set("Sum Name", n.res.phi);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);
    sum_names->push_back(n.res.phi + n.op.laplacian);
    sum_names->push_back(n.res.phi + n.op.src);
    if (addTrapCharge)  sum_names->push_back(n.res.phi + "_TRAP_SOURCE_OP");

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }


  // ***************************************************************************
  // Construct Electron Continuity Equation based on the CVFEM-SG formulation
  // ***************************************************************************

 if (solveElectron == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    // integrate dn_dt over subcv, \int_cv dn_dt dV, where dn_dt = \sum_BFi dni_dt BFi
    // interpolate dn_dt at primary nodes to the center point of subcv and integrate over subcv
    ParameterList p("Electron Transient Residual");
    p.set("Residual Name", n.res.edensity+n.op.trans);
    p.set("Value Name", n.dxdt.edensity);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute electron CVFEM-SG edge current density at the midpoints of primary edges
   {
    ParameterList p("Electron CVFEM-SG Edge Current Density");
    p.set("Carrier Type", "Electron");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_EdgeCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // interpolate electron CVFEM-SG edge current density to the midpoints of subcv
   // edges(2D)/faces(3D) using edge basis functions (vectors)
   {
    ParameterList p("Electron CVFEM-SG SubCV Current Density");
    p.set("Carrier Type", "Electron");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", hcurl_side_cvfem);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_SubCVCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

/*
   // TESTING
   // interpolate electron CVFEM-SG edge current density to centroids
   // of subcontrol volumes using edge basis functions (vectors)
   {
    ParameterList p("Electron CVFEM-SG SubCV Centroid Vector");
    p.set("Carrier Type", "Electron");
    p.set("Vector Name", "SubCV Centroid Vector");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", hcurl_vol_cvfem);
    p.set("IR", cvfem_vol_ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_SubCVCentroidVector<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
*/

   // integrate electron CVFEM-SG current density at the midpoints of subcv
   // edges(2D)/faces(3D) over the subcv boundary
   {
    ParameterList p("Electron Convection Diffusion Residual");
    p.set("Residual Name", n.res.edensity+n.op.conv_diff);
    p.set("Flux Name", n.field.elec_curr_dens_cvedge);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", cvfem_side_ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   if(withAvaGen) 
   {
     // interpolate electron CVFEM-SG edge current density to centroids
     // of subcontrol volumes using edge basis functions (vectors)
     ParameterList p("Electron CVFEM-SG SubCV Centroid Current Density");
     p.set("Carrier Type", "Electron");
     p.set("Vector Name", n.field.elec_curr_density);
     p.set<RCP<const charon::Names> >("Names", m_names);
     p.set("Basis", hcurl_vol_cvfem);
     p.set("IR", cvfem_vol_ir);

     RCP< PHX::Evaluator<panzer::Traits> > op =
       rcp(new charon::SGCVFEM_CentroidCurrDens<EvalT,panzer::Traits>(p));
     fm.template registerEvaluator<EvalT>(op);
   }

   if(withAvaGen || withTrapSRH)  // need driving force for both ava. gen. and trap srh
   {
     // compute electron driving force at the centroids of the subcontrol volumes
     ParameterList p("Electron Centroid Driving Force");
     p.set("Carrier Type", "Electron");
     p.set("Driving Force", drForce);
     p.set("Scaling Parameters", scaleParams);
     p.set<RCP<const charon::Names> >("Names", m_names);
     p.set("HCurlBasis", hcurl_vol_cvfem);
     p.set("HGradBasis", hgrad_vol_cvfem);
     p.set("IR", cvfem_vol_ir);

     RCP< PHX::Evaluator<panzer::Traits> > op = 
       rcp(new charon::SGCVFEM_CentroidDriveForce<EvalT,panzer::Traits>(p));
     fm.template registerEvaluator<EvalT>(op);
   }

   // integrate the total source term at the centroid of the subcontrol volume
   // over subcv
   if (haveSource)
     {
       ParameterList p("Electron Source Residual");
       p.set("Residual Name", n.res.edensity+n.op.src);
       p.set("Value Name", n.field.total_recomb);
       p.set("Basis", hgrad_vol_cvfem);
       p.set("IR", cvfem_vol_ir);
       p.set("Multiplier", 1.0);
       p.set("WithInterpolation", false);

       RCP< PHX::Evaluator<panzer::Traits> > op =
         rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
       fm.template registerEvaluator<EvalT>(op);
     }
  }

  // Use a sum operator to form the overall residual for the equation
  {

    ParameterList p;
    p.set("Sum Name", n.res.edensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.edensity + n.op.conv_diff);

    if (this->buildTransientSupport())
      sum_names->push_back(n.res.edensity + n.op.trans);

    if (haveSource)
      sum_names->push_back(n.res.edensity + n.op.src);

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveElectron == "True")


 // ***************************************************************************
 // Construct Hole Continuity Equation based on the CVFEM-SG formulation
 // ***************************************************************************

 if (solveHole == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    // integrate dp_dt over subcv, \int_cv dp_dt dV, where dp_dt = \sum_BFi dpi_dt BFi
    // interpolate dp_dt at primary nodes to the center point of subcv and integrate over subcv
    ParameterList p("Hole Transient Residual");
    p.set("Residual Name", n.res.hdensity+n.op.trans);
    p.set("Value Name", n.dxdt.hdensity);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute hole CVFEM-SG edge current density at the midpoints of primary edges
   {
    ParameterList p("Hole CVFEM-SG Edge Current Density");
    p.set("Carrier Type", "Hole");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_EdgeCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // interpolate hole CVFEM-SG edge current density to the midpoints of subcv
   // edges(2D)/faces(3D) using edge basis functions (vectors)
   {
    ParameterList p("Hole CVFEM-SG SubCV Current Density");
    p.set("Carrier Type", "Hole");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", hcurl_side_cvfem);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_SubCVCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // integrate hole CVFEM-SG current density at the midpoints of subcv
   // edges(2D)/faces(3D) over the subcv boundary
   {
    ParameterList p("Hole Convection Diffusion Residual");
    p.set("Residual Name", n.res.hdensity+n.op.conv_diff);
    p.set("Flux Name", n.field.hole_curr_dens_cvedge);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", cvfem_side_ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   if(withAvaGen) 
   {
     // interpolate hole CVFEM-SG edge current density to centroids
     // of subcontrol volumes using edge basis functions (vectors)
     ParameterList p("Hole CVFEM-SG SubCV Centroid Current Density");
     p.set("Carrier Type", "Hole");
     p.set("Vector Name", n.field.hole_curr_density);
     p.set<RCP<const charon::Names> >("Names", m_names);
     p.set("Basis", hcurl_vol_cvfem);
     p.set("IR", cvfem_vol_ir);

     RCP< PHX::Evaluator<panzer::Traits> > op =
       rcp(new charon::SGCVFEM_CentroidCurrDens<EvalT,panzer::Traits>(p));
     fm.template registerEvaluator<EvalT>(op);
   }

   if(withAvaGen || withTrapSRH)  // need driving force for both ava. gen. and trap srh
   {
     // compute hole driving force at the centroids of the subcontrol volumes
     ParameterList p("Hole Centroid Driving Force");
     p.set("Carrier Type", "Hole");
     p.set("Driving Force", drForce);
     p.set("Scaling Parameters", scaleParams);
     p.set<RCP<const charon::Names> >("Names", m_names);
     p.set("HCurlBasis", hcurl_vol_cvfem);
     p.set("HGradBasis", hgrad_vol_cvfem);
     p.set("IR", cvfem_vol_ir);

     RCP< PHX::Evaluator<panzer::Traits> > op = 
       rcp(new charon::SGCVFEM_CentroidDriveForce<EvalT,panzer::Traits>(p));
     fm.template registerEvaluator<EvalT>(op);
   }

   // integrate the total source term at the centroid of the subcontrol volume over subcv
   if (haveSource)
   {
     ParameterList p("Hole Source Residual");
     p.set("Residual Name", n.res.hdensity+n.op.src);
     p.set("Value Name", n.field.total_recomb);
     p.set("Basis", hgrad_vol_cvfem);
     p.set("IR", cvfem_vol_ir);
     p.set("Multiplier", 1.0);
     p.set("WithInterpolation", false);

     RCP< PHX::Evaluator<panzer::Traits> > op =
       rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
     fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Use a sum operator to form the overall residual for the equation
  {
    ParameterList p;
    p.set("Sum Name", n.res.hdensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.hdensity + n.op.conv_diff);

    if (this->buildTransientSupport())
      sum_names->push_back(n.res.hdensity + n.op.trans);

    if (haveSource)
      sum_names->push_back(n.res.hdensity+ n.op.src);

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveHole == "True")

}

// ***********************************************************************

#endif
