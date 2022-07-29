
#ifndef CHARON_EQUATIONSET_DDLATTICE_IMPL_HPP
#define CHARON_EQUATIONSET_DDLATTICE_IMPL_HPP

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
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Sum.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PotentialFlux.hpp"
#include "Charon_PoissonSource.hpp"

// FEM SUPG-related evaluators
#include "Charon_SUPG_Peclet.hpp"
#include "Charon_PDE_Residual_DD.hpp"
#include "Charon_SUPG_Tau_Linear.hpp"
#include "Charon_SUPG_Tau_Tanh.hpp"

// FEM vector field evaluators
#include "Charon_DDLattice_ElectricField.hpp"
#include "Charon_DDLattice_CurrentDensity.hpp"
#include "Charon_DDLattice_HeatGeneration.hpp"
#include "Charon_DDLattice_VacuumPotential.hpp"
#include "Charon_FEM_Velocity.hpp"
#include "Charon_Effective_Velocity.hpp"
#include "Charon_FEM_ElectricField.hpp"

// Include the file to compute current density using the EFFPG stabilization
#include "Charon_EFFPG_DDIonLattice_CurrentDensity.hpp"


// Implement the coupled DDLattice equations (Poisson+Continuity+Heating)

// ***********************************************************************
template <typename EvalT>
charon::EquationSet_DDLattice<EvalT>::
EquationSet_DDLattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
      "Determine if users want to solve the Electron Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Hole",
      "True",
      "Determine if users want to solve the Hole Continuity equation",
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
      "Diffusion",  // default
      "Use the Diffusion method for FD implementation",
      Teuchos::tuple<std::string>("Diffusion"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determine if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determine if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determine if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determine if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SUPG Stabilization",
      "On",
      "Enable or disable SUPG stabilization in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_E",
      "Tanh",
      "Determine the Tau_E model in the Electron equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_H",
      "Tanh",
      "Determines the Tau_H model in the Hole equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Length Scale",
      "Stream",
      "Determine the Length Scale model in the SUPG stabilization",
      Teuchos::tuple<std::string>("Stream","Shakib"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Source Term",
      "Off",
      "Determine if users want to add the source term in the SUPG stabilization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determine if users want to turn on the BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Electric Field Model",
      "Potential Gradient",
      "Determine the model for the electric field calculation",
      Teuchos::tuple<std::string>("Potential Gradient"),
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
      "Discretization Method",
      "FEM-SUPG",
      "Determine the discretization method",
      Teuchos::tuple<std::string>("FEM-SUPG", "FEM-EFFPG", "CVFEM-SG"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Driving Force",
      "GradPotential",
      "Determines driving force for different models in the continuity equation(s)",
      Teuchos::tuple<std::string>("EffectiveField","GradQuasiFermi","GradPotential"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Optical Generation",
      "Off",
      "Determine if users want to turn on the optical generation",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    params->validateParametersAndSetDefaults(valid_parameters);
  }

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH");
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");
    const std::string& opt_gen = params->sublist("Options").get<std::string>("Optical Generation"); 

    // Must solve both continuity equations when any recombination model is on.
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || 
        (ava_gen == "On") || (opt_gen == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
      haveSource = true;
    }
    else
      haveSource = false;  // NO recomb. model

    supg_stab = params->sublist("Options").get<std::string>("SUPG Stabilization");
    tau_e_type = params->sublist("Options").get<std::string>("Tau_E");
    tau_h_type = params->sublist("Options").get<std::string>("Tau_H");
    ls_type = params->sublist("Options").get<std::string>("Length Scale");

    const std::string& source_stab = params->sublist("Options").get<std::string>("Add Source Term");
    if (source_stab == "On")
    {
      if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || 
          (ava_gen == "On") || (opt_gen == "On"))
        add_source_stab = true;
      else
        add_source_stab = false;
    }
    else
      add_source_stab = false;

    // Electric Field Model = Potential Gradient by default, and can be expanded as needed
    field_model = params->sublist("Options").get<std::string>("Electric Field Model");

    // Set the flag for BGN
    const std::string& bgn = params->sublist("Options").get<std::string>("Band Gap Narrowing");
    if (bgn == "On")
      haveBGN = true;
    else
      haveBGN = false;

    // Determine the discretization method
    discMethod = params->sublist("Options").get<std::string>("Discretization Method");

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

  // Always solve the Poisson equation (no transients required)
  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);

  // Solve the electron continuity equation by default and users can disable it
  if (solveElectron == "True")
  {
    this->addDOF(m_names->dof.edensity,basis_type,basis_order,integration_order,m_names->res.edensity);
    this->addDOFGrad(m_names->dof.edensity,m_names->grad_dof.edensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.edensity,m_names->dxdt.edensity);
  }

  // Solve the hole continuity equation by default and users can disable it
  if (solveHole == "True")
  {
    this->addDOF(m_names->dof.hdensity,basis_type,basis_order,integration_order,m_names->res.hdensity);
    this->addDOFGrad(m_names->dof.hdensity,m_names->grad_dof.hdensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.hdensity,m_names->dxdt.hdensity);
  }

  // Always solve the lattice temperature equation
  this->addDOF(m_names->dof.latt_temp,basis_type,basis_order,integration_order,m_names->res.latt_temp);
  this->addDOFGrad(m_names->dof.latt_temp,m_names->grad_dof.latt_temp);
  if (this->buildTransientSupport())
    this->addDOFTimeDerivative(m_names->dof.latt_temp,m_names->dxdt.latt_temp);

  this->addClosureModel(model_id);
  this->setupDOFs();

  // Create an EquationSet_DDIonLattice object to call its computeSUPGStabResidual(.) member function
  this->eqnSetDDIonLatt = Teuchos::rcp(new charon::EquationSet_DDIonLattice<EvalT>(params, default_integration_order, cell_data, global_data,build_transient_support,"dummy"));

}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_DDLattice<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary&  fl,
                                      const Teuchos::ParameterList& user_data) const
{
  using panzer::EvaluatorStyle;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Traits;
  using PHX::Evaluator;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const charon::Names& n = *m_names;
  bool includeSoret = true;

  // Get scaling parameters
  RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // For now, assume that all dofs use the same basis
  RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_names->dof.phi);
  RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_names->dof.phi);


  // ***************************************************************************
  // Construct the Poisson equation (solved with the continuity equations)
  // ***************************************************************************
/*
  // Compute the vacuum potential
  {
    ParameterList p("Vacuum Potential");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Data Layout", ir->dl_scalar);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_VacuumPotential<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }
  {
    ParameterList p("Vacuum Potential");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Data Layout", basis->functional);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_VacuumPotential<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Compute the gradient of vacuum potential at IPs
  {
    ParameterList p("Vacuum Potential Gradient");
    p.set("Name", n.field.vac_pot);
    p.set("Gradient Name", n.field.grad_vac_pot);
    p.set("Basis", basis);
    p.set("IR", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }
*/

  // Laplacian Operator: \int_{\Omega} \lambda2 \epsilon_r \grad_vacpot \cdot \grad_basis d\Omega
  {
    ParameterList p("Laplacian Residual");
    p.set("Residual Name", n.res.phi);
    // p.set("Flux Name", n.field.grad_vac_pot);
    p.set("Flux Name", n.grad_dof.phi);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", scaleParams->scale_params.Lambda2);

    Teuchos::RCP<std::vector<std::string> > fms = Teuchos::rcp(new std::vector<std::string>);
    fms->push_back(n.field.rel_perm);
    p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Source Operator
  {
    // Compute the Poisson source term: (p-n+Nd-Na) @ IP
    {
      ParameterList p("Poisson Source");
      p.set("Source Name", n.field.psrc);
      p.set("Data Layout", ir->dl_scalar);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Solve Electron", solveElectron);
      p.set("Solve Hole", solveHole);
      p.set("Scaling Parameters", scaleParams);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::PoissonSource<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the source term
    {
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
        n.res.phi, n.field.psrc, *basis, *ir, -1));
      fm.template registerEvaluator<EvalT>(op);
    }
  }


  // ***************************************************************************
  // Construct the electron continuity equation
  // ***************************************************************************

 if (solveElectron == "True")
 {
  // Transient Operator
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.edensity, n.dxdt.edensity, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Construct the residual of current density dot grad_basis
  {
   // Compute the electric field whose expression depends on the model of choice.
   // Currently it is trivially set to -\grad_phi, but can be expanded to include other models.
    /*   {
    ParameterList p("Electron Electric Field");
    p.set("Carrier Type", "Electron");
    p.set("Electric Field Model", field_model);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Band Gap Narrowing", haveBGN);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
    }*/

   // Computing the electron effective electric field: Fn,eff at IPs,
   // valid for BGN = On and Off (dEg=0 when BGN = Off).
   {
    ParameterList p("Electron Electric Field");
    p.set("Carrier Type", "Electron");
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // Current density calculation depends on the discretization method
   if (discMethod == "FEM-SUPG")
   {
    // Computing electron current density at IPs: Jn = n*mun*Fn + Dn*\grad_n + n*mun*\grad_T
    ParameterList p("Electron Current Density");
    p.set("Carrier Type", "Electron");
    p.set("Current Name", m_names->field.elec_curr_density);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_CurrentDensity<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   else if (discMethod == "FEM-EFFPG")
   {
    // Compute electron current density using the EFFPG stabilization
    ParameterList p("Electron EFFPG Current Density");
    p.set("Carrier Type", "Electron");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set<bool>("Temperature Gradient", true);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::EFFPG_DDIonLattice_CurrentDensity<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // The current density dot grad_basis residual
   {
    ParameterList p("Electron Convection Diffusion Residual");
    p.set("Residual Name", n.res.edensity+n.op.conv_diff);
    p.set("Flux Name", n.field.elec_curr_density);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

  } // end of the residual of current density dot grad_basis

  // Source Operator: Integrate the total source term
  if (haveSource)
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.edensity, n.field.total_recomb, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  if ((discMethod == "FEM-SUPG") && (supg_stab == "On"))
    eqnSetDDIonLatt->computeSUPGStabResidual(fm, fl, user_data, m_names, ir, basis, ls_type, tau_e_type, add_source_stab, "Electron", 0, includeSoret);

  // Use a sum operator to form the overall residual for the equation
  {
    ParameterList p;
    p.set("Sum Name", n.res.edensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.edensity + n.op.conv_diff);

    if ((discMethod == "FEM-SUPG") && (supg_stab == "On"))
    {
      sum_names->push_back(n.res.edensity + n.op.supg_conv);
      if (add_source_stab) sum_names->push_back(n.res.edensity + n.op.supg_src);
    }

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveElectron == "True")


  // ***************************************************************************
  // Construct the hole continuity equation
  // ***************************************************************************

 if (solveHole == "True")
 {
  // Transient Operator
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.hdensity, n.dxdt.hdensity, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Construct the convection-diffusion residual
  {
   // Compute the electric field whose expression depends on the model of choice.
   // Currently it is trivially set to -\grad_phi, but can be expanded to include other models.
    /*   {
    ParameterList p("Hole Electric Field");
    p.set("Carrier Type", "Hole");
    p.set("Electric Field Model", field_model);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Band Gap Narrowing", haveBGN);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
    }*/

   // Computing the hole effective electric field: Fp,eff at IPs,
   // valid for BGN = On and Off (dEg=0 when BGN = Off).
   {
    ParameterList p("Hole Electric Field");
    p.set("Carrier Type", "Hole");
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // Current density calculation depends on the discretization method
   if (discMethod == "FEM-SUPG")
   {
    // Computing hole current density at IPs: Jp = p*mup*Fp - Dp*\grad_p - p*mup*\grad_T
    ParameterList p("Hole Current Density");
    p.set("Carrier Type", "Hole");
    p.set("Current Name", m_names->field.hole_curr_density);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_CurrentDensity<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   else if (discMethod == "FEM-EFFPG")
   {
    // compute hole current density using the EFFPG stabilization
    ParameterList p("Hole EFFPG Current Density");
    p.set("Carrier Type", "Hole");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set<bool>("Temperature Gradient", true);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::EFFPG_DDIonLattice_CurrentDensity<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // The convection-diffusion residual
   {
    ParameterList p("Hole Convection Diffusion Residual");
    p.set("Residual Name", n.res.hdensity + n.op.conv_diff);
    p.set("Flux Name", n.field.hole_curr_density);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

  }  // end of constructing convection-diffusion residual

  // Source Operator: Integrate the total source term
  if (haveSource)
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.hdensity, n.field.total_recomb, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  if ((discMethod == "FEM-SUPG") && (supg_stab == "On"))
    eqnSetDDIonLatt->computeSUPGStabResidual(fm, fl, user_data, m_names, ir, basis, ls_type, tau_h_type, add_source_stab, "Hole", 0, includeSoret);

  // Use a sum operator to form the overall residual for the equation
  {
    ParameterList p;
    p.set("Sum Name", n.res.hdensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.hdensity + n.op.conv_diff);

    if ((discMethod == "FEM-SUPG") && (supg_stab == "On"))
    {
      sum_names->push_back(n.res.hdensity + n.op.supg_conv);
      if (add_source_stab) sum_names->push_back(n.res.hdensity + n.op.supg_src);
    }

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveHole == "True")


  // ***************************************************************************
  // Construct the lattice temperature equation
  // ***************************************************************************
  // Transient Operator: \int_{\Omega} cL * \frac {\partial_T} {\partial_t} * basis d\Omega
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.latt_temp, n.dxdt.latt_temp, *basis, *ir, 1, {n.field.heat_cap}));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Laplacian Operator: \int_{\Omega} kL * \grad_T \cdot \grad_basis d\Omega
  {
    ParameterList p("Lattice Temperature Laplacian Residual");
    p.set("Residual Name", n.res.latt_temp);
    p.set("Flux Name", n.grad_dof.latt_temp);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    Teuchos::RCP<std::vector<std::string> > fms = Teuchos::rcp(new std::vector<std::string>);
    fms->push_back(n.field.kappa);
    p.set< Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers",fms);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Source Operator: \int_{\Omega} H * basis d\Omega
  {
   // Obtain the heat generation
   {
    bool enableElectron = false;
    bool enableHole = false;
    if (solveElectron == "True") enableElectron = true;
    if (solveHole == "True") enableHole = true;

    ParameterList p("Heat Generation");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set<bool>("Solve Electron", enableElectron);
    p.set<bool>("Solve Hole", enableHole);
    p.set<bool>("Solve Ion", false);
    p.set<bool>("Have Source", haveSource);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_HeatGeneration<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // Integrate the source term
   {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.latt_temp, n.field.heat_gen, *basis, *ir, -1));
    fm.template registerEvaluator<EvalT>(op);
   }
  }
}

// ***********************************************************************

#endif
