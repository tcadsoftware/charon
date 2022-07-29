
#ifndef CHARON_EQUATIONSET_EFFPG_DDIONLATTICE_IMPL_HPP
#define CHARON_EQUATIONSET_EFFPG_DDIONLATTICE_IMPL_HPP

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
#include "Charon_EFFPG_DDIonLattice_CurrentDensity.hpp"
#include "Charon_DDLattice_ElectricField.hpp"
#include "Charon_DDLattice_HeatGeneration.hpp"
#include "Charon_DDLattice_VacuumPotential.hpp"




// ***********************************************************************
template <typename EvalT>
charon::EquationSet_EFFPG_DDIonLattice<EvalT>::
EquationSet_EFFPG_DDIonLattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
                               const int& default_integration_order,
                               const panzer::CellData& cell_data,
                               const Teuchos::RCP<panzer::GlobalData>& global_data,
                               const bool build_transient_support) :
  charon::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support)
{
  // Options
  {
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);
    valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
    valid_parameters.set("Prefix","","Prefix for using multiple instatiations of their equation set");
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
      "Solve Ion",
      "True",
      "Determines if users want to solve the Ion Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    opt.set<int>("Ion Charge", 1, "Integer number of ion charge");

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
      "Band Gap Narrowing",
      "Off",
      "Determines if users want to turn on BGN effect due to heavy doping",
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

    params->validateParametersAndSetDefaults(valid_parameters);

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");
    solveIon = params->sublist("Options").get<std::string>("Solve Ion");
    // Get the value of Ion Charge
    ion_charge = params->sublist("Options").get<int>("Ion Charge");

    // Must solve at least one continuity equation
    if ((solveElectron == "False") && (solveHole == "False") && (solveIon == "False"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: must solve at least one"
           << " continuity equation !" << "\n");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH");
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");

    // Electric Field Model = Potential Gradient by default, and can be expanded as needed
    field_model = params->sublist("Options").get<std::string>("Electric Field Model");

    // Must solve both continuity equations when any recomb./gen. model is on.
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || (ava_gen == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True") &&(solveIon == "True"));
      haveSource = true;
    }
    else
      haveSource = false;  // NO recomb. model

    // Set the flag for BGN
    const std::string& bgn = params->sublist("Options").get<std::string>("Band Gap Narrowing");
    if (bgn == "On")
      haveBGN = true;
    else
      haveBGN = false;

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
  std::string discfields = params->isParameter("Discontinuous Fields") ? params->get<std::string>("Discontinuous Fields") : "";
  std::string discsuffix = params->isParameter("Discontinuous Suffix") ? params->get<std::string>("Discontinuous Suffix") : "";
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

  // Solve the Ion equation by default and users can disable it
  if (solveIon == "True")
  {
    this->addDOF(m_names->dof.iondensity,basis_type,basis_order,integration_order,m_names->res.iondensity);
    this->addDOFGrad(m_names->dof.iondensity,m_names->grad_dof.iondensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.iondensity,m_names->dxdt.iondensity);
  }

  // Always solve the lattice temperature equation
  this->addDOF(m_names->dof.latt_temp,basis_type,basis_order,integration_order,m_names->res.latt_temp);
  this->addDOFGrad(m_names->dof.latt_temp,m_names->grad_dof.latt_temp);
  if (this->buildTransientSupport())
    this->addDOFTimeDerivative(m_names->dof.latt_temp,m_names->dxdt.latt_temp);

  this->addClosureModel(model_id);

  this->setupDOFs();

}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_EFFPG_DDIonLattice<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
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
//  This flag isn't needed, since it will be set to default
//  bool includeSoret = true;  // set "Include Soret Effect"

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // For now, assume that all dofs use the same basis
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_names->dof.phi);
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_names->dof.phi);

  // ***************************************************************************
  // Construct Poisson Equation (solved with continuity equations)
  // The same code as in Charon_EquationSet_DDIonLattice_impl.hpp
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
    // Compute the Poisson source term: (p-n+i+Nd-Na) @ IP
    {
      ParameterList p("Poisson Source");
      p.set("Source Name", n.field.psrc);
      p.set("Data Layout", ir->dl_scalar);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Solve Electron", solveElectron);
      p.set("Solve Hole", solveHole);
      p.set<bool>("Solve Ion", true);  // always true for DDIon
      p.set<int>("Ion Charge", ion_charge);
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
  // Construct Electron Continuity Equation
  // ***************************************************************************

 if (solveElectron == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.edensity, n.dxdt.edensity, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute electron EFFPG DDIonLattice current density
   {
    ParameterList p("Electron EFFPG DDIonLattice Current Density");
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

   // integrate convection(/advection/drift)-diffusion term (current density)
   {
    ParameterList p("Electron Convection Diffusion Residual");
    p.set("Residual Name", n.res.edensity);
    p.set("Flux Name", n.field.elec_curr_density);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   // integrate the total source term (exclude avalanche gen.)
   if (haveSource)
   {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.edensity, n.field.total_recomb, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
   }
  }
 }  // end of if (solveElectron == "True")


  // ***************************************************************************
  // Construct Hole Continuity Equation
  // ***************************************************************************

 if (solveHole == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.hdensity, n.dxdt.hdensity, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute hole EFFPG DDIonLattice current density
   {
    ParameterList p("Hole EFFPG DDIonLattice Current Density");
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

   // integrate convection(advection/drift)-diffusion term (current density)
   {
    ParameterList p("Hole Convection Diffusion Residual");
    p.set("Residual Name", n.res.hdensity);
    p.set("Flux Name", n.field.hole_curr_density);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   // integrate the total source term (exclude avalanche gen.)
   if (haveSource)
   {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.hdensity, n.field.total_recomb, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
   }
  }
 }  // end of if (solveHole == "True")


  // ***************************************************************************
  // Construct Ion Continuity Equation
  // ***************************************************************************

 if (solveIon == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.iondensity, n.dxdt.iondensity, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute ion EFFPG DDIonLattice current density
   {
    ParameterList p("Ion EFFPG DDIonLattice Current Density");
    p.set("Carrier Type", "Ion");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set<bool>("Temperature Gradient", true);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::EFFPG_DDIonLattice_CurrentDensity<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // integrate convection(/advection/drift)-diffusion term (current density)
   {
    ParameterList p("Ion Convection Diffusion Residual");
    p.set("Residual Name", n.res.iondensity);
    p.set("Flux Name", n.field.ion_curr_density);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }
  }
 }  // end of if (solveIon == "True")


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
   // obtain electron electric field at IPs, used in computing the heat generation H
   {
    ParameterList p("Electron Electric Field");
    p.set("Carrier Type", "Electron");
    p.set("Band Gap Narrowing", haveBGN);
    p.set("Electric Field Model", field_model);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // obtain hole electric field at IPs, used in computing the heat generation H
   {
    ParameterList p("Hole Electric Field");
    p.set("Carrier Type", "Hole");
    p.set("Band Gap Narrowing", haveBGN);
    p.set("Electric Field Model", field_model);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // obtain ion electric field at IPs, used in computing the heat generation H
   {
    ParameterList p("Ion Electric Field");
    p.set("Carrier Type", "Ion");
    p.set("Band Gap Narrowing", false);
    p.set("Electric Field Model", field_model);
    p.set< RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLattice_ElectricField<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
   }

   // Obtain the heat generation
   {
    bool enableElectron = false;
    bool enableHole = false;
    bool enableIon = false;
    if (solveElectron == "True") enableElectron = true;
    if (solveHole == "True") enableHole = true;
    if (solveIon == "True") enableIon = true;

    ParameterList p("Heat Generation");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("IR", ir);
    p.set<bool>("Solve Electron", enableElectron);
    p.set<bool>("Solve Hole", enableHole);
    p.set<bool>("Solve Ion", enableIon);
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
