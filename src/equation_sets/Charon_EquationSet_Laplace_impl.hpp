
#ifndef CHARON_EQUATIONSET_LAPLACE_IMPL_HPP
#define CHARON_EQUATIONSET_LAPLACE_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Sum.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PotentialFlux.hpp"

#include "Charon_DisplacementCurrentDensity.hpp"
#include "Charon_PrevPotentialGrad.hpp"

// ***********************************************************************
template <typename EvalT>
charon::EquationSet_Laplace<EvalT>::
EquationSet_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
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

    // valid_parameters.sublist("Options");
    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");
    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "TID",
      "Off",
      "Determine if users want to add TID models in an insulator region",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    // determine if a user wants to include bulk fixed charges
    addFixCharge = false; 
    if (params->sublist("Options").get<std::string>("Fixed Charge") == "True")
      addFixCharge = true; 
    addTID = false; 
    if (params->sublist("Options").get<std::string>("TID") == "On")
      addTID = true; 
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

  // ********************
  // Assemble DOF names and Residual names
  // ********************
  m_names = Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix));

  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_names));

  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);
  if (this->buildTransientSupport())
    this->addDOFTimeDerivative(m_names->dof.phi,m_names->dxdt.phi);

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_Laplace<EvalT>::
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

  Teuchos::RCP<panzer::IntegrationRule> ir = (this->isFreqDom ? this->base_ir : this->getIntRuleForDOF(m_names->dof.phi));
  Teuchos::RCP<panzer::BasisIRLayout> basis = (this->isFreqDom ? this->base_basis : this->getBasisIRLayoutForDOF(m_names->dof.phi));

  // ********************
  // Laplace Equation
  // ********************
/*
  // Transient Operator
  if (this->buildTransientSupport())
  {
    using panzer::EvaluatorStyle;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Traits;
    using PHX::Evaluator;
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.phi, n.dxdt.phi, *basis, *ir));
    fm.template registerEvaluator<EvalT>(op);
  }
*/
  // Laplacian Operator
  // \int_{\Omega} \lambda2 \epsilon_r \grad_phi \cdot \grad_basis d\Omega
  {
    // Compute the potential flux: Lambda2*epsilon_r*grad_phi @ IP
    {
      ParameterList p("Potential Flux");
      p.set("Flux Name", n.field.phi_flux);
      p.set("Gradient Name", n.grad_dof.phi);
      p.set("IR", ir);
      p.set("Scaling Parameters", scaleParams);
      p.set< RCP<const charon::Names> >("Names", m_names);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::PotentialFlux<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the Laplacian term
    {
      ParameterList p("Laplacian Residual");
      p.set("Residual Name", n.res.phi);
      p.set("Flux Name", n.field.phi_flux);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::Integrator_GradBasisDotVector<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Add the fixed charges to the equation residual
  if (this->isFreqDom ? this->base_addFixCharge : this->addFixCharge)
  {
    using panzer::EvaluatorStyle;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Traits;
    using PHX::Evaluator;
    RCP<Evaluator<Traits>> op = rcp(new
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
      n.res.phi, n.field.fixed_charge, *basis, *ir, -1));
    fm.template registerEvaluator<EvalT>(op);
  }  
  
  // Add the TID trapped hole charges to the equation residual
  if (this->addTID && !this->isFreqDom)
  {
    {
      using panzer::EvaluatorStyle;
      using panzer::Integrator_BasisTimesScalar;
      using panzer::Traits;
      using PHX::Evaluator;
      RCP<Evaluator<Traits>> op = rcp(new
	Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
        n.res.phi, n.field.ins_htrappedcharge, *basis, *ir, -1));
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Computing displacement current density: Jd = -eps*d(grad(psi))/dt at IPs
  if (this->buildTransientSupport()) {
    {
      ParameterList p("Prev Potential Gradient");
      p.set("Current Name", m_names->field.grad_phi_prev);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Scaling Parameters", scaleParams);
      p.set("IR", ir);
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new charon::PrevPotentialGrad<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op); 
    }

    {
      ParameterList p("Displacement Current Density");
      p.set("Current Name", m_names->field.displacement_curr_density);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Scaling Parameters", scaleParams);
      p.set("IR", ir);

      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new charon::DisplacementCurrentDensity<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
  }
}

// ***********************************************************************

#endif
