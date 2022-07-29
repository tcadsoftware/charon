
#ifndef CHARON_EQUATIONSET_NLPOISSON_IMPL_HPP
#define CHARON_EQUATIONSET_NLPOISSON_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Charon_Material_Properties.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Sum.hpp"
#include "Charon_Integrator_GradBasisDotGradDOF.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_NLPoissonSource.hpp"
#include "Charon_PotentialFlux.hpp"


// ***********************************************************************
template <typename EvalT>
charon::EquationSet_NLPoisson<EvalT>::
EquationSet_NLPoisson(const Teuchos::RCP<Teuchos::ParameterList>& params,
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
      "Fermi Dirac",
      "False",
      "Determine if users want to use the Fermi-Dirac statistics for the source term",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int integration_order = params->get<int>("Integration Order");

  UseFD = params->sublist("Options").get<std::string>("Fermi Dirac");

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
void charon::EquationSet_NLPoisson<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const charon::Names& n = *m_names;

  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_names->dof.phi);
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_names->dof.phi);

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // ***************************************************************************
  // Thermal Equilibrium Poisson Equation (all bias is at 0)
  // ***************************************************************************

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

  // Source Operator
  {
    // Compute the nonlinear Poisson source term @ IP
    {
      ParameterList p("Nonlinear Poisson Source");
      p.set("Source Name", n.field.nlprho);
      p.set("Data Layout", ir->dl_scalar);
      p.set("Scaling Parameters", scaleParams);
      p.set< RCP<const charon::Names> >("Names", m_names);

      p.set("Fermi Dirac",UseFD);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::NLPoissonSource<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the source term
    {
      using panzer::EvaluatorStyle;
      using panzer::Integrator_BasisTimesScalar;
      using panzer::Traits;
      using PHX::Evaluator;
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::CONTRIBUTES,
        n.res.phi, n.field.nlprho, *basis, *ir, -1));
      fm.template registerEvaluator<EvalT>(op);
    }
  }
}

// ***********************************************************************

#endif
