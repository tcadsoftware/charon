
#ifndef CHARON_EQUATIONSET_LATTICE_IMPL_HPP
#define CHARON_EQUATIONSET_LATTICE_IMPL_HPP

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


// SUPG-FEM formulation of the lattice-temperature heat equation

// ***********************************************************************
template <typename EvalT>
charon::EquationSet_Lattice<EvalT>::
EquationSet_Lattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
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

    valid_parameters.set("Model ID","","Closure model id associated with this equaiton set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");

    Teuchos::setStringToIntegralParameter<int>(
      "Solve DD",
      "False",
      "Determine if the Poisson+DD eqns are solved together with the lattice eqn",
      Teuchos::tuple<std::string>("False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Heat Generation",
      "Analytic",
      "Determine the type of heat generation",
      Teuchos::tuple<std::string>("Analytic"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);
  }

  solveDD = params->sublist("Options").get<std::string>("Solve DD");
  heatGenType = params->sublist("Options").get<std::string>("Heat Generation");

  if ( (solveDD == "False") && (heatGenType == "Joule Heating") )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
      "Error: Joule Heating can be used only when the Poisson+DD eqns are solved with the lattice eqn !");

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
void charon::EquationSet_Lattice<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Traits;
  using PHX::Evaluator;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const charon::Names& n = *m_names;

  // For now, assume that all dofs use the same basis
  RCP<IntegrationRule> ir = this->getIntRuleForDOF(m_names->dof.latt_temp);
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_names->dof.latt_temp);


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
   // If heatGenType = "Analytic", the heat generation is computed by charon::Analytic_HeatGeneration
   // instantiated in the Charon_ClosureModel_Factory_impl.hpp file.

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
