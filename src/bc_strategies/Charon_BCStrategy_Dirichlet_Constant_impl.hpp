
#ifndef CHARON_BCSTRATEGY_DIRICHLET_CONSTANT_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_CONSTANT_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PureBasis.hpp"

// Evaluators
#include "Panzer_Constant.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_Constant<EvalT>::
BCStrategy_Dirichlet_Constant(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Constant");
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_Constant<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  string dof_name = this->m_bc.equationSetName();

  // need the dof value to form the residual
  this->required_dof_names.push_back(dof_name);

  // unique residual name
  this->residual_name = "Residual_" + dof_name;

  // map residual to dof
  this->residual_to_dof_names_map[residual_name] = dof_name;

  // map residual to target field
  this->residual_to_target_field_map[residual_name] = "Constant_" + dof_name;

  // find the basis for this dof
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
    dofs.begin(); dof_it != dofs.end(); ++dof_it)
  {
    if (dof_it->first == dof_name) this->basis = dof_it->second;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
                     "Error the name \"" << this->m_bc.equationSetName()
                     << "\" is not a valid DOF for the boundary condition:\n"
                     << this->m_bc << "\n");
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_Constant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& /* pb */,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
                           const Teuchos::ParameterList& /* models */,
                           const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // provide a constant target value to map into residual
  {
    ParameterList p("BC Constant Dirichlet");
    p.set("Name", "Constant_" + this->m_bc.equationSetName());
    p.set("Data Layout", basis->functional);
    p.set("Value", this->m_bc.params()->template get<double>("Value"));

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Constant<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

#endif
