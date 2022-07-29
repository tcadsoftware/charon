#ifndef CHARON_SOLVER_STEADYSTATE_CONSTRAINT_HPP
#define CHARON_SOLVER_STEADYSTATE_CONSTRAINT_HPP

#include "Teuchos_RCP.hpp"
#include "Charon_CurrentConstraintList.hpp"
#include <vector>
#include <ostream>

namespace Teuchos { class ParameterList; }
namespace Thyra {
  template<typename T> class ModelEvaluator;
}
namespace panzer_stk {
  template<typename T> class ModelEvaluatorFactory;
  class NOXObserverFactory;
}

namespace charon {

  /** Solve the problem using the LOCA Householder Bordering solver.

      \param[in] inputParameters Top level input ParameterList.
      \param[in] completePhysicsModel Physics ModelEvaluator wrapped in a charon::CurrentConstraintModelEvaluatorLOCA decorator.
      \param[in] currentConstraints CurrentConstraintList object with all active current constraints from BCs.
      \param[in] contact_sides Required for correctly creating the response names for current constraints.
      \param[in] serial_ostream Output stream that prints from one MPI process.
      \param[in] parallel_ostream Output stream that prints from every MPI process.
      \param[in] me_factory ModelEvalautor factory.
      \param[in] nox_observer_factory For creating and registering steady state solver observers.
  */
  void solve_steadystate_constraint_problem(const Teuchos::RCP<Teuchos::ParameterList>& input_parameters,
                                            const Teuchos::RCP<Thyra::ModelEvaluator<double>>& thyra_me,
                                            charon::CurrentConstraintList& current_constraints,
                                            const std::vector<std::pair<std::string, std::string>>& contact_sides, // for response names
                                            const Teuchos::RCP<std::ostream>& serial_ostream,
                                            const Teuchos::RCP<std::ostream>& parallel_ostream,
                                            const panzer_stk::ModelEvaluatorFactory<double>& me_factory,
                                            const panzer_stk::NOXObserverFactory& nox_observer_factory);

} // namespace charon

#endif
