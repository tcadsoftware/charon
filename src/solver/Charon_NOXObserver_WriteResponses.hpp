
#ifndef CHARON_NOXOBSERVER_WRITERESPONSES_HPP
#define CHARON_NOXOBSERVER_WRITERESPONSES_HPP

#include "Teuchos_FancyOStream.hpp"

#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_LineSearchBased.H"

namespace charon {

/**
 * \brief Class to output "responses" to the screen and a file at the
 * end of a Nonlinear iteration.
 *
 * This class optionally outputs the responses to the screen and a
 * file. Primarily this is used for outputting scalar electric currents
 * at the end of a steady-state solve or a LOCA continuation problem.
 */
class NOXObserver_WriteResponses : public NOX::Abstract::PrePostOperator {

public:

  /**
   * \brief Constructor.
   *
   * This is instantiated via this ctor from charon::NOXObserverFactory
   */
  NOXObserver_WriteResponses(std::vector<std::string> const& in_response_names,
                             Teuchos::RCP<Thyra::ModelEvaluator<double> > in_full_me,
                             bool writeToScreen,
                             bool writeToFile,
                             std::string const& outFilename,
                             bool isLOCASolver,
                             bool writeOnSolverFailure);

  /**
   * \brief This is the interface executed at the end of a nonlinear
   * solve.
   */
  void runPostSolve(NOX::Solver::Generic const& solver);

private:

  /**
   * \brief any specific steps for LOCA runs
   */
  Teuchos::RCP<NOX::Abstract::Vector const> get_loca_x(NOX::Solver::Generic const& solver);

  /**
   * \brief any specific steps for generic NOX runs
   */
  Teuchos::RCP<NOX::Abstract::Vector const> get_nox_x(NOX::Solver::Generic const& solver);

  /**
   * \brief for output to the screen
   */
  Teuchos::FancyOStream m_out;

  /**
   * \brief the response names
   */
  std::vector<std::string> response_names;

  /**
   * \brief a Thyra::ModelEvalutor that contains full set of constraints.
   */
  Teuchos::RCP<Thyra::ModelEvaluator<double> > full_me;

  bool write_to_screen;
  bool write_to_file;
  std::string out_filename;
  bool out_file_header_written;
  bool loca_solver;
  bool write_on_solver_failure;

};

}

#endif
