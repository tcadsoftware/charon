
#include <fstream>

// Charon
#include "Charon_CurrentConstraintList.hpp"
#include "Charon_CurrentConstraintModelEvaluator.hpp"
#include "Charon_NOXObserver_WriteResponses.hpp"

// LOCA
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_Thyra_Group.H"

// NOX
#include "NOX_Abstract_PrePostOperator.H"

// Piro
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"

// Teuchos
#include "Teuchos_dyn_cast.hpp"

// Thyra
#include "Thyra_DefaultProductVector_def.hpp"
#include "Thyra_SpmdVectorBase.hpp"

// System
#include <fstream>

namespace charon {

//***********************************************************************************************************
// ctor
//***********************************************************************************************************
NOXObserver_WriteResponses::NOXObserver_WriteResponses(std::vector<std::string> const& in_response_names,
                                                       Teuchos::RCP<Thyra::ModelEvaluator<double> > in_full_me,
                                                       bool writeToScreen,
                                                       bool writeToFile,
                                                       std::string const& outFilename,
                                                       bool isLOCASolver,
                                                       bool writeOnSolverFailure) :
  m_out(Teuchos::rcpFromRef(std::cout)),
  response_names(in_response_names),
  full_me(in_full_me),
  write_to_screen(writeToScreen),
  write_to_file(writeToFile),
  out_filename(outFilename),
  out_file_header_written(false),
  loca_solver(isLOCASolver),
  write_on_solver_failure(writeOnSolverFailure)
{

  // If this is a LOCA run and the responses are going to be output to a
  // file then set things up here.
  if (loca_solver && write_to_file)
  {
    if (out_filename == "")
      out_filename = "currents-loca.dat"; // default output file

  }
  if (write_to_screen)
    m_out.setOutputToRootOnly(0);

}

//***********************************************************************************************************
// NOX interface implementation to be run at the end of a nonlinear
// solve.
//***********************************************************************************************************
void
NOXObserver_WriteResponses::
runPostSolve(
  const NOX::Solver::Generic& solver)
{
  using charon::CurrentConstraintList;
  using LOCA::MultiContinuation::AbstractStrategy;
  using std::endl;
  using std::max;
  using std::setw;
  using std::size_t;
  using std::string;
  using std::vector;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::ProductVectorBase;
  using Thyra::SpmdVectorBase;
  using Thyra::VectorBase;
  using CCME = charon::CurrentConstraintModelEvaluator<double>;

  // Unless overridden, return without writing anything in the case of a
  // solver failure.
  //
  // GLH: Don't ask me why "getStatus()" isn't declared const so that I
  // could avoid the "const_cast" below.
  if (const_cast<NOX::Solver::Generic&>(solver).getStatus() == NOX::StatusTest::Failed && !write_on_solver_failure)
  {
    return;
  }

  RCP<NOX::Abstract::Vector const> x;
  if (loca_solver)
    x = get_loca_x(solver);
  else
    x = get_nox_x(solver);

  NOX::Thyra::Vector const& nThX = Teuchos::dyn_cast<NOX::Thyra::Vector const>(*x);
  RCP<VectorBase<double> const> const& thX = nThX.getThyraRCPVector();

  // Physics model-evaluator based evaluation of responses
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = full_me->createInArgs();
  inArgs.set_x(thX);

  // Allocate vectors for each of the responses
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = full_me->createOutArgs();
  for(int i=0; i < outArgs.Ng(); ++i)
    outArgs.set_g(i,Thyra::createMember(*(full_me->get_g_space(i))));

  // Evaluate the model (computing the responses)
  full_me->evalModel(inArgs, outArgs);

  // Output the responses.
  if ((write_to_screen) or (write_to_file))
  {
    int Ng = outArgs.Ng();
    size_t precision(8), value_width(precision + 7);

    m_out << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
    if (write_to_screen)
      m_out << "Responses after nonlinear solve:" << endl;

    // If this is a LOCA solve output the associated value of the
    // parameter, a voltage at the moment
    if (loca_solver)
    {
      auto locaStrategy = rcp_dynamic_cast<const AbstractStrategy>(
        solver.getSolutionGroupPtr());
      string continuationParameterName(
        locaStrategy->getContinuationParameterName());
      double continuationParameterValue(
        locaStrategy->getContinuationParameter());

      if (write_to_screen)
        m_out << "  Continuation parameter: " << continuationParameterName
              << "; Value: " << continuationParameterValue << endl;

      // If the output is going to a file
      if (write_to_file)
      {

        Teuchos::RCP<std::ofstream> fileOutput = Teuchos::rcp(new std::ofstream());
        fileOutput->open(out_filename.c_str(), std::ofstream::out | std::ofstream::app);

        Teuchos::RCP<Teuchos::FancyOStream> f_out = Teuchos::rcp(new Teuchos::FancyOStream(fileOutput));

        // Only output to processor 0
        f_out->setOutputToRootOnly(0);
        f_out->setTabIndentStr(" |   ");

        *f_out << std::scientific
               << std::showpoint
               << std::setprecision(precision)
               << std::right;

        // Determine the widths of the columns to be printed below.
        vector<size_t> column(Ng + 1);
        column[0] = max(continuationParameterName.size(),
          value_width);
        for (int i(0); i < Ng; ++i)
          column[i + 1] = max(response_names[i].size(), value_width);

        // This just writes the header out, voltage and currents, the first
        // time through.
        if (!out_file_header_written)
        {
          out_file_header_written = true;

          *f_out << setw(column[0]) << continuationParameterName << " ";
          for (int i(0); i < Ng - 1; ++i)
            *f_out << setw(column[i + 1]) << response_names[i] << " ";
          *f_out << setw(column[Ng]) << response_names[Ng - 1] << endl;
        }

        // The parameter (voltage) value
        *f_out << setw(column[0]) << continuationParameterValue << " ";

        // The responses (currents)
        for (int i=0; i < Ng-1; ++i)
        {
          RCP<VectorBase<double> > g = outArgs.get_g(i);

          // For total correctness this should loop over the dimensions
          // of the parameters but we only calculate a single current
          // value right now. See the screen output below for an example
          // of how we might write out values for more than one
          // dimension here.
          *f_out << setw(column[i + 1]) << Thyra::get_ele(*g, 0) << " ";
        }
        Teuchos::RCP<Thyra::VectorBase<double> > g = outArgs.get_g(Ng-1);
        *f_out << setw(column[Ng]) << Thyra::get_ele(*g, 0) << std::endl;

        fileOutput->close();
      }

    }

    // Output the responses (currents) to screen.
    if (write_to_screen)
    {
      // Loop over the responses...
      for (int i(0); i < Ng; ++i)
      {
        RCP<VectorBase<double>> g = outArgs.get_g(i);
        m_out << "    " << response_names[i] << " = ";
        for (Thyra::Ordinal k(0); k < g->space()->dim(); ++k)
          m_out << Thyra::get_ele(*g, k) << " ";
        m_out << endl;
      } // end loop over the responses

      // Check to see if we have any current constraints.
      RCP<const CCME> ccme = rcp_dynamic_cast<const CCME>(full_me);
      if (not ccme.is_null())
      {
        RCP<const VectorBase<double>> volt =
          rcp_dynamic_cast<const ProductVectorBase<double>>(thX, true)->
          getVectorBlock(1);
        ArrayRCP<const double> vData;
        rcp_dynamic_cast<const SpmdVectorBase<double>>(volt, true)->
          getLocalData(ptrFromRef(vData));
        m_out << "Voltage values:" << endl;
        const CurrentConstraintList& constraints = ccme->constraints();
        for (int i(0); i < vData.size(); ++i)
          m_out << "  " << constraints[i]->sidesetId()
                << constraints[i]->type() << "Voltage = " << vData[i] << endl;
      } // end if we have any current constraints
    } // end if (write_to_screen)
  } // end if ((write_to_screen) or (write_to_file))
} // end of runPostSolve()

//***********************************************************************************************************
// There are slight differences between a NOX-only run and LOCA one. This is for NOX-specific things
//***********************************************************************************************************
Teuchos::RCP<NOX::Abstract::Vector const>
NOXObserver_WriteResponses::get_nox_x(NOX::Solver::Generic const& solver)
{

  NOX::Abstract::Group const& grp = solver.getSolutionGroup();
  Teuchos::RCP<NOX::Abstract::Vector const> x = grp.getXPtr();

  return x;
}

//***********************************************************************************************************
// There are slight differences between a NOX-only run and LOCA one. This is for LOCA-specific things
//***********************************************************************************************************
Teuchos::RCP<NOX::Abstract::Vector const>
NOXObserver_WriteResponses::get_loca_x(NOX::Solver::Generic const& solver)
{

  // Get the solution vector for the original, non-constraint,
  // non-continuation problem when LOCA is being used
  NOX::Abstract::Group const& grp = solver.getSolutionGroup();
  NOX::Abstract::Vector const& x = grp.getX();
  LOCA::MultiContinuation::ExtendedVector const& l_ext_x =
    Teuchos::dyn_cast<LOCA::MultiContinuation::ExtendedVector const>(x);

  Teuchos::RCP<NOX::Abstract::Vector const> n_x = l_ext_x.getXVec();

  return n_x;

}

}
