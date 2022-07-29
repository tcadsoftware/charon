
#include <fstream>

// Charon
#include "Charon_CurrentConstraintList.hpp"
#include "Charon_CurrentConstraintModelEvaluator.hpp"
#include "Charon_CurrentConstraintModelEvaluatorLOCA.hpp"
#include "Charon_NOXObserver_WriteResponses.hpp"

// LOCA
#include "LOCA_MultiContinuation_ExtendedVector.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
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
						       Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                                                       bool isLOCASolver,
                                                       bool writeOnSolverFailure) :
  m_out(Teuchos::rcpFromRef(std::cout)),
  response_names(in_response_names),
  full_me(in_full_me),
  write_to_screen(writeToScreen),
  write_to_file(writeToFile),
  out_filename(outFilename),
  parameterLibrary_(parameterLibrary),
  out_file_header_written(false),
  loca_solver(isLOCASolver),
  write_on_solver_failure(writeOnSolverFailure),
  loca_householder_solver(false)
{
  // Check for new LOCA Householder solver
  {
    auto check = Teuchos::rcp_dynamic_cast<charon::CurrentConstraintModelEvaluatorLOCA<double>>(full_me,false);
    if (nonnull(check))
      loca_householder_solver = true;
  }

  // If this is a LOCA run and the responses are going to be output to a
  // file then set things up here.
  if (loca_solver && write_to_file)
  {
    if (out_filename == "")
      out_filename = "currents-loca.dat"; // default output file

  }
  if (write_to_screen)
    m_out.setOutputToRootOnly(0);

  //Get contact voltage names ; create permutation for output
  pPE = Teuchos::rcp(new panzerParameterExtractor(parameterLibrary_));
  contactVoltages_ = pPE->get_VoltageParameters();
  std::string  cpName = pPE->getContinuationParameterName();
  cvPermutation.resize(contactVoltages_.size());
  cvNames.resize(contactVoltages_.size());
  int cvCounter=0;
  std::map<std::string,double>::iterator pLit = contactVoltages_.begin();
  while (pLit != contactVoltages_.end())
    {
      cvPermutation[cvCounter] = cvCounter;
      cvNames[cvCounter] = pLit->first;
      if (pLit->first == cpName) // perform permute
	{
	  cvPermutation[0] = cvCounter;
	  cvPermutation[cvCounter] = 0;
	  std::string tempName = cvNames[0];
	  cvNames[0] = pLit->first;
	  cvNames[cvCounter] = tempName;
	}
      ++cvCounter;
      ++pLit;
    }

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
  using CCME_LOCA = charon::CurrentConstraintModelEvaluatorLOCA<double>;

  // Unless overridden, return without writing anything in the case of a
  // solver failure.
  if (solver.getStatus() == NOX::StatusTest::Failed && !write_on_solver_failure)
  {
    return;
  }
  
  //Extract contact voltages from the parameter library.
  contactVoltages_ = pPE->get_VoltageParameters();

  RCP<NOX::Abstract::Vector const> x;
  if (loca_solver || loca_householder_solver)
    x = get_loca_x(solver);
  else
    x = get_nox_x(solver);

  NOX::Thyra::Vector const& nThX = Teuchos::dyn_cast<NOX::Thyra::Vector const>(*x);
  RCP<VectorBase<double> const> const& thX = nThX.getThyraRCPVector();

  // The LOCA Householder ME wrapper converts the current responses
  // into the residual equation responses. To get the true current
  // values, we need to use the underlying ME.
  const RCP<const CCME_LOCA> ccme_loca = rcp_dynamic_cast<const CCME_LOCA>(full_me,false);
  if (loca_householder_solver && (nonnull(ccme_loca)))
    full_me = ccme_loca->getInternalPhysicsME();

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

      std::string checkName = pPE->getContinuationParameterName();
      int cvStartIndex = 0;
      cvToWrite = cvNames.size();
      if ( checkName != "NoName")
	{
	  continuationParameterName = checkName;
	  cvStartIndex = 1;
	  cvToWrite -= 1;
	}
      //put contact values into a vector for output.
      contactVoltages_ = pPE->get_VoltageParameters();
      std::vector<double> tempVal, cvVal;
      tempVal.resize(cvNames.size());
      cvVal.resize(cvNames.size());
      int cvCounter=0;
      std::map<std::string,double>::iterator pLit = contactVoltages_.begin();
      while (pLit != contactVoltages_.end())
	{
	  cvVal[cvPermutation[cvCounter]] = pLit->second;
	  ++pLit;
	  ++cvCounter;
	}

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
        vector<size_t> column(Ng + 1 + cvToWrite);
        column[0] = max(continuationParameterName.size(),
          value_width);
	for (size_t i=cvStartIndex ; i<cvNames.size() ; ++i)
	  column[i+1] = max(cvNames[i].size(),value_width);
        for (int i(0); i < Ng; ++i)
          column[i + 1 + cvToWrite] = max(response_names[i].size(), value_width);

        // This just writes the header out, voltage and currents, the first
        // time through.
        if (!out_file_header_written)
        {
          out_file_header_written = true;

          *f_out << setw(column[0]) << continuationParameterName << " ";
	  for (size_t i=cvStartIndex ; i<cvNames.size() ; ++i)
	    *f_out<< setw(column[i+1]) << cvNames[i] << " ";

          for (int i(0); i < Ng - 1; ++i)
            *f_out << setw(column[i + 1 + cvToWrite]) << response_names[i] << " ";
          *f_out << setw(column[Ng]) << response_names[Ng - 1] << endl;
        }

        // The parameter (voltage) value

	if (continuationParameterName.find("Varying Voltage")==0)
	  *f_out << setw(column[0]) << cvVal[0] << " ";
	else
	  *f_out << setw(column[0]) << continuationParameterValue << " ";

	//write the other contact voltages
	for(size_t i=cvStartIndex ; i<cvNames.size() ; ++i)
	  *f_out <<setw(column[i+1]) << cvVal[i] << " ";

        // The responses (currents)
        for (int i=0; i < Ng-1; ++i)
        {
          RCP<VectorBase<double> > g = outArgs.get_g(i);

          // For total correctness this should loop over the dimensions
          // of the parameters but we only calculate a single current
          // value right now. See the screen output below for an example
          // of how we might write out values for more than one
          // dimension here.
          *f_out << setw(column[i + 1 + cvToWrite]) << Thyra::get_ele(*g, 0) << " ";
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

      std::map<std::string,double>::iterator pLit = contactVoltages_.begin();
      while (pLit != contactVoltages_.end())
	{
	  m_out << "    " << pLit->first << " = ";
          m_out << pLit->second << " ";
	  m_out << endl;
	  ++pLit;
	}

      // Check to see if we have any current constraints.
      RCP<const CCME> ccme = rcp_dynamic_cast<const CCME>(full_me);
      if (nonnull(ccme))
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
      }
      else if (nonnull(ccme_loca))
      {
        const auto& nox_grp = solver.getSolutionGroup();
        const auto* loca_mc_grp = dynamic_cast<const LOCA::MultiContinuation::AbstractGroup*>(&nox_grp);
        if (loca_mc_grp == nullptr) {
          const auto* check = dynamic_cast<const LOCA::Extended::MultiAbstractGroup*>(&nox_grp);
          if (check != nullptr) {
            loca_mc_grp = check->getBaseLevelUnderlyingGroup().getRawPtr();
          }
        }
        m_out << "Voltage values:" << endl;
        const CurrentConstraintList& constraints = ccme_loca->constraints();
        for (int i(0); i < constraints.size(); ++i)
          m_out << "  " << constraints[i]->sidesetId()
                << constraints[i]->type() << "Voltage = " << loca_mc_grp->getParam(constraints[i]->parameterName()) << endl;
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

  // If using Householder constraints with continuation, the vector is
  // double nested.
  const auto* test = dynamic_cast<const LOCA::MultiContinuation::ExtendedVector*>(n_x.getRawPtr());
  if (test != nullptr)
    return test->getXVec();

  return n_x;

}

}
