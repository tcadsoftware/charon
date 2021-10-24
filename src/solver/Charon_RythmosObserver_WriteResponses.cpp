
#include <algorithm>

// Charon
#include "Charon_CurrentConstraintModelEvaluator.hpp"
#include "Charon_RythmosObserver_WriteResponses.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_SpmdVectorBase.hpp"

namespace charon {

RythmosObserver_WriteResponses::
RythmosObserver_WriteResponses(const std::vector<std::string> & response_names,
                               const std::vector<std::string> & parameter_names,
                               const Teuchos::RCP<charon::Scaling_Parameters> & scale_params,
                               bool write_to_screen,
                               const std::string & filename) :
  m_out(Teuchos::rcpFromRef(std::cout)),
  m_file_output(false),
  m_write_to_screen(write_to_screen),
  m_out_field_size(32), // Initialized to the inimum field size for
                        // output file
  m_scale_params(scale_params),
  m_response_names(response_names),
  m_parameter_names(parameter_names),
  m_filename(filename)
{
  TEUCHOS_ASSERT(m_scale_params!=Teuchos::null);

  // setup fancy ostream for output
  m_out.setOutputToRootOnly(0);
  m_out.setTabIndentStr(" |   ");

  if(filename != "") {
    m_file_output = true;
  }
}

void
RythmosObserver_WriteResponses::
observeStartTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */)
{
  if(m_response_names.size()==0) {
    m_file_output = false;
    m_write_to_screen = false;
  }

  // setup file output for responses
  if(m_file_output) {

    Teuchos::RCP<std::ofstream> fileOutput = Teuchos::rcp(new std::ofstream());
    fileOutput->open(m_filename.c_str(), std::ofstream::out | std::ofstream::app);

    Teuchos::RCP<Teuchos::FancyOStream> fout = Teuchos::rcp(new Teuchos::FancyOStream(fileOutput));

    fout->setOutputToRootOnly(0);
    fout->setTabIndentStr(" |   ");

    // Find the maximum string length and size the output fields
    // accordingly
    for (std::size_t i=0; i < m_response_names.size(); ++i)
      m_out_field_size = std::max(m_out_field_size, m_response_names[i].size());

    for (std::size_t i=0; i < m_parameter_names.size(); ++i)
      m_out_field_size = std::max(m_out_field_size, m_parameter_names[i].size());

    // this writes out the header for the table
    fout->width(m_out_field_size);
    *fout << "Time";

    // output all responses of functional type (with their name)

    for(std::size_t i=0;i<m_response_names.size();i++) {
      fout->width(m_out_field_size);

      char const* spcchar = "\t";
      if ((i+1) == m_response_names.size())
        spcchar = ""; // No tab at the end of the line

      *fout << m_response_names[i] << spcchar;
    }

    for(std::size_t i=0;i<m_parameter_names.size();i++) {
      fout->width(m_out_field_size);

      char const* spcchar = "\t";
      if ((i+1) == m_response_names.size())
        spcchar = ""; // No tab at the end of the line

      *fout << m_parameter_names[i] << spcchar;
    }

    *fout << std::endl;

    fileOutput->close();
  }
}

void
RythmosObserver_WriteResponses::
writeToScreen(double t)
{
  m_out << "Responses at t=" << t << "\n";
  m_out.pushTab();

  // output all responses of functional type (with their name)
  for(std::size_t i=0;i<m_response_names.size();i++)
    m_out << m_response_names[i] << " = " << m_current_responses[i] << std::endl;

  for(std::size_t i=0;i<m_parameter_names.size();i++)
    m_out << m_parameter_names[i] << " = " << m_parameter_values[i] << std::endl;

  m_out.popTab();
  m_out << std::endl;
}

void
RythmosObserver_WriteResponses::
writeToFile(double t)
{

  // Open and close the file for each write to avoid issues with
  // buffering so that the user can see values as the simulation
  // progresses.
  Teuchos::RCP<std::ofstream> fileOutput = Teuchos::rcp(new std::ofstream());
  fileOutput->open(m_filename.c_str(), std::ofstream::out | std::ofstream::app);

  Teuchos::RCP<Teuchos::FancyOStream> fout = Teuchos::rcp(new Teuchos::FancyOStream(fileOutput));

  fout->setOutputToRootOnly(0);
  fout->setTabIndentStr(" |   ");

  fout->width(m_out_field_size);
  *fout <<std::scientific<<std::setprecision(10);
  *fout << t;

  // output all responses of functional type (with their name)
  for(std::size_t i=0;i<m_response_names.size();i++) {
    fout->width(m_out_field_size);
    *fout << m_current_responses[i];
  }

  // output all responses of functional type (with their name)
  for(std::size_t i=0;i<m_parameter_names.size();i++) {
    fout->width(m_out_field_size);
    *fout << m_parameter_values[i];
  }

  *fout << std::endl;

  fileOutput->close();
}

void
RythmosObserver_WriteResponses::
observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                         const Rythmos::StepControlInfo<double> & /* stepCtrlInfo */,
                         const int /* timeStepIter */)
{
  // compute responses, we will write them in a minute
  this->executeResponses(stepper);

  const Rythmos::StepStatus<double> status = stepper.getStepStatus();

  // time scaling
  double t = status.time*m_scale_params->scale_params.t0;

  //if(m_write_to_screen)
    writeToScreen(t);

  if(m_file_output)
    writeToFile(t);
}

bool
RythmosObserver_WriteResponses::
executeResponses(const Rythmos::StepperBase<double> &stepper)
{
  using   std::size_t;
  using   Teuchos::ArrayRCP;
  using   Teuchos::ptrFromRef;
  using   Teuchos::RCP;
  using   Teuchos::rcp_dynamic_cast;
  using   Thyra::createMember;
  using   Thyra::ModelEvaluator;
  using   Thyra::ProductVectorBase;
  using   Thyra::SpmdVectorBase;
  using   Thyra::VectorBase;
  typedef charon::CurrentConstraintModelEvaluator<double> CCME;
  typedef Thyra::ModelEvaluatorBase::OutArgs<double>      OutArgs;
  typedef Thyra::ModelEvaluatorBase::InArgs<double>       InArgs;

  RCP<const VectorBase<double>> thX, thXDot;
  double time(0.0);

  // Clear out any previous values.
  m_current_responses.clear();
  m_parameter_values.clear();

  time = stepper.getStepStatus().time;

  thX = stepper.getStepStatus().solution;
  if (thX == Teuchos::null)
    thX = Rythmos::get_x(stepper,time);
  TEUCHOS_ASSERT(thX != Teuchos::null);

  thXDot = stepper.getStepStatus().solutionDot;
  if (thXDot == Teuchos::null)
    thXDot = Rythmos::get_xdot(stepper,time);                                    // ???:  This is why Rythmos makes me angry.  I
  TEUCHOS_ASSERT(thXDot != Teuchos::null);                                       //       can't use
                                                                                 //       stepper.getStepStatus().solution.  Was
  // Get the model and create an OutArgs object.                                 //       this the intended approach?  Who would
  RCP<const ModelEvaluator<double>> model   = stepper.getModel();                //       ever  want to use the time derivative
  OutArgs                           outArgs = model->createOutArgs();            //       after solving?  This last question is
                                                                                 //       sarcasm.
  // Compute the responses only if there are any to compute.
  RCP<const CCME> ccme;
  if (outArgs.Ng() > 0)
  {
    // Build an InArgs object to evaluate the responses.
    InArgs basePoint = stepper.getInitialCondition();
    InArgs inArgs    = model->createInArgs();
    inArgs.setArgs(basePoint);
    inArgs.set_x(thX);
    inArgs.set_x_dot(thXDot);
    inArgs.set_t(time);

    // Populate the OutArgs with the response vectors.
    for (int i(0); i < outArgs.Ng(); ++i)
      outArgs.set_g(i, createMember(*model->get_g_space(i)));

    // Evaluate the model.
    model->evalModel(inArgs, outArgs);

    // Distribute the response vectors back to the single response vector for
    // output purposes.
    m_current_responses.resize(outArgs.Ng());
    for (int i(0); i < outArgs.Ng(); ++i)
    {
      ArrayRCP<const double> gData;
      rcp_dynamic_cast<SpmdVectorBase<double>>(outArgs.get_g(i), true)->
        getLocalData(ptrFromRef(gData));
      TEUCHOS_ASSERT(gData.size() == 1);
      m_current_responses[i] = gData[0];
    } // end loop over outArgs

    // Extract the voltage values as well if we have any current constraints.
    ccme = rcp_dynamic_cast<const CCME>(model);
    if (not ccme.is_null())
    {
      RCP<const VectorBase<double>> volt =
        rcp_dynamic_cast<const ProductVectorBase<double>>(thX, true)->
        getVectorBlock(1);
      ArrayRCP<const double> vData;
      rcp_dynamic_cast<const SpmdVectorBase<double>>(volt, true)->
        getLocalData(ptrFromRef(vData));
      TEUCHOS_ASSERT(vData.size() ==
        static_cast<unsigned int>(m_parameter_names.size()));
      m_parameter_values.resize(m_parameter_names.size());
      for (size_t i(0); i < m_parameter_names.size(); ++i)
        m_parameter_values[i] = vData[i];
    } // end if we have any current constraints
  } // end if (outArgs.Ng() > 0)

  // A simple precaution/sanity check.
  TEUCHOS_ASSERT(m_response_names.size() == m_current_responses.size());

  // Return true if we were dealing with a current constraint case.
  if (not ccme.is_null())
    return true;
  else
    return false;
} // end of executeResponses()

}
