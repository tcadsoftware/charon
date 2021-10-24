
#ifndef CHARON_RYTHMOS_OBSERVER_WRITERESPONSES_HPP
#define CHARON_RYTHMOS_OBSERVER_WRITERESPONSES_HPP

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Panzer_GlobalIndexer.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"

#include "Panzer_STK_Utilities.hpp"
#include "Charon_Scaling_Parameters.hpp"

#include <fstream>

namespace charon {

  /** This class is a rythmos observer that outputs all the respones to either the screen
    * and/or a file. This is done every time step. It makes use of the response_library.
    */
  class RythmosObserver_WriteResponses :
    public Rythmos::IntegrationObserverBase<double> {

  public:

    RythmosObserver_WriteResponses(const std::vector<std::string> & response_names,
                                   const std::vector<std::string> & parameter_names,
                                   const Teuchos::RCP<charon::Scaling_Parameters> & scale_params,
                                   bool write_to_screen=true,
                                   const std::string & filename="");

    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    cloneIntegrationObserver() const
    { return Teuchos::rcp(new RythmosObserver_WriteResponses(m_response_names,m_parameter_names,m_scale_params)); }

    void resetIntegrationObserver(const Rythmos::TimeRange<double> & /* integrationTimeDomain */) {}

    void observeStartTimeIntegration(const Rythmos::StepperBase<double> &stepper);
    void writeToScreen(double t);
    void writeToFile(double t);
    void observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                                  const Rythmos::StepControlInfo<double> &stepCtrlInfo,
                                  const int timeStepIter);

    void observeEndTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */) {}

  private:

    bool executeResponses(const Rythmos::StepperBase<double> &stepper);

  private:

    // For outputting time step information
    Teuchos::FancyOStream m_out;

    bool m_file_output;

    bool m_write_to_screen;

    // Maximum size of the fields. This is the longest string that is in
    // either m_response_names or m_parameter_names.
    size_t m_out_field_size;

    // For unscaling the times to values in seconds
    Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;

    // Local storage for response names
    std::vector<double> m_current_responses;
    std::vector<std::string> m_response_names;

    // local storage for voltage parameters
    std::vector<double> m_parameter_values;
    std::vector<std::string> m_parameter_names;

    std::string m_filename;
  };

}

#endif
