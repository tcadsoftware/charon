
#ifndef CHARON_RYTHMOS_OBSERVER_CLUSTEROBSERVER_HPP
#define CHARON_RYTHMOS_OBSERVER_CLUSTEROBSERVER_HPP

#include <Charon_config.hpp>

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Panzer_GlobalIndexer.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_GlobalData.hpp"

#include "Panzer_STK_Utilities.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PanzerParameterExtractor.hpp"

#include <fstream>


namespace charon {

  /** This class is a rythmos observer that manages cluster operations. It makes use of the response_library.
    */
  class RythmosObserver_ClusterObserver :
    public Rythmos::IntegrationObserverBase<double> {

  public:

   //Add vector of global contact voltages to write to exodus
    std::map<std::string,double> contactVoltages_;
    Teuchos::RCP<panzerParameterExtractor>  pPE;
    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;
 
 
    RythmosObserver_ClusterObserver(const std::vector<std::string> & response_names,
                                    const std::vector<std::string> & parameter_names,
				    Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                                   const Teuchos::RCP<charon::Scaling_Parameters> & scale_params);


    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    cloneIntegrationObserver() const
    { return Teuchos::rcp(new RythmosObserver_ClusterObserver(m_response_names,m_parameter_names,parameterLibrary_,m_scale_params)); }

    void resetIntegrationObserver(const Rythmos::TimeRange<double> &integrationTimeDomain);

    void observeStartTimeIntegration(const Rythmos::StepperBase<double> &stepper);
    void observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                                  const Rythmos::StepControlInfo<double> &stepCtrlInfo,
                                  const int timeStepIter);

    void observeEndTimeIntegration(const Rythmos::StepperBase<double> &stepper);

  private:

    bool executeResponses(const Rythmos::StepperBase<double> &stepper);

  private:

    // For unscaling the times to values in seconds
    Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;

    // Local storage for response names
    std::vector<double> m_current_responses;
    std::vector<std::string> m_response_names;

    // local storage for voltage parameters
    std::vector<double> m_parameter_values;
    std::vector<std::string> m_parameter_names;

    //Keep track of where we are in break points
    size_t breakPointIndex;


  };

}

#endif  // CHARON_RYTHMOS_OBSERVER_CLUSTEROBSERVER_HPP
