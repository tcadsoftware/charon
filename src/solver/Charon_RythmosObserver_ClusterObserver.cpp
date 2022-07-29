
// C++
#include <iostream>

// Charon
#include "Charon_RythmosObserver_ClusterObserver.hpp"

// Rythmos
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_SpmdVectorBase.hpp"

namespace charon {

RythmosObserver_ClusterObserver::
RythmosObserver_ClusterObserver(const std::vector<std::string> & response_names,
                                const std::vector<std::string> & parameter_names,
				Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                               const Teuchos::RCP<charon::Scaling_Parameters> & scale_params) :
  parameterLibrary_(parameterLibrary),
  m_scale_params(scale_params),
  m_response_names(response_names),
  m_parameter_names(parameter_names)
{
  TEUCHOS_ASSERT(m_scale_params!=Teuchos::null);
}





void
RythmosObserver_ClusterObserver::
observeStartTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */)
{
}


void
RythmosObserver_ClusterObserver::
resetIntegrationObserver(const Rythmos::TimeRange<double> & /* integrationTimeDomain */)
{

}

void
RythmosObserver_ClusterObserver::
observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                         const Rythmos::StepControlInfo<double> & /* stepCtrlInfo */,
                         const int /* timeStepIter */)
{


  const Rythmos::StepStatus<double> status = stepper.getStepStatus();



}


void
RythmosObserver_ClusterObserver::
observeEndTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */)
{

}



}
