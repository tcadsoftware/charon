
#ifndef CHARON_RYTHMOS_OBSERVER_WRITETOEXODUS_HPP
#define CHARON_RYTHMOS_OBSERVER_WRITETOEXODUS_HPP

// Charon
#include "Charon_Scaling_Parameters.hpp"

// Panzer
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_Utilities.hpp"

// Rythmos
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_TimeRange.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_SpmdVectorBase.hpp"

namespace charon {

  class RythmosObserver_WriteToExodus :
    public Rythmos::IntegrationObserverBase<double> {

  public:

    RythmosObserver_WriteToExodus(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                                   const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
                                   const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof,
                                   const int write_after_this_many_steps,
                                   const bool write_initial_condition,
                                   const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & response_library,
                                   const Teuchos::RCP<charon::Scaling_Parameters> & scale_params) :
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_current_step(0),
      m_last_step_written(0),
      m_write_every_n_steps(write_after_this_many_steps),
      m_write_initial_condition(write_initial_condition),
      m_response_library(response_library),
      m_out(Teuchos::rcpFromRef(std::cout)),
      m_scale_params(scale_params)
    {
      TEUCHOS_ASSERT(m_lof!=Teuchos::null);
      TEUCHOS_ASSERT(m_response_library!=Teuchos::null);

      // get all element blocks and add them to the list
      std::vector<std::string> eBlocks;
      mesh->getElementBlockNames(eBlocks);

      panzer_stk::RespFactorySolnWriter_Builder builder;
      builder.mesh = mesh;
      m_response_library->addResponse("Main Field Output",eBlocks,builder);

      // setup fancy ostream for output
      m_out.setOutputToRootOnly(0);
      m_out.setTabIndentStr(" |   ");
    }

    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    cloneIntegrationObserver() const
    {
      return Teuchos::rcp(new RythmosObserver_WriteToExodus(m_mesh, m_dof_manager, m_lof, m_write_every_n_steps, m_write_initial_condition,m_response_library, m_scale_params));
    }

    void resetIntegrationObserver(const Rythmos::TimeRange<double> & /* integrationTimeDomain */)
    {
      m_current_step = 0;
      m_last_step_written = 0;
    }

    void observeStartTimeIntegration(const Rythmos::StepperBase<double> &stepper)
    {
      const Rythmos::StepStatus<double> status = stepper.getStepStatus();

      // time scaling
      double t0 = m_scale_params->scale_params.t0;

      m_out << std::endl;
      m_out << "Begin time integration" << std::endl;
      m_out.pushTab();
      m_out << "Time:      " << status.time *t0 << " in seconds" << std::endl;
      m_out << "Step Size: " << status.stepSize *t0 << " in seconds" << std::endl;

      if (m_write_initial_condition) {
        m_out << "Writing solution to exodus" << std::endl;
        m_out << "   ... ";

        this->writeSolutionToExodus(stepper, true);

        m_out << "complete" << std::endl;
      }

      m_out.popTab();
      m_out << std::endl;
    }

    void observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                                  const Rythmos::StepControlInfo<double> & /* stepCtrlInfo */,
                                  const int timeStepIter)
    {
      m_current_step = timeStepIter;

      const Rythmos::StepStatus<double> status = stepper.getStepStatus();

      // time scaling
      double t0 = m_scale_params->scale_params.t0;

      m_out << std::endl;
      m_out << "Completed time step" << std::endl;
      m_out.pushTab();
      m_out << "Time:      " << status.time *t0 << " in seconds" << std::endl;
      m_out << "Step:      " << timeStepIter << std::endl;
      m_out << "Step Size: " << status.stepSize *t0 << " in seconds" << std::endl;

      if (((timeStepIter +1)  % m_write_every_n_steps) == 0) {
        m_out << "Writing solution to exodus" << std::endl;
        m_out << "   ... ";

        this->writeSolutionToExodus(stepper);
        m_last_step_written = timeStepIter;

        m_out << "complete" << std::endl;
      }

      m_out.popTab();
      m_out << std::endl;
    }

    void observeEndTimeIntegration(const Rythmos::StepperBase<double> &stepper)
    {
      if (m_last_step_written != m_current_step)
      this->writeSolutionToExodus(stepper);
    }

  private:

    void writeSolutionToExodus(const Rythmos::StepperBase<double> &stepper, bool use_initial_condition = false)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;

      TEUCHOS_ASSERT(m_lof!=Teuchos::null);

      // time scaling
      double t0 = m_scale_params->scale_params.t0;

      Teuchos::RCP<const Thyra::VectorBase<double> > th_x;
      double time = 0.0;
      if (use_initial_condition) {
        th_x = stepper.getInitialCondition().get_x();
        time = stepper.getInitialCondition().get_t();
      }
      else {
        th_x = stepper.getStepStatus().solution;
        time = stepper.getStepStatus().time *t0;
      }

      // If we're using a blocked system (for instance, for the constant
      // current or resistor contact cases), we need to extract only the PDE
      // degrees of freedom, which are contained in the first vector block.
      const RCP<const Thyra::ProductVectorBase<double> > blo_x =
        Thyra::castOrCreateProductVectorBase(th_x);
      if (blo_x->productSpace()->numBlocks() > 1)
        th_x = blo_x->getVectorBlock(0);

      // initialize the assembly container
      panzer::AssemblyEngineInArgs ae_inargs;
      ae_inargs.container_ = m_lof->buildLinearObjContainer();
      ae_inargs.ghostedContainer_ = m_lof->buildGhostedLinearObjContainer();
      ae_inargs.alpha = 0.0;
      ae_inargs.beta = 1.0;
      ae_inargs.evaluate_transient_terms = false;

      // initialize the ghosted container
      m_lof->initializeGhostedContainer(panzer::LinearObjContainer::X,*ae_inargs.ghostedContainer_);

      Teuchos::MpiComm<int> comm = m_lof->getComm();

      {
         // initialize the x vector
         const Teuchos::RCP<panzer::ThyraObjContainer<double> > thyraContainer
            = Teuchos::rcp_dynamic_cast<panzer::ThyraObjContainer<double> >(ae_inargs.container_,true);
         thyraContainer->set_x_th(Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(th_x));
      }

      // force an output to STK mesh
      m_response_library->addResponsesToInArgs<panzer::Traits::Residual>(ae_inargs);
      m_response_library->evaluate<panzer::Traits::Residual>(ae_inargs);

      // write to disk
      m_mesh->writeToExodus(time);
    }

  private:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<const panzer::GlobalIndexer> m_dof_manager;
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > m_lof;

    /** Stores the current step */
    int m_current_step;

    /** Last time step that was written.  Used to determine if at the
        end of the simulation that the final time step is written.
        Used in conjunction with the output every N steps feature.
    */
    int m_last_step_written;

    /** Instead of writing every time step to exodus, this parameter
        lets a user write only every Nth time step.
    */
    int m_write_every_n_steps;

    bool m_write_initial_condition;

    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > m_response_library;

    // For outputting time step information
    Teuchos::FancyOStream m_out;

    // For unscaling the times to values in seconds
    Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;

  };

}

#endif
