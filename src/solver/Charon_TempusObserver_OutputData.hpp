#ifndef CHARON_TEMPUS_OBSERVER_OUTPUTDATA_HPP
#define CHARON_TEMPUS_OBSERVER_OUTPUTDATA_HPP

// Charon
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PanzerParameterExtractor.hpp"

// Panzer
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_STK_Utilities.hpp"
#include "Panzer_GlobalData.hpp"

// Tempus
#include "Tempus_IntegratorObserverBasic.hpp"
#include "Tempus_StepperBackwardEuler_decl.hpp"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_SpmdVectorBase.hpp"

namespace charon {

  class TempusObserver_OutputData :
    public Tempus::IntegratorObserverBasic<double> {

  public:

    //Add vector of global contact voltages to write to exodus
    std::map<std::string,double> contactVoltages_;
    std::map<std::string,double> norms_;
    Teuchos::RCP<panzerParameterExtractor>  pPE;
    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;

    TempusObserver_OutputData(Teuchos::RCP<panzer_stk::STK_Interface> const& mesh,
                              Teuchos::RCP<const panzer::GlobalIndexer> const& dof_manager,
                              Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > const& lof,
                              Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > const& response_library,
                              Teuchos::RCP<Thyra::ModelEvaluator<double>> const& full_me,
                              std::vector<std::string> const& response_names,
                              std::vector<std::string> const& param_names,
                              Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                              Teuchos::RCP<charon::Scaling_Parameters> const& scale_params,
                              bool output_exodus,
                              bool output_responses,
                              bool output_txt_response_file,
                              bool output_exo_responses,
                              int write_after_this_many_steps,
                              std::string const& txt_response_filename,
                              bool write_initial_condition) :
      parameterLibrary_(parameterLibrary),
      m_mesh(mesh),
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_response_library(response_library),
      m_full_me(full_me),
      m_response_names(response_names),
      m_param_names(param_names),
      m_scale_params(scale_params),
      m_output_exodus(output_exodus),
      m_output_responses(output_responses),
      m_output_txt_response_file(output_txt_response_file),
      m_output_exo_responses(output_exo_responses),
      m_write_every_n_steps(write_after_this_many_steps),
      m_txt_response_filename(txt_response_filename),
      m_write_initial_condition(write_initial_condition),
      m_out_field_size(8),
      m_current_step(0),
      m_last_step_written(0)
    {
      TEUCHOS_ASSERT(m_lof!=Teuchos::null);
      TEUCHOS_ASSERT(m_response_library!=Teuchos::null);

      // get all element blocks and add them to the list
      std::vector<std::string> eBlocks;
      mesh->getElementBlockNames(eBlocks);

      panzer_stk::RespFactorySolnWriter_Builder builder;
      builder.mesh = mesh;
      m_response_library->addResponse("Main Field Output",eBlocks,builder);

      // time scaling
      m_t0 = m_scale_params->scale_params.t0;

      //create the panzer parameter extractor object
      pPE = Teuchos::rcp(new panzerParameterExtractor(parameterLibrary_));

    }

    void observeStartIntegrator(Tempus::Integrator<double> const& l_integrator)
    {

      Teuchos::RCP<Teuchos::FancyOStream> const out = l_integrator.getOStream();
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(out,0,"ScreenOutput");
      *out << "  Stepper = " << l_integrator.getStepper()->description() << "\n"
           << "  Simulation Time Range  [" << l_integrator.getTimeStepControl()->getInitTime() * m_t0
           << ", " << l_integrator.getTimeStepControl()->getFinalTime() * m_t0 << "]\n"
           << "  Simulation Index Range [" << l_integrator.getTimeStepControl()->getInitIndex()
           << ", " << l_integrator.getTimeStepControl()->getFinalIndex() << "]\n"
           << std::endl;

      if (m_output_responses)
        calculate_responses(l_integrator, true);

      if (m_write_initial_condition)
      {
        *out << "Writing solution to exodus" << std::endl;
        *out << "   ... ";

        this->writeSolutionToExodus(l_integrator, m_t0*l_integrator.getTime(), true);

        *out << "complete" << std::endl;
      }

      *out << std::endl;

      // setup text file for response output, if requested
      if (m_output_txt_response_file)
      {

        Teuchos::RCP<std::ofstream> fileOutput = Teuchos::rcp(new std::ofstream());
        fileOutput->open(m_txt_response_filename.c_str(), std::ofstream::out | std::ofstream::app);

        Teuchos::RCP<Teuchos::FancyOStream> fout = Teuchos::rcp(new Teuchos::FancyOStream(fileOutput));

        fout->setOutputToRootOnly(0);
        fout->setTabIndentStr(" |   ");

        // Find the maximum string length and size the output fields
        // accordingly

        for (std::size_t i=0; i < m_response_names.size(); ++i)
          m_out_field_size = std::max(m_out_field_size, m_response_names[i].size());

        for (std::size_t i=0; i < m_param_names.size(); ++i)
          m_out_field_size = std::max(m_out_field_size, m_param_names[i].size());

        // this writes out the header for the table
        fout->width(m_out_field_size);
        *fout << "Time";

        //output contact voltage names
        contactVoltages_ = pPE->get_VoltageParameters();
        std::map<std::string,double>::iterator cVit = contactVoltages_.begin();
        while (cVit != contactVoltages_.end())
        {
          fout->width(m_out_field_size);
          *fout << cVit->first;
          ++cVit;
        }

        // output all responses of functional type (with their name)
        for(std::size_t i=0;i<m_response_names.size();i++)
        {
          fout->width(m_out_field_size);

          char const* spcchar = "\t";
          if ((i+1) == m_response_names.size())
            spcchar = ""; // No tab at the end of the line

          *fout << m_response_names[i] << spcchar;
        }

        for(std::size_t i=0;i<m_param_names.size();i++)
        {
          fout->width(m_out_field_size);

          char const* spcchar = "\t";
          if ((i+1) == m_response_names.size())
            spcchar = ""; // No tab at the end of the line

          *fout << m_param_names[i] << spcchar;
        }

        *fout << std::endl;

        fileOutput->close();
      }

      // Output initial responses
      if(m_output_responses)
        write_to_screen(out, m_t0*l_integrator.getTime());

      if (m_output_txt_response_file)
        write_to_file(0.0);

    }

    void observeEndTimeStep(Tempus::Integrator<double> const& l_integrator)
    {
      Teuchos::RCP<Teuchos::FancyOStream> const out = l_integrator.getOStream();

      // Determine if the time step failed or succeeded via a test of
      // the working state, if it's null proceed as if the step were
      // successful.
      if (Teuchos::is_null(l_integrator.getSolutionHistory()->getWorkingState(false)))
      {
        auto cs = l_integrator.getSolutionHistory()->getCurrentState();
        
        m_current_step = cs->getIndex();

        *out << std::endl;
        *out << "Completed Tempus time step" << std::endl;
        *out << "Time:      " << cs->getTime() * m_t0 << " in seconds" << std::endl;
        *out << "Step:      " << m_current_step << std::endl;
        *out << "Step Size: " << cs->getTimeStep() * m_t0 << " in seconds" << std::endl;

        if (m_output_responses)
          calculate_responses(l_integrator, true);
        
        if (((m_current_step+1) % m_write_every_n_steps) == 0) {
          *out << "Writing solution to exodus" << std::endl;
          *out << "   ... ";
          
          this->writeSolutionToExodus(l_integrator, m_t0*cs->getTime());
          m_last_step_written = m_current_step;
          
          *out << "complete" << std::endl;
        }
        
        *out << std::endl;
        
        if(m_output_responses)
          write_to_screen(out, m_t0*cs->getTime());
        
        if (m_output_txt_response_file)
          write_to_file(m_t0*cs->getTime());
      }
      else // step failed
      {
        *out << "Completed Tempus time step" << std::endl;
        *out << "Warning, a time step failure occurred" << std::endl;
      }
    }

    void observeEndIntegrator(Tempus::Integrator<double> const& l_integrator)
    {
      auto cs = l_integrator.getSolutionHistory()->getCurrentState();

      if (m_output_responses)
        calculate_responses(l_integrator, true);

      if (m_last_step_written != m_current_step)
        this->writeSolutionToExodus(l_integrator, m_t0*cs->getTime());
    }

  private:

    void write_to_screen(Teuchos::RCP<Teuchos::FancyOStream> const out, double const s_time)
    {
      *out << "Responses at t=" << s_time << "\n";

      // output all responses of functional type (with their name)
      for(std::size_t i=0;i<m_response_names.size();i++)
        *out << m_response_names[i] << " = " << m_response_values[i] << std::endl;

      for(std::size_t i=0;i<m_param_names.size();i++)
        *out << m_param_names[i] << " = " << m_param_values[i] << std::endl;

      *out << std::endl;
    }

    void write_to_file(double const s_time)
    {
      // Open and close the file for each write to avoid issues with
      // buffering so that the user can see values as the simulation
      // progresses.
      Teuchos::RCP<std::ofstream> fileOutput = Teuchos::rcp(new std::ofstream());
      fileOutput->open(m_txt_response_filename.c_str(), std::ofstream::out | std::ofstream::app);

      Teuchos::RCP<Teuchos::FancyOStream> fout = Teuchos::rcp(new Teuchos::FancyOStream(fileOutput));

      fout->setOutputToRootOnly(0);
      fout->setTabIndentStr(" |   ");

      fout->width(m_out_field_size);
      *fout <<std::scientific<<std::setprecision(10);
      *fout << s_time;

      //Add contact voltages to csv file
      //Get contact voltage names and values
      contactVoltages_ = pPE->get_VoltageParameters();
      std::map<std::string,double>::iterator cVit = contactVoltages_.begin();
      while (cVit != contactVoltages_.end())
      {
        fout->width(m_out_field_size);
        *fout << cVit->second;
        ++cVit;
      }

      // output all responses of functional type (with their name)
      for(std::size_t i=0;i<m_response_names.size();i++)
      {
        fout->width(m_out_field_size);
        *fout << m_response_values[i];
      }

      // output all responses of functional type (with their name)
      for(std::size_t i=0;i<m_param_names.size();i++)
      {
        fout->width(m_out_field_size);
        *fout << m_param_values[i];
      }

      *fout << std::endl;

      fileOutput->close();
    }

    // calculate the responses for all output techniques
    void calculate_responses(Tempus::Integrator<double> const& l_integrator,
                             bool use_initial_condition=false)
    {


      m_response_values.clear();
      m_param_values.clear();

      auto cs = l_integrator.getSolutionHistory()->getCurrentState();
      double s_time = cs->getTime();

      Teuchos::RCP<Thyra::VectorBase<double> const> thX = cs->getX();
      TEUCHOS_ASSERT(thX != Teuchos::null);

      Teuchos::RCP<Thyra::VectorBase<double> const> thXDot = cs->getXDot();
      TEUCHOS_ASSERT(thXDot != Teuchos::null);

      // Get the model and create an OutArgs object.
      Teuchos::RCP<Thyra::ModelEvaluator<double> const> model = l_integrator.getStepper()->getModel();
      Thyra::ModelEvaluatorBase::OutArgs<double> out_args = model->createOutArgs();

      // Compute the responses only if there are any to compute.
      Teuchos::RCP<charon::CurrentConstraintModelEvaluator<double> const> ccme;
      if (out_args.Ng() > 0)
      {
        // Build an InArgs object to evaluate the responses.
        // Thyra::ModelEvaluatorBase::InArgs<double> base_point = stepper.getInitialCondition();
        Thyra::ModelEvaluatorBase::InArgs<double> base_point = model->createInArgs();
        base_point.set_x((*l_integrator.getSolutionHistory())[0]->getX());
        base_point.set_x_dot((*l_integrator.getSolutionHistory())[0]->getXDot());
        base_point.set_t((*l_integrator.getSolutionHistory())[0]->getTime());

        Thyra::ModelEvaluatorBase::InArgs<double> in_args = model->createInArgs();
        in_args.setArgs(base_point);
        in_args.set_x(thX);
        in_args.set_x_dot(thXDot);
        in_args.set_t(s_time);

        // Populate the OutArgs with the response vectors.
        for (int i=0; i < out_args.Ng(); ++i)
          out_args.set_g(i, createMember(*model->get_g_space(i)));

        // Evaluate the model.
        model->evalModel(in_args, out_args);

        // Distribute the response vectors back to the single response vector for
        // output purposes.
        m_response_values.resize(out_args.Ng());
        for (int i=0; i < out_args.Ng(); ++i)
        {
          Teuchos::ArrayRCP<const double> gData;
          Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double>>(out_args.get_g(i), true)->
            getLocalData(ptrFromRef(gData));
          TEUCHOS_ASSERT(gData.size() == 1);
          m_response_values[i] = gData[0];
        } // end loop over outArgs

        // Extract the voltage values as well if we have any current constraints.
        ccme = Teuchos::rcp_dynamic_cast<charon::CurrentConstraintModelEvaluator<double> const>(model);
        if (not ccme.is_null())
        {
          Teuchos::RCP<const Thyra::VectorBase<double>> volt =
            Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double>>(thX, true)->
            getVectorBlock(1);
          Teuchos::ArrayRCP<const double> vData;
          Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double>>(volt, true)->
            getLocalData(ptrFromRef(vData));
          TEUCHOS_ASSERT(vData.size() ==
                         static_cast<unsigned int>(m_param_names.size()));
          m_param_values.resize(m_param_names.size());
          for (size_t i=0; i < m_param_names.size(); ++i)
            m_param_values[i] = vData[i];

        } // end if we have any current constraints

      } // end if (out_args.Ng() > 0)

    }

    void writeSolutionToExodus(Tempus::Integrator<double> const& l_integrator,
                               double s_time,
                               bool use_initial_condition = false)
    {

      Teuchos::RCP<const Thyra::VectorBase<double> > th_x;

      if (use_initial_condition)
        th_x = (*l_integrator.getSolutionHistory())[0]->getX();
      else
        th_x = l_integrator.getSolutionHistory()->getCurrentState()->getX();

      // If we're using a blocked system (for instance, for the constant
      // current or resistor contact cases), we need to extract only the PDE
      // degrees of freedom, which are contained in the first vector block.
      const Teuchos::RCP<const Thyra::ProductVectorBase<double> > blo_x =
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

      // add globals to exodus output
      if (m_output_responses)
      {
        for (size_t i=0; i < m_response_names.size(); ++i)
        {
          m_mesh->addGlobalToExodus(m_response_names[i], m_response_values[i]);
        }
      }

      //Add contact voltages to exodus
      //Get contact voltage names and values
      contactVoltages_ = pPE->get_VoltageParameters();
      std::map<std::string,double>::iterator cVit = contactVoltages_.begin();
      while (cVit != contactVoltages_.end())
      {
        m_mesh->addGlobalToExodus(cVit->first, cVit->second);
        ++cVit;
      }

      //Add norms to exodus
      //Get norm names and values
      norms_ = pPE->get_NormParameters();
      std::map<std::string,double>::iterator norm_it = norms_.begin();
      while (norm_it != norms_.end()){
        m_mesh->addGlobalToExodus(norm_it->first, norm_it->second);
        ++norm_it;
      }

      // write to disk
      m_mesh->writeToExodus(s_time);
    }

  private:

    Teuchos::RCP<panzer_stk::STK_Interface> m_mesh;
    Teuchos::RCP<const panzer::GlobalIndexer> m_dof_manager;
    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits>> m_lof;
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> m_response_library;
    Teuchos::RCP<Thyra::ModelEvaluator<double>> m_full_me;

    std::vector<std::string> m_response_names;
    std::vector<double> m_response_values;

    std::vector<std::string> m_param_names;
    std::vector<double> m_param_values;

    Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;

    bool m_output_exodus;
    bool m_output_responses;
    bool m_output_txt_response_file;
    bool m_output_exo_responses;

    /** Instead of writing every time step to exodus, this parameter
        lets a user write only every Nth time step.
    */
    int m_write_every_n_steps;

    std::string m_txt_response_filename;

    bool m_write_initial_condition;

    size_t m_out_field_size;

    /** Stores the current step */
    int m_current_step;

    /** Last time step that was written.  Used to determine if at the
        end of the simulation that the final time step is written. Used
        in conjunction with the output every N steps feature.
    */
    int m_last_step_written;

    bool m_do_precalc;

    double m_t0;
  };

}

#endif
