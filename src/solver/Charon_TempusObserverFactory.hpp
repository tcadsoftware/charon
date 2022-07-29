#ifndef CHARON_TEMPUS_OBSERVER_FACTORY
#define CHARON_TEMPUS_OBSERVER_FACTORY

#include "Charon_config.hpp"
#include "Charon_Scaling_Parameters.hpp"

// concrete tempus observers
#include "Charon_TempusObserver_OutputData.hpp"

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "Tempus_IntegratorObserverComposite_decl.hpp"

#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_ParameterLibrary.hpp"

#include "Panzer_STK_TempusObserverFactory.hpp"

namespace charon {

  class TempusObserverFactory : public panzer_stk::TempusObserverFactory,
                                public Teuchos::ParameterListAcceptorDefaultBase  {

  private:

    mutable Teuchos::RCP<Teuchos::ParameterList> valid_params_;
    bool use_nox_observer_;

    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    //! Store Response names for pretty printing of currents, etc...
    std::vector<std::string> responseNames_;

    Teuchos::RCP<Thyra::ModelEvaluator<double> > full_model_evaluator_;

    //! Store parameter names for pretty printing of voltages, etc...
    std::vector<std::string> parameterNames_;

    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;

    Teuchos::RCP<charon::Scaling_Parameters> scaleParams_;



  public:
    TempusObserverFactory() {}

    TempusObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                           const std::vector<std::string> & responseNames, 
			   Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                           const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams)
       : stkIOResponseLibrary_(stkIOResponseLibrary), responseNames_(responseNames),	
	 parameterLibrary_(parameterLibrary), scaleParams_(scaleParams) {}

    TempusObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
			   Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                           const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams)
       : stkIOResponseLibrary_(stkIOResponseLibrary), parameterLibrary_(parameterLibrary), scaleParams_(scaleParams) {}


    void setResponseNames(const std::vector<std::string> & responseNames)
    { responseNames_ = responseNames; }

    /** Set the name of the parameters to be output, this is only enabled if the constant
      * current calculation is on.
      */
    void setParameterNames(const std::vector<std::string> & parameterNames)
    { parameterNames_ = parameterNames; }

    //! Use NOX observer in solver
    bool useNOXObserver() const
    { return use_nox_observer_; }

    Teuchos::RCP<Tempus::IntegratorObserver<double> >
    buildTempusObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                      const Teuchos::RCP<const panzer::GlobalIndexer> & dof_manager,
                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      // note: Composite observer loops in the order added
      Teuchos::RCP<Tempus::IntegratorObserverComposite<double> > composite_observer =
        Teuchos::rcp(new Tempus::IntegratorObserverComposite<double>());

      TEUCHOS_ASSERT(nonnull(this->getMyParamList()));
      const Teuchos::ParameterList& p = *this->getMyParamList();

      // Add concrete observers below based on parameter list entries

      // clip first
      if (p.get<std::string>("Clip Solution Variables") == "ON") {
        std::vector<std::string> clip_dof_names;
        panzer::StringTokenizer(clip_dof_names,p.get<std::string>("Clipped Variable Names"),",",true);

        Teuchos::RCP<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > ep_lof
          = Teuchos::rcp_dynamic_cast<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >(lof,true);
        Teuchos::RCP<const panzer::GlobalIndexer> ep_dof_manager
          = Teuchos::rcp_dynamic_cast<const panzer::GlobalIndexer>(dof_manager,true);
        // Teuchos::RCP<charon::TempusObserver_ClipSolution> clip_solution_observer
        //   = Teuchos::rcp(new charon::TempusObserver_ClipSolution(ep_dof_manager,ep_lof,clip_dof_names));

        // composite_observer->addObserver(clip_solution_observer);
      }

      // Set some variables for output of responses
      bool output_responses = (p.get<std::string>("Output Responses") == "ON");
      bool output_exodus = (p.get<std::string>("Write Solution to Exodus File") == "ON");
      bool output_txt_response_file = (p.get<std::string>("Write Response File") != "");
      std::string txt_response_filename = p.get<std::string>("Write Response File");
      bool output_exo_responses = false;
      bool write_initial_condition = false;
      int write_after_this_many_steps = 1;
      if (output_exodus) {
        
        write_after_this_many_steps = p.get<int>("Time Step Interval for Writing Solution");
        
        if (p.get<std::string>("Write Initial Condition") == "TRUE")
          write_initial_condition = true;

        // If the user chose to output responses always output them to
        // exodus
        if (output_responses)
          output_exo_responses = true;
      }

      // Last sanity check, if there are no responses turn their output off
      if (!(parameterNames_.size()>0 || responseNames_.size()>0)) {
        output_responses = false;
        output_txt_response_file = false;
        output_exo_responses = false;
      }

      //-----------------------------------------------------------------------
      // Output data to various entities
      // ----------------------------------------------------------------------
      Teuchos::RCP<charon::TempusObserver_OutputData> response_observer =
        Teuchos::rcp(new charon::TempusObserver_OutputData(mesh,
                                                           dof_manager,
                                                           lof,
                                                           stkIOResponseLibrary_,
                                                           full_model_evaluator_,
                                                           responseNames_,
                                                           parameterNames_,
                                                           parameterLibrary_,
                                                           scaleParams_,
                                                           output_exodus,
                                                           output_responses,
                                                           output_txt_response_file,
                                                           output_exo_responses,
                                                           write_after_this_many_steps,
                                                           txt_response_filename,
                                                           write_initial_condition));
      
      composite_observer->addObserver(response_observer);

      //----------------------------------------------------------------------------------------------
      // Cluster observer
      //----------------------------------------------------------------------------------------------

      return composite_observer;
    }

    /** \name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
    {
      paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
      setMyParamList(paramList);

      use_nox_observer_ = (paramList->get<std::string>("Use NOX Observer") == "TRUE");
    }


    //--------------------------------------------------------------------------
    // Get valid parameters
    //--------------------------------------------------------------------------

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
    {
      if (valid_params_.is_null()) {

       valid_params_ = Teuchos::rcp(new Teuchos::ParameterList);


       //Write response file

       Teuchos::setStringToIntegralParameter<int>(
          "Write Solution to Exodus File",
          "ON",
          "Enables or disables writing of solution to Exodus file at end of NOX solve",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );

       // output responses. if ON then responses will be written to the
       // screen and optionally to a response file
       Teuchos::setStringToIntegralParameter<int>(
         "Output Responses",
         "OFF",
         "Should responses be output, always to the screen and possibly to an output file",
         Teuchos::tuple<std::string>("ON", "OFF"),
         valid_params_.get()
         );

       valid_params_->set<std::string>("Write Response File", "", "Writes table of responses to a file.");

       //clip solution

       Teuchos::RCP<Teuchos::EnhancedNumberValidator<int> > validator =
         Teuchos::rcp(new Teuchos::EnhancedNumberValidator<int>);
       validator->setMin(1);
       valid_params_->set<int>("Time Step Interval for Writing Solution", 1, "Writes solution to Exodus file after taking this many time steps.",validator);

       Teuchos::setStringToIntegralParameter<int>(
          "Write Initial Condition",
          "TRUE",
          "If set to TRUE, write the initial conditions to the Exodus file.",
          Teuchos::tuple<std::string>("TRUE","FALSE"),
          valid_params_.get()
          );

       Teuchos::setStringToIntegralParameter<int>(
          "Use NOX Observer",
          "FALSE",
          "If set to TRUE, then also install the NOX observer",
          Teuchos::tuple<std::string>("TRUE","FALSE"),
          valid_params_.get()
          );

       Teuchos::setStringToIntegralParameter<int>(
          "Clip Solution Variables",
          "OFF",
          "Enables or disables clipping of solution below 0",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );

        valid_params_->set<std::string>("Clipped Variable Names","","Variables names to clip");

        //Cluster observer

        valid_params_->set<std::string>("Cluster Observer","","");


      }
      return valid_params_;
    }

    void setModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > const full_me)
    {
      full_model_evaluator_ = full_me;
    }

    //@}

  };
}

#endif // CHARON_TEMPUS_OBSERVER_FACTORY
