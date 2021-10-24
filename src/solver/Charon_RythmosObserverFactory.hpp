
#ifndef CHARON_RYTHMOS_OBSERVER_FACTORY_HPP
#define CHARON_RYTHMOS_OBSERVER_FACTORY_HPP

#include <Charon_config.hpp>

#include "Panzer_STK_RythmosObserverFactory.hpp"
#include "Rythmos_CompositeIntegrationObserver.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_String_Utilities.hpp"


// Concrete rythmos observers
#include "Charon_RythmosObserver_WriteToExodus.hpp"
#include "Charon_RythmosObserver_WriteResponses.hpp"
#include "Charon_RythmosObserver_ClipSolution.hpp"
#include "Charon_Scaling_Parameters.hpp"


namespace charon {

  class RythmosObserverFactory : public panzer_stk::RythmosObserverFactory,
                             public Teuchos::ParameterListAcceptorDefaultBase  {

  private:

    mutable Teuchos::RCP<Teuchos::ParameterList> valid_params_;
    bool use_nox_observer_;

    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    //! Store Response names for pretty printing of currents, etc...
    std::vector<std::string> responseNames_;

    //! Store parameter names for pretty printing of voltages, etc...
    std::vector<std::string> parameterNames_;

    Teuchos::RCP<charon::Scaling_Parameters> scaleParams_;



  public:
    RythmosObserverFactory() {}

    RythmosObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                           const std::vector<std::string> & responseNames,
                           const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams)
       : stkIOResponseLibrary_(stkIOResponseLibrary), responseNames_(responseNames), scaleParams_(scaleParams) {}

    RythmosObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                           const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams)
       : stkIOResponseLibrary_(stkIOResponseLibrary), scaleParams_(scaleParams) {}


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

    Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
    buildRythmosObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                      const Teuchos::RCP<const panzer::GlobalIndexer> & dof_manager,
                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_dynamic_cast;

      // note: Composite observer loops in the order added
      RCP<Rythmos::CompositeIntegrationObserver<double> > composite_observer =
          Rythmos::createCompositeIntegrationObserver<double>();

      TEUCHOS_ASSERT(nonnull(this->getMyParamList()));
      const Teuchos::ParameterList& p = *this->getMyParamList();

      // Add concrete observers below based on parameter list entries

      // clip first
      if (p.get<std::string>("Clip Solution Variables") == "ON") {
        std::vector<std::string> clip_dof_names;
        panzer::StringTokenizer(clip_dof_names,p.get<std::string>("Clipped Variable Names"),",",true);

        RCP<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > ep_lof
          = rcp_dynamic_cast<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >(lof,true);
        RCP<const panzer::GlobalIndexer> ep_dof_manager
          = rcp_dynamic_cast<const panzer::GlobalIndexer>(dof_manager,true);
        RCP<charon::RythmosObserver_ClipSolution> clip_solution_observer
          = rcp(new charon::RythmosObserver_ClipSolution(ep_dof_manager,ep_lof,clip_dof_names));

        composite_observer->addObserver(clip_solution_observer);
      }

      // write solution last
      if (p.get<std::string>("Write Solution to Exodus File") == "ON") {
        int write_after_this_many_steps = p.get<int>("Time Step Interval for Writing Solution");
        bool write_initial_condition = false;

        if (p.get<std::string>("Write Initial Condition") == "TRUE")
          write_initial_condition = true;

        RCP<charon::RythmosObserver_WriteToExodus> exodus_observer =
          rcp(new charon::RythmosObserver_WriteToExodus(mesh,dof_manager,lof,write_after_this_many_steps,write_initial_condition,stkIOResponseLibrary_,scaleParams_));

        composite_observer->addObserver(exodus_observer);
      }

      //----------------------------------------------------------------------------------------------
      // write out responses (only if you have a response names vector and the user requests responses
      //----------------------------------------------------------------------------------------------
      if(   (parameterNames_.size()>0 || responseNames_.size()>0)
         && (p.get<std::string>("Output Responses") == "ON"
         || p.get<std::string>("Write Response File")!="")) {
        RCP<charon::RythmosObserver_WriteResponses> response_observer =
          rcp(new charon::RythmosObserver_WriteResponses(responseNames_,parameterNames_,scaleParams_,
                                                         p.get<std::string>("Output Responses")=="ON",
                                                         p.get<std::string>("Write Response File")));

        composite_observer->addObserver(response_observer);
      }


      //----------------------------------------------------------------------------------------------
      // Cluster observer
      //----------------------------------------------------------------------------------------------

      return composite_observer;
    }

    /** \name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;

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

       Teuchos::setStringToIntegralParameter<int>(
          "Output Responses",
          "ON",
          "Enables or disables writing of responses (like the current) to the screen.",
          Teuchos::tuple<std::string>("ON","OFF"),
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

    //@}

  };

}

#endif
