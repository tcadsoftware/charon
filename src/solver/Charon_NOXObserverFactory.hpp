
#ifndef CHARON_NOX_OBSERVER_FACTORY_HPP
#define CHARON_NOX_OBSERVER_FACTORY_HPP

#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include "NOX_PrePostOperator_Vector.H"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Assert.hpp"

#include "Charon_Scaling_Parameters.hpp"

// Concrete nox observers
#include "Charon_NOXObserver_EorTpetraToExodus.hpp"
#include "Charon_NOXObserver_EorTpetraOutput.hpp"
#include "Charon_NOXObserver_WriteResponses.hpp"
#include "Charon_PanzerParameterExtractor.hpp"

namespace charon {

  class NOXObserverFactory :
    public panzer_stk::NOXObserverFactory,
    public Teuchos::ParameterListAcceptorDefaultBase {

  private:

    //! Store STK IO response library...be careful, it will be modified externally
    Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > stkIOResponseLibrary_;

    mutable Teuchos::RCP<Teuchos::ParameterList> valid_params_;

    Teuchos::RCP<std::map<std::string,double> > const& scaleFactors_;

    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;

    //! Store Response names for pretty printing of currents, etc...
    std::vector<std::string> responseNames_;

    bool const isLOCASolver_;

    bool const writeOnSolveFail_;

    Teuchos::RCP<Thyra::ModelEvaluator<double> > fullModelEvaluator_;

  public:

    NOXObserverFactory(const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits> > & stkIOResponseLibrary,
                       Teuchos::RCP<std::map<std::string,double> > const& scaleFactors,
		       Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
                       bool isLOCASolver,
                       bool writeOnSolveFail)
      : stkIOResponseLibrary_(stkIOResponseLibrary),
        scaleFactors_(scaleFactors),
	parameterLibrary_(parameterLibrary),
        isLOCASolver_(isLOCASolver),
        writeOnSolveFail_(writeOnSolveFail) {}


    Teuchos::RCP<NOX::Abstract::PrePostOperator>
    buildNOXObserver(const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
                     const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
                     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> >& lof) const
    {
      using std::string;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::rcp_dynamic_cast;

      TEUCHOS_ASSERT( nonnull(this->getParameterList()) );

      RCP<const panzer::GlobalIndexer> epetraDOFManager = rcp_dynamic_cast<const panzer::GlobalIndexer>(dof_manager);

      RCP<NOX::PrePostOperatorVector> observer = rcp(new NOX::PrePostOperatorVector);

      // solution output has been canceled
      if (this->getParameterList()->get<std::string>("Write Solution to Exodus File") == "ON") {
        Teuchos::RCP<NOX::Abstract::PrePostOperator> solution_writer =
          Teuchos::rcp(new charon::NOXObserver_EorTpetraToExodus(mesh,
                                                                 lof,
                                                                 stkIOResponseLibrary_,
                                                                 scaleFactors_,
                                                                 writeOnSolveFail_,
                                                                 responseNames_,
                                                                 fullModelEvaluator_,
								 parameterLibrary_,
                                                                this->getParameterList()->get<bool>("Output Responses")));
        observer->pushBack(solution_writer);
      }

      // solution output has been canceled
      if (this->getParameterList()->get<std::string>("Write Linear System") == "ON") {
        TEUCHOS_ASSERT(epetraDOFManager!=Teuchos::null);

        Teuchos::RCP<NOX::Abstract::PrePostOperator> linear_sys_writer
          = rcp(new charon::NOXObserver_EorTpetraOutput());

        observer->pushBack(linear_sys_writer);
      }

      // Output any responses. Presently this is just the scalar
      // electric current, but others can be added in the future.
      if (this->getParameterList()->get<bool>("Output Responses"))
      {
        string responseFileName("");
        if (this->getParameterList()->isParameter("Output Responses File"))
          responseFileName =
            this->getParameterList()->get<string>("Output Responses File");
        Teuchos::RCP<NOX::Abstract::PrePostOperator> response_writer =
          Teuchos::rcp(new charon::NOXObserver_WriteResponses(responseNames_,
                                                              fullModelEvaluator_,
                                                              true,
                                                              true,
                                                              responseFileName,
							      parameterLibrary_,
                                                              isLOCASolver_,
                                                              writeOnSolveFail_));
        observer->pushBack(response_writer);
      }

      return observer;
    }

    void setResponseNames(std::vector<std::string> const& responseNames)
    {
      responseNames_ = responseNames;
    }

    /**
     * \brief pass in the model evaluator to ease/facilitate access
     */
    void setModelEvaluator(Teuchos::RCP<Thyra::ModelEvaluator<double> > const full_me)
    {
      fullModelEvaluator_ = full_me;
    }

    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
    {
      using Teuchos::RCP;
      using Teuchos::rcp;

      paramList->validateParametersAndSetDefaults(*(this->getValidParameters()));
      setMyParamList(paramList);
    }

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
    {
      using std::string;

      if (valid_params_.is_null())
      {

        valid_params_ = Teuchos::rcp(new Teuchos::ParameterList);

        Teuchos::setStringToIntegralParameter<int>(
          "Write Solution to Exodus File",
          "ON",
          "Enables or disables writing of solution to Exodus file at end of NOX solve",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );

        Teuchos::setStringToIntegralParameter<int>(
          "Write Linear System",
          "OFF",
          "Enables or disables writing of linear system to matrix market file at the end of the NOX solve",
          Teuchos::tuple<std::string>("ON","OFF"),
          valid_params_.get()
          );

        valid_params_->set<bool>("Output Responses", false);
        valid_params_->set<string>("Output Responses File", "currents-loca.dat");

        valid_params_->set<bool>("Output Responses to File", false);

      }
      return valid_params_;
    }

  };

}

#endif
