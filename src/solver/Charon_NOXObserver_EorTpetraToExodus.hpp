
#ifndef CHARON_NOX_OBSERVER_EORTPETRA_TO_EXODUS_HPP
#define CHARON_NOX_OBSERVER_EORTPETRA_TO_EXODUS_HPP

// Charon
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_CurrentConstraintModelEvaluatorLOCA.hpp"

// Epetra
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

// LOCA
#include "LOCA_MultiContinuation_ArcLengthGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"

// NOX
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Epetra_Vector.H"

// Panzer
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"

#include "Panzer_TpetraLinearObjFactory_decl.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"

#include "Panzer_ResponseLibrary.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter.hpp"
#include "Panzer_GlobalIndexer.hpp"

// Piro
#include "Piro_ConfigDefs.hpp"
#include "Piro_NOXSolver.hpp"

// Teuchos
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"

// LOCA
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_Thyra_Group.H"

#include "Panzer_GlobalData.hpp"
#include "Charon_PanzerParameterExtractor.hpp"

namespace charon
{
  /**
   *  \brief Description.                                                        // JMG:  Fill this out.
   *                                                                             //
   *  Detailed description.                                                      //
   */
  class NOXObserver_EorTpetraToExodus
    :
    public NOX::Abstract::PrePostOperator
  {
    public:

      /**
       *  \brief Default Constructor.
       *
       *  Detailed description.                                                  // JMG:  Fill this out.
       *                                                                         //
       *  \param[in] mesh         Description.                                   //
       *  \param[in] lof          Description.                                   //
       *  \param[in] rLib         Description.                                   //
       *  \param[in] scaleFactors Description.                                   //
       *                                                                         //
       *  \throws x Description.                                                 //
       *                                                                         //
       *  \returns Something.                                                    //
       */
 
    //Add vector of global contact voltages to write to exodus
    std::map<std::string,double> contactVoltages_;
    std::map<std::string,double> norms_;
    Teuchos::RCP<panzerParameterExtractor>  pPE;
    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;
 

     NOXObserver_EorTpetraToExodus(
        const Teuchos::RCP<panzer_stk::STK_Interface>& mesh,
        const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits>>& lof,
        const Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>& rLib,
        const Teuchos::RCP<std::map<std::string, double>>& scaleFactors,
        bool writeOnSolveFail,
        const std::vector<std::string>& responseNames,
        const Teuchos::RCP<Thyra::ModelEvaluator<double>>& fullME,
	Teuchos::RCP<panzer::ParamLib> const& parameterLibrary,
        const bool& outputResponses)
        :
        parameterLibrary_(parameterLibrary),
        mesh_(mesh),
        lof_(lof),
        responseLibrary_(rLib),
        scaleFactors_(scaleFactors),
        writeOnSolveFail_(writeOnSolveFail),
        responseNames_(responseNames),
        fullME_(fullME),
        outputResponses_(outputResponses),
        stepCount_(0)
      {
        using panzer::Traits;
        using panzer_stk::RespFactorySolnWriter_Builder;
        using std::string;
        using std::vector;
        using Teuchos::rcp_dynamic_cast;
        using BELOF   =  panzer::BlockedEpetraLinearObjFactory<Traits, int>;
        using TLOF    = panzer::TpetraLinearObjFactory<Traits, double, int, panzer::GlobalOrdinal>;
        using MapIter = std::map<std::string, double>::const_iterator;

        TEUCHOS_ASSERT(not lof_.is_null())
        vector<string> eBlocks;
        mesh_->getElementBlockNames(eBlocks);
        RespFactorySolnWriter_Builder builder;
        builder.mesh = mesh_;
        for (MapIter i = (*scaleFactors_).begin(); i != (*scaleFactors_).end(); ++i)
        {
          builder.scaleField(i->first, i->second);
        }
        
        responseLibrary_->addResponse("Main Field Output", eBlocks, builder);
        
        useEpetra_ = not rcp_dynamic_cast<const BELOF>(lof_).is_null();
        useTpetra_ = not rcp_dynamic_cast<const TLOF>(lof_).is_null();

      } // end of Default Constructor

      /**
       *  \brief Things to be done before each solver iteration.
       *
       *  This routine doesn't do anything.
       *
       *  \param[in] solver This input isn't used.
       */
      void
      runPreIterate(
        const NOX::Solver::Generic& /* solver */)
      {
      } // end of runPreIterate()

      /**
       *  \brief Things to be done after each solver iteration.
       *
       *  This routine doesn't do anything.
       *
       *  \param[in] solver This input isn't used.
       */
      void
      runPostIterate(
        const NOX::Solver::Generic& /* solver */)
      {
      } // end of runPostIterate()

      /**
       *  \brief Things to be done before the solve.
       *
       *  This routine doesn't do anything.
       *
       *  \param[in] solver This input isn't used.
       */
      void
      runPreSolve(
        const NOX::Solver::Generic& /* solver */)
      {
      } // end of runPreSolve

      /**
       *  \brief Things to be done after the solve.
       *
       *  Detailed description.                                                  // JMG:  Fill this out.
       *                                                                         //
       *  \param[in] solver Description.                                         //
       *                                                                         //
       *  \throws x Description.                                                 //
       */
      void
      runPostSolve(
        const NOX::Solver::Generic& solver)
      {
        using LOCA::MultiContinuation::ArcLengthGroup;
        using LOCA::MultiContinuation::ExtendedVector;
        using LOCA::MultiContinuation::NaturalGroup;
        using panzer::AssemblyEngineInArgs;
        using panzer::EpetraLinearObjContainer;
        using panzer::TpetraLinearObjContainer;
        using panzer::Traits;
        using std::runtime_error;
        using Teuchos::RCP;
        using Teuchos::rcp_const_cast;
        using Teuchos::rcp_dynamic_cast;
        using Thyra::castOrCreateProductVectorBase;
        using Thyra::createMember;
        using Thyra::get_ele;
        using Thyra::ModelEvaluatorBase;
        using Thyra::ProductVectorBase;
        using Thyra::VectorBase;
        typedef NOX::Abstract::Group       AbsGroup;
        typedef NOX::Abstract::Vector      AbsVec;
        typedef NOX::Thyra::Group          ThGroup;
        typedef NOX::Thyra::Vector         ThVec;
        typedef panzer::LinearObjContainer LOC;

        TEUCHOS_ASSERT(not lof_.is_null())

        if (const_cast<NOX::Solver::Generic&>(solver).getStatus() == NOX::StatusTest::Failed &&
            !writeOnSolveFail_)
        {
          return;
        }

        // Attempt to cast the solver.getSolutionGroup() to a
        // NOX::Thyra::Group.  This should work if we're using NOX as the
        // solver.  We then cast it to a NOX::Abstract::Group because that's
        // the common ancestor between this and the other option below.
        const AbsGroup*       solnGroup(NULL);
        const NaturalGroup*   locaNGroup(NULL);
        const ArcLengthGroup* locaALGroup(NULL);
        const ThGroup*        thyraGroup =
          dynamic_cast<const ThGroup*>(&(solver.getSolutionGroup()));

        // Solution and linear system for LOCA Householder constraint
        // solver has no augmented rows so it can be used just like a
        // noraml nox thyra group. Set the values coreectly here.
        if (!thyraGroup) {
          const auto* locaMCCG = dynamic_cast<const LOCA::MultiContinuation::ConstrainedGroup*>(&(solver.getSolutionGroup()));
          if (locaMCCG) {
            auto locaThyraGroup = Teuchos::rcp_dynamic_cast<const LOCA::Thyra::Group>(locaMCCG->getUnderlyingGroup());
            TEUCHOS_ASSERT(nonnull(locaThyraGroup));
            thyraGroup = locaThyraGroup.get();
            solnGroup = locaThyraGroup.get();
          }
        }

        if (thyraGroup)
        {
          solnGroup = dynamic_cast<const AbsGroup*>(thyraGroup);

          // If our attempt at casting to what we needed ultimately failed,
          // we can't continue.
          TEUCHOS_TEST_FOR_EXCEPTION(solnGroup == NULL, runtime_error,
            "Failed to dynamic_cast &(solver.getSolutionGroup()) from " \
            "NOX::Thyra::Group* to NOX::Abstract::Group*.")
        }

        // If the cast to a NOX::Thyra::Group doesn't work, then let's try
        // casting it to a LOCA::MultiContinuation::NaturalGroup instead.  That
        // would be the case if we're using LOCA as the solver instead of NOX.
        // We then cast it to a NOX::Abstract::Group because that's the common
        // ancestor between this and the other option above.
        else
        {
          locaNGroup =
            dynamic_cast<const NaturalGroup*>(&(solver.getSolutionGroup()));
          if (locaNGroup == NULL)
          {
            locaALGroup = dynamic_cast<const ArcLengthGroup*>(
              &(solver.getSolutionGroup()));
            TEUCHOS_TEST_FOR_EXCEPTION(locaALGroup == NULL, runtime_error,
              "Failed to dynamic_cast &(solver.getSolutionGroup()) to "       \
              "either a LOCA::MultiContinuation::NaturalGroup* or "           \
              "ArcLengthGroup*.")
            solnGroup = dynamic_cast<const AbsGroup*>(locaALGroup);

            // If our attempt at casting to what we needed ultimately failed,
            // we can't continue.
            TEUCHOS_TEST_FOR_EXCEPTION(solnGroup == NULL, runtime_error,
              "Failed to dynamic_cast &(solver.getSolutionGroup()) from "     \
              "LOCA::MultiContinuation::ArcLengthGroup* to "                  \
              "NOX::Abstract::Group*.")
          }
          else // if (locaNGroup == NULL)
            solnGroup = dynamic_cast<const AbsGroup*>(locaNGroup);

          // If our attempt at casting to what we needed ultimately failed,
          // we can't continue.
          TEUCHOS_TEST_FOR_EXCEPTION(solnGroup == NULL, runtime_error,
            "Failed to dynamic_cast &(solver.getSolutionGroup()) from "       \
            "LOCA::MultiContinuation::NaturalGroup* to NOX::Abstract::Group*.")
        }

        // Grab the solution vector from the solution group and attempt to cast
        // it to a NOX::Thyra::Vector, which should work if NOX is your solver.
        const AbsVec& x    = solnGroup->getX();
        const ThVec*  nThX = dynamic_cast<const ThVec*>(&x);

        // If that cast doesn't work, try casting to a
        // LOCA::MultiContinuation::ExtendedVector instead, which should work
        // if you're using LOCA as the solver.
        if (not nThX)
        {
          const ExtendedVector* locaVec =
            dynamic_cast<const ExtendedVector*>(&x);
          TEUCHOS_TEST_FOR_EXCEPTION(locaVec == NULL, runtime_error,
            "Failed to dynamic_cast x from a NOX::Abstract::Vector to a " \
            "LOCA::MultiContinuation::ExtendedVector.")

          // Then you can get the solution vector out of the extended vector,
          // cast that to a NOX::Thyra::Vector, and be on your merry way.
          nThX = dynamic_cast<const ThVec*>(locaVec->getXVec().getRawPtr());

          // If using Housholder constrained loca continuation, the
          // solution is double nested. Have to go one layer deeper.
          if (not nThX) {
            const ExtendedVector* locaVecLevel2 =
              dynamic_cast<const ExtendedVector*>(locaVec->getXVec().getRawPtr());
            TEUCHOS_TEST_FOR_EXCEPTION(locaVecLevel2 == NULL, runtime_error,
                                       "Failed to dynamic_cast x from a NOX::Abstract::Vector to a " \
                                       "LOCA::MultiContinuation::ExtendedVector.")
            nThX = dynamic_cast<const ThVec*>(locaVecLevel2->getXVec().getRawPtr());
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(nThX == NULL, runtime_error,
                                   "Failed to dynamic_cast to NOX::Thyra::Vector!");
        
        RCP<const VectorBase<double>> thX = nThX->getThyraRCPVector();
        RCP<const VectorBase<double>> origThX = thX;

        // If we're using a blocked system (for instance, for the constant
        // current or resistor contact cases), we need to extract only the PDE
        // degrees of freedom, which are contained in the first vector block.
        const RCP<const ProductVectorBase<double>> bloX =
          castOrCreateProductVectorBase(thX);
        if (bloX->productSpace()->numBlocks() > 1)
          thX = bloX->getVectorBlock(0);

        // Initialize the assembly container.
        AssemblyEngineInArgs aeInArgs;
        aeInArgs.container_ = lof_->buildLinearObjContainer();
        aeInArgs.ghostedContainer_ = lof_->buildGhostedLinearObjContainer();
        aeInArgs.alpha = 0.0;
        aeInArgs.beta = 1.0;
        aeInArgs.evaluate_transient_terms = false;

        // Initialize the ghosted container.
        lof_->initializeGhostedContainer(LOC::X, *aeInArgs.ghostedContainer_);

        if (useEpetra_)
        {
           // Initialize the x vector.
           const RCP<EpetraLinearObjContainer> epGlobalContainer =
             rcp_dynamic_cast<EpetraLinearObjContainer>(aeInArgs.container_,
             true);
           epGlobalContainer->set_x_th(
             rcp_const_cast<VectorBase<double>>(thX));
        }
        else if (useTpetra_)
        {
          const RCP<TpetraLinearObjContainer<double, int, panzer::GlobalOrdinal> > tpGlobalContainer =
            rcp_dynamic_cast<TpetraLinearObjContainer<double, int, panzer::GlobalOrdinal> >(aeInArgs.container_,
                                                       true);
          tpGlobalContainer->set_x_th(rcp_const_cast<VectorBase<double> >(thX));
        }
        else // Unknown
        {
           std::ostringstream err_msg("");
           err_msg << "Unknown container in " << __FILE__ << " at line " << __LINE__;
           TEUCHOS_TEST_FOR_EXCEPTION(true, runtime_error, err_msg.str());
          
        }

        responseLibrary_->addResponsesToInArgs<Traits::Residual>(aeInArgs);
        responseLibrary_->evaluate<Traits::Residual>(aeInArgs);

        // Write to disk.
        if (locaNGroup)
          mesh_->addGlobalToExodus("ContinuationParameter",
            locaNGroup->getContinuationParameter());
        else if (locaALGroup)
          mesh_->addGlobalToExodus("ContinuationParameter",
            locaALGroup->getContinuationParameter());
	//Add contact voltages to exodus
	//Get contact voltage names and values
	pPE = Teuchos::rcp(new panzerParameterExtractor(parameterLibrary_));
	contactVoltages_ = pPE->get_VoltageParameters();
	std::map<std::string,double>::iterator cVit = contactVoltages_.begin();
	while (cVit != contactVoltages_.end())
	  {
	    mesh_->addGlobalToExodus(cVit->first, cVit->second);
	    ++cVit;
	  }

  //Add norms to exodus
	//Get norm names and values
	//pPE = Teuchos::rcp(new panzerParameterExtractor(parameterLibrary_));
	norms_ = pPE->get_NormParameters();
	std::map<std::string,double>::iterator norm_it = norms_.begin();
	while (norm_it != norms_.end())
	  {
	    mesh_->addGlobalToExodus(norm_it->first, norm_it->second);
	    ++norm_it;
	  }


        if (outputResponses_)
        {
          // LOCA Householder model evalautor returns residual
          // equations in place of response current values. To write
          // the correct responses, we need to pull out the underlying
          // model evalautor.
          auto tmpME = fullME_;
          auto checkME = Teuchos::rcp_dynamic_cast<charon::CurrentConstraintModelEvaluatorLOCA<double>>(fullME_);
          if (nonnull(checkME))
            tmpME = checkME->getInternalPhysicsME();

          ModelEvaluatorBase::InArgs<double> inArgs = tmpME->createInArgs();
          inArgs.set_x(origThX);
          ModelEvaluatorBase::OutArgs<double> outArgs = tmpME->createOutArgs();
          int Ng(outArgs.Ng());
          for (int i(0); i < Ng; ++i)
            outArgs.set_g(i, createMember(*(tmpME->get_g_space(i))));
          tmpME->evalModel(inArgs, outArgs);
          for (int i(0); i < Ng; ++i)
          {
            RCP<VectorBase<double>> response = outArgs.get_g(i);
            mesh_->addGlobalToExodus(responseNames_[i], get_ele(*response, 0));
          } // end loop over the responses
        } // end if (outputResponses_)
        double fake_time = stepCount_++;
        mesh_->writeToExodus(fake_time);
      } // end of runPostSolve()

    protected:

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       *                                                                         //
       *  \param[in, out] os  Description.                                       //
       *  \param[in]      src Description.                                       //
       */
      void
      writeToScreen(
        std::ostream&                    os,
        const Thyra::VectorBase<double>& src)
      {
        using std::endl;
        using Teuchos::ArrayRCP;
        using Teuchos::dyn_cast;
        using Teuchos::ptrFromRef;
        using Thyra::SpmdVectorBase;
        const SpmdVectorBase<double>& spmdSrc =
          dyn_cast<const SpmdVectorBase<double>>(src);

        // Get access to the data.
        ArrayRCP<const double> srcData;
        spmdSrc.getLocalData(ptrFromRef(srcData));
        os << "Local Size = " << srcData.size() << endl;
        for (int i(0); i < srcData.size(); ++i)
           os << "   " << srcData[i] << endl;
      } // end of writeToScreen()

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       */
      Teuchos::RCP<panzer_stk::STK_Interface> mesh_;

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       */
      Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits>> lof_;

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       */
      Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> responseLibrary_;

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       */
      Teuchos::RCP<std::map<std::string, double>> const& scaleFactors_;

      /**
       * \brief Should output be done even in the event of a solver failure
       */
      bool writeOnSolveFail_;

      /**
       *  \brief A list of the names of all the responses.
       */
      std::vector<std::string> responseNames_;

      /**
       *  \brief The `ModelEvaluator` used to set up and solve the system.
       */
      Teuchos::RCP<Thyra::ModelEvaluator<double>> fullME_;

      /**
       *  \brief Whether or not we should output the current responses to the
       *         Exodus file.
       */
      bool outputResponses_;

      /**
       * \brief keep track of the step count for LOCA so that it can be used as the value for "time" in the case of a continuation run.
       */
      size_t stepCount_;

      /**
       *  \brief Description.                                                    // JMG:  Fill this out.
       *                                                                         //
       *  Detailed description.                                                  //
       */
      bool useEpetra_;

      bool useTpetra_;

  }; // end of class NOXObserver_EorTpetraToExodus
} // end of namespace charon

#endif
