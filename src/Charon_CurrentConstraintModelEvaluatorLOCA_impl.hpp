
///////////////////////////////////////////////////////////////////////////////
//
//  Charon_CurrentConstraintModelEvaluator_impl.hpp
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __Charon_CurrentConstraintModelEvaluator_impl_hpp__
#define __Charon_CurrentConstraintModelEvaluator_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Teuchos
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

// Thyra
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultSpmdVector.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_MultiVectorStdOps_def.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_VectorStdOps_def.hpp"

// Mixed mode requirements
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

#include "NOX_TpetraTypedefs.hpp" // For underlying types for LOCA Bordering
#include "Thyra_TpetraThyraWrappers.hpp" // Convert Thyra to Tpetra
#include "Thyra_DefaultSpmdVector_decl.hpp" // Convert Thyra to Tpetra

namespace charon
{
  // *********************************************************
  template <typename Scalar>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::
  CurrentConstraintModelEvaluatorLOCA(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
    MPI_Comm                                           comm,
    const charon::CurrentConstraintList&               constraints,
    const int&                                         dimension,
    const Teuchos::RCP<Teuchos::ParameterList>&        parameters,
    const bool                                         print_debug)
    :
    Thyra::ModelEvaluatorDelegatorBase<Scalar>(physics),
    physics_(physics),
    comm_(Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(comm)))),
    constraints_(constraints),
    dimension_(dimension),
    parameters_(parameters),
    print_debug_(print_debug)
  {
    TEUCHOS_ASSERT(nonnull(physics_));

    // both p and g are single locally replicated double values
    auto pMap = Teuchos::rcp(new const NOX::TMap(1, 0, comm_, Tpetra::LocallyReplicated));
    pSpace_ = ::Thyra::createVectorSpace<Scalar, NOX::LocalOrdinal, NOX::GlobalOrdinal, NOX::NodeType>(pMap);
    gSpace_ = pSpace_;

    // We need to swap the Thyra::Spmd to Tpetra::MultiVector for p,
    // g, and DgDp in nominal values coming from underlying physics
    // model.
    nominalValues_ = physics_->getNominalValues();
    for (int i=0; i < constraints_.size(); ++i) {
      using tpetra_extract = ::Thyra::TpetraOperatorVectorExtraction<Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
      int pIdx = constraints_[i]->parameterIndex();
      auto new_p = Thyra::createMember(pSpace_);
      auto old_p = nominalValues_.get_p(pIdx);
      Scalar val = Thyra::get_ele(*old_p,0);
      auto new_p_tpetra = tpetra_extract::getTpetraMultiVector(new_p);
      auto new_p_host = new_p_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
      new_p_host(0,0) = val;
      nominalValues_.set_p(pIdx,new_p);
    }
  }

  // *********************************************************
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluatorLOCA<Scalar>::
  assignValueTpetraToSpmd(const Teuchos::RCP<const ::Thyra::VectorBase<Scalar>>& source_tpetra_value,
                          const Teuchos::RCP<::Thyra::VectorBase<Scalar>>& target_simd_value) const
  {
    using tpetra_extract = ::Thyra::TpetraOperatorVectorExtraction<Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
    auto src_tpetra = tpetra_extract::getConstTpetraMultiVector(source_tpetra_value);
    auto src_tpetra_host = src_tpetra->getLocalViewHost(Tpetra::Access::ReadOnly);
    Thyra::assign(target_simd_value.ptr(),src_tpetra_host(0,0));

    Scalar val = Thyra::get_ele(*target_simd_value,0);

    if (print_debug_)
      std::cout << "CurrentConstraintModelEvalautorLOCA::assignValueTpetraToSpmd p=" << std::setprecision(10) << val << std::endl;
  }

  // *********************************************************
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluatorLOCA<Scalar>::
  assignValueSpmdToTpetra(const Teuchos::RCP<const ::Thyra::MultiVectorBase<Scalar>>& source_simd_value,
                          const Teuchos::RCP<::Thyra::MultiVectorBase<Scalar>>& target_tpetra_value) const
  {
    auto src_spmd_vector = Teuchos::rcp_dynamic_cast<const ::Thyra::DefaultSpmdVector<Scalar>>(source_simd_value,true);
    Teuchos::ArrayRCP<const Scalar> src_vals;
    src_spmd_vector->getLocalData(Teuchos::ptr(&src_vals));

    using tpetra_extract = ::Thyra::TpetraOperatorVectorExtraction<Scalar,NOX::LocalOrdinal,NOX::GlobalOrdinal,NOX::NodeType>;
    auto target_tpetra = tpetra_extract::getTpetraMultiVector(target_tpetra_value);
    auto target_tpetra_host = target_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
    target_tpetra_host(0,0) = src_vals[0];

    if (print_debug_)
      std::cout << "CurrentConstraintModelEvaluatorLOCA::assignValueSpmdToTpetra g=" << std::setprecision(10) << src_vals[0] << std::endl;
  }

  // *********************************************************
  template <typename Scalar>
  Teuchos::RCP<::Thyra::ModelEvaluator<Scalar>>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::getInternalPhysicsME() const
  {
    return physics_;
  }

  // *********************************************************
  template <typename Scalar>
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::get_p_space(int l) const
  {
    for (int i=0; i < constraints_.size(); ++i)
      if (constraints_[i]->parameterIndex() ==  l)
        return pSpace_;

    // If not a constraint, just pass through to underlying ME.
    return physics_->get_p_space(l);
  }

  // *********************************************************
  template <typename Scalar>
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::get_g_space(int j) const
  {
    for (int i=0; i < constraints_.size(); ++i)
      if (constraints_[i]->responseIndex() ==  j)
        return gSpace_;

    // If not a constraint, just pass through to underlying ME.
    return physics_->get_g_space(j);
  }

  // *********************************************************
  template<typename Scalar>
  std::string
  CurrentConstraintModelEvaluatorLOCA<Scalar>::description() const
  {
    return physics_->description();
  }

  // *********************************************************
  template<typename Scalar>
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::createInArgs() const
  {
    return this->getNominalValues();
  }

  // *********************************************************
  template <typename Scalar>
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::getNominalValues() const
  {
    return nominalValues_;
  }

  // *********************************************************
  template <typename Scalar>
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  CurrentConstraintModelEvaluatorLOCA<Scalar>::create_DgDp_op( int j, int l ) const
  {
    return ::Thyra::createMembers(gSpace_,1,"CurrentConstraint:DgDp");
  }

  // *********************************************************
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluatorLOCA<Scalar>::
  evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
  {
    using LO = NOX::LocalOrdinal;
    using GO = NOX::GlobalOrdinal;
    using NT = NOX::NodeType;
    using tpetra_extract = ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,NT>;

    // Create a copy of both InArgs and OutArgs to pass to the
    // underlying physics model evaluator. Panzer uses raw Thyra::Spmd
    // vectors for p, g, and dgdp. LOCA requires Tpetra multivectors
    // for p, g, and dgdp. This ME wrapper translates the parameters
    // between LOCA and Panzer. Long term, Panzer should just use the
    // Tpetra version and we should dump Spmd.
    typename ::Thyra::ModelEvaluatorBase::InArgs<Scalar> physicsInArgs(inArgs);
    typename ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> physicsOutArgs(outArgs);

    for (int i=0; i < constraints_.size(); ++i) {
      auto& c = *constraints_[i];

      // Replace constraint parameters (Tpetra MultiVectors) with
      // underlying panzer parameter type (Thyra::SpmdVector).
      auto p_index = c.parameterIndex();
      auto p = inArgs.get_p(p_index);
      if (nonnull(p)) {
        auto physics_p = ::Thyra::createMember(physics_->get_p_space(p_index),"physics_p");
        this->assignValueTpetraToSpmd(p,physics_p);
        physicsInArgs.set_p(p_index,physics_p);
      }

      // Replace g and dgdp responses (Tpetra MultiVectors) with
      // underlying panzer response type (Thyra::SpmdVector).
      auto g_index = c.responseIndex();
      auto g = outArgs.get_g(g_index);

      if (nonnull(g)) {
        auto physics_g = ::Thyra::createMember(physics_->get_g_space(g_index),"physics_g");
        // NOTE: may need to initialize to incoming value
        ::Thyra::assign(physics_g.ptr(),0.0);
        physicsOutArgs.set_g(g_index,::Thyra::ModelEvaluatorBase::Evaluation<::Thyra::VectorBase<double>>(physics_g));
      }
      for (int j=0; j < constraints_.size(); ++j) {
        auto dgdp = outArgs.get_DgDp(g_index,constraints_[j]->parameterIndex()).getMultiVector();
        if (nonnull(dgdp)) {
          // NOTE: we fill in dgdp in this ME, don't pass on dgdp to
          // underlying physics ME. Rather set a dummy object.
          ::Thyra::ModelEvaluatorBase::Derivative<Scalar> dummy;
          physicsOutArgs.set_DgDp(g_index,constraints_[j]->parameterIndex(),dummy);
        }
        // auto physicsDgDp = physics_->create_DgDp_op(g_index,p_index);
        // auto physicsDgDpMV = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(physicsDgDp);
        // // NOTE: may need to initialize to incoming value
        // ::Thyra::assign(physicsDgDpMV.ptr(),0.0);
        // physicsOutArgs.set_DgDp(g_index,p_index,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(physicsDgDpMV,::Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM));
      }
    }

    // Evalaute the underlying model
    physics_->evalModel(physicsInArgs,physicsOutArgs);

    // Get the current response values and transform them in place
    // into the constraint equations for LOCA Householder Bordering
    // constraints. If this works, long term we should just add the
    // extra parameter separately to support both constraints and
    // currents at the same time. It's extra work since we will have
    // to support a new outArgs and swap arguments in outArgs.

    // Loop over responses and update values/derivatives if requested
    for (int i=0; i < constraints_.size(); ++i) {

      auto& c = *constraints_[i];

      auto p_index = c.parameterIndex();
      auto p = inArgs.get_p(p_index);
      if (p.is_null())
        p = this->getNominalValues().get_p(p_index);

      auto g_index = c.responseIndex();
      auto g = outArgs.get_g(g_index);
      auto dgdx = outArgs.get_DgDx(g_index).getMultiVector();

      // Copy internal physics ME Spmd to Tpetra
      if (nonnull(g)) {
        auto g_spmd = physicsOutArgs.get_g(g_index);
        this->assignValueSpmdToTpetra(g_spmd,g);
      }
      for (int j=0; j < constraints_.size(); ++j) {
        auto dgdp = outArgs.get_DgDp(g_index,constraints_[j]->parameterIndex()).getMultiVector();
        if (nonnull(dgdp)) {
          // auto dgdp_spmd = physicsOutArgs.get_DgDp(g_index,p_index).getMultiVector();
          // this->assignValueSpmdToTpetra(dgdp_spmd,dgdp);
        }
      }

      if (constraints_[i]->isConstantCurrent()) {
        // If we have a constant current case, then in 3-D our
        // constraint equation is relatively straightforward: I =
        // Ia, where I is the current at the contact, Ia is the
        // applied current, and both I and Ia are in Amps.  We can
        // rewrite that as g = I - Ia = 0.  Our parameter is the
        // voltage, that is, p = V, so dg/dp = 0.  In 2-D the
        // current we're working with is actually in Amps per
        // centimeter, so we need to modify our constraint equation
        // using the simulation contact length, L, and the device
        // contact area, A, which have units of micrometers and
        // square micrometers, respectively.  The constraint
        // equation is therefore (I * A / L * 1e-4) = Ia, which we
        // can rewrite as g = I - Ia * (L / A * 1e4) = 0, where the
        // factor of 1e4 comes from the conversion between
        // micrometers and centimeters.  Our parameter is still the
        // voltage, so dg/dp = 0 in 2-D as well.

        // Evaluate g
        if (nonnull(g)) {
          // after evalModel(), g contains I
          // transform to g = I - Ia * (L / A * 1e4)
          // g should be locally replicated scalar
          auto g_tpetra = tpetra_extract::getTpetraMultiVector(g);
          TEUCHOS_ASSERT(nonnull(g_tpetra));

          Scalar adjustment = -c.currentValue();
          if (dimension_ == 2)
            adjustment *= c.contactLength() / c.contactArea() * 1.0e4;

          auto g_host = g_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
          g_host(0,0) += adjustment;
        }

        // Evaluate DgDp
        for (int j=0; j < constraints_.size(); ++j) {
          auto dgdp = outArgs.get_DgDp(g_index,constraints_[j]->parameterIndex()).getMultiVector();
          if (nonnull(dgdp)) {
            auto dgdp_tpetra = tpetra_extract::getTpetraMultiVector(dgdp);

            if(parameters_->isParameter("Mixed Mode")){
              Teuchos::RCP<panzer::ParamLib> paramLib = 
                (parameters_->get<Teuchos::RCP<panzer::GlobalData > >("Global Data"))->pl;
	      double xyceSensitivity = paramLib->getValue<panzer::Traits::Residual>("Xyce Coupling Sensitivity"+std::to_string(j));
	      // for mixed mode, the contact current is scaled; this may be useful debugging info
              //std::cout << "Mixed Mode contact current scaling in dg/dp: " << c.contactLength() / c.contactArea() * 1.0e4 << std::endl;
              double scaleFactor = 1.0;
	      if (dimension_ == 2)
	        scaleFactor = c.contactLength() / c.contactArea() * 1.0e4;
              dgdp_tpetra->putScalar(-1.0*xyceSensitivity*scaleFactor);
	      // check this sign so that it is consistent with in/out-flowing current
	      // is there supposed to be a second -1.0 sign? dg/dp = d(Icharon - Ia*scaleFactor )/dp
	    }
	    else
	      dgdp_tpetra->putScalar(0.0);
          }
        }

        // Evaluate DgDx
        if (nonnull(dgdx)) {
          // no adjustment needed. Just pass through
        }
      }
      else { // if (constraints_[i]->isResistorContact())
        // If we have a resistor contact case, then in 3-D our constraint
        // equation is V - Va = I * R, where V is the voltage at the contact,
        // Va is the applied voltage, both in Volts, I is the current at the
        // contact in Amps, and R is the resistance in Ohms.  We can rewrite
        // that as g = I - (V - Va) / R = 0.  Taking the derivative with
        // respect to voltage again, we have dg/dp = -1 / R.
        //   In 2-D the current is in Amps per centimeter again, so we need
        // to modify the equation as above:  V - Va = (I * A / L * 1e-4) * R,
        // which we can rewrite as g = I - (V - Va) / R * (L / A * 1e4) = 0.
        // The derivative with respect to voltage is therefore
        // dg/dp = (-1 / R) * (L / A * 1e4).

        // value = -1.0 / constraints_[i]->resistorValue();
        // if (dimension_ == 2)
        //   value *= constraints_[i]->contactLength() /
        //     constraints_[i]->contactArea() * 1e4;

        // need parameter for g evaluation
        auto p_tpetra = tpetra_extract::getConstTpetraMultiVector(p);
        TEUCHOS_ASSERT(nonnull(p_tpetra));
        Scalar p_val = (p_tpetra->getLocalViewHost(Tpetra::Access::ReadOnly))(0,0);

        // Evaluate g
        if (nonnull(g)) {
          // after evalModel(), g contains I
          // transform to g = I - (V - Va) / R
          // g should be locally replicated scalar
          auto g_tpetra = tpetra_extract::getTpetraMultiVector(g);
          TEUCHOS_ASSERT(nonnull(g_tpetra));

          Scalar adjustment = -1.0 * (p_val - c.appliedVoltage()) / c.resistorValue();
          if (dimension_ == 2)
            adjustment *= c.contactLength() / c.contactArea() * 1.0e4;

          auto g_host = g_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
          g_host(0,0) += adjustment;
        }

        // Evaluate DgDp
        for (int j=0; j < constraints_.size(); ++j) {
          auto dgdp = outArgs.get_DgDp(g_index,constraints_[j]->parameterIndex()).getMultiVector();
          if (nonnull(dgdp)) {
            auto dgdp_tpetra = tpetra_extract::getTpetraMultiVector(dgdp);
            auto dgdp_host = dgdp_tpetra->getLocalViewHost(Tpetra::Access::ReadWrite);
            if (i == j) {
              dgdp_host(0,0) = -1.0 / c.resistorValue();
              if (dimension_ == 2)
                dgdp_host(0,0) *= c.contactLength() / c.contactArea() * 1.0e4;
            }
            else {
              dgdp_host(0,0) = 0.0;
            }
          }
        }

        // Evaluate DgDx
        if (nonnull(dgdx)) {
          // no adjustment needed. Just pass through
        }

      }
    }
  }

} // end of namespace charon

#endif // __Charon_CurrentConstraintModelEvaluatorLOCA_impl_hpp__
