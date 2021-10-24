
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

namespace charon
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Physics Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  CurrentConstraintModelEvaluator<Scalar>::
  CurrentConstraintModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
    MPI_Comm                                           comm,
    const charon::CurrentConstraintList&               constraints,
    const int&                                         dimension)
    :
    Thyra::ModelEvaluatorDelegatorBase<Scalar>(physics),
    physics_(physics),
    comm_(Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(
      Teuchos::opaqueWrapper(comm)))),
    constraints_(constraints),
    dimension_(dimension)
  {
    // Ensure we were actually given a physics ModelEvaluator.
    TEUCHOS_ASSERT(not physics_.is_null())

    // Build the necessary vector spaces.
    initialize();
  } // end of Physics Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  get_x_space()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  get_x_space() const
  {
    return xSpace_;
  } // end of get_x_space()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  get_f_space()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  get_f_space() const
  {
    return fSpace_;
  } // end of get_f_space()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  create_W_op()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  create_W_op() const
  {
    using std::size_t;
    using std::vector;
    using Teuchos::Comm;
    using Teuchos::DefaultComm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::DefaultSpmdVectorSpace;
    using Thyra::defaultSpmdVectorSpace;
    using Thyra::LinearOpBase;
    using Thyra::createMembers;
    using Thyra::nonconstBlock2x2;
    using Thyra::nonconstTranspose;
    using Thyra::productVectorSpace;
    using Thyra::VectorSpaceBase;

    // Grab the tangent stiffness matrix for the unconstrained system (df/du).
    RCP<LinearOpBase<Scalar>> A00  = physics_->create_W_op();

    // Create the transpose of dg/du.
    RCP<LinearOpBase<Scalar>> A10T =
      createMembers(physics_->get_x_space(), constraintVS_);

    // Create df/dv.
    RCP<LinearOpBase<Scalar>> A01  =
      createMembers(physics_->get_f_space(), constraintVS_);

    // Create dg/dv.
    int numResponses(constraints_.size());
    RCP<DefaultSpmdVectorSpace<Scalar>> responseSpace =
      defaultSpmdVectorSpace<Scalar>(comm_, numResponses, numResponses);
    RCP<LinearOpBase<Scalar>> A11 = createMembers(rcp_dynamic_cast<const
      VectorSpaceBase<Scalar>>(responseSpace), constraintVS_);

    // Create the blocked Jacobian for the system.
    return nonconstBlock2x2<Scalar>(A00, A01, nonconstTranspose(A10T), A11);
  } // end of create_W_op()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  get_W_factory()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  get_W_factory() const
  {
    return physics_->get_W_factory();
  } // end of get_W_factory()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getNominalValues()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
  CurrentConstraintModelEvaluator<Scalar>::
  getNominalValues() const
  {
    using   Teuchos::ArrayRCP;
    using   Teuchos::RCP;
    using   Teuchos::rcp;
    using   Teuchos::rcp_dynamic_cast;
    using   Thyra::assign;
    using   Thyra::DefaultSpmdVector;
    using   Thyra::SpmdVectorSpaceBase;
    using   Thyra::VectorBase;
    typedef Thyra::ModelEvaluatorBase MEB;

    // Get the nominal values from the physics ModelEvaluator.
    MEB::InArgs<Scalar> physNomValues = physics_->getNominalValues();
    RCP<const VectorBase<Scalar>> x, xDot, p, xFull, xFullDot;
    if (physNomValues.supports(MEB::IN_ARG_x))
      x = physNomValues.get_x();
    if (physNomValues.supports(MEB::IN_ARG_x_dot))
      xDot = physNomValues.get_x_dot();
    ArrayRCP<Scalar> initialValues(constraints_.size());
    for (int i(0); i < constraints_.size(); ++i)
      initialValues[i] = rcp_dynamic_cast<const DefaultSpmdVector<Scalar>>
        (physNomValues.get_p(constraints_[i]->parameterIndex()))->
        getRCPtr()[0];
    p = rcp(new DefaultSpmdVector<Scalar>(
      rcp_dynamic_cast<const SpmdVectorSpaceBase<Scalar>>(constraintVS_),
      initialValues, /* stride = */ 1));
    RCP<VectorBase<Scalar>> pDot = p->clone_v();
    assign(pDot.ptr(), 0.0);
    MEB::InArgs<Scalar> nomValues = this->createInArgs();
    nomValues.setArgs(physNomValues);
    if (not x.is_null())
      xFull = buildXVector(x, p);
    if (not xDot.is_null())
      xFullDot = buildXVector(xDot, pDot);
    if (physNomValues.supports(MEB::IN_ARG_x))
      nomValues.set_x(xFull);
    if (physNomValues.supports(MEB::IN_ARG_x_dot))
      nomValues.set_x_dot(xFullDot);
    return nomValues;
  } // end of getNominalValues()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evalModelImpl()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluator<Scalar>::
  evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>&  inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
  {
    using std::string;
    using std::vector;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::add_scalar;
    using Thyra::assign;
    using Thyra::BlockedLinearOpBase;
    using Thyra::createMember;
    using Thyra::DefaultSpmdVector;
    using Thyra::get_ele;
    using Thyra::LinearOpBase;
    using Thyra::MultiVectorBase;
    using Thyra::ProductVectorBase;
    using Thyra::scale;
    using Thyra::ScaledAdjointLinearOpBase;
    using Thyra::set_ele;
    using Thyra::VectorBase;
    using MEB        = Thyra::ModelEvaluatorBase;
    using Derivative = MEB::Derivative<Scalar>;
    MEB::EDerivativeMultiVectorOrientation
      Gradient(MEB::DERIV_MV_GRADIENT_FORM);
    MEB::EDerivativeMultiVectorOrientation
      Jacobian(MEB::DERIV_MV_JACOBIAN_FORM);

    // Get the inArg information.
    bool isTransient = inArgs.supports(MEB::IN_ARG_x_dot);
    RCP<const ProductVectorBase<Scalar>> xFull =
      rcp_dynamic_cast<const ProductVectorBase<Scalar>>(inArgs.get_x(), true);
    RCP<const ProductVectorBase<Scalar>> xDotFull;
    if (isTransient)
      xDotFull = rcp_dynamic_cast<const
        ProductVectorBase<Scalar>>(inArgs.get_x_dot(), true);
    RCP<const VectorBase<Scalar>> x, xDot, parameters;
    if (not xFull.is_null())
    {
      x          = xFull->getVectorBlock(0);
      parameters = xFull->getVectorBlock(1);
    }
    TEUCHOS_ASSERT(not x.is_null())
    TEUCHOS_ASSERT(not parameters.is_null())
    if (isTransient)
    {
      if (not xDotFull.is_null())
        xDot = xDotFull->getVectorBlock(0);
      TEUCHOS_ASSERT(not xDot.is_null())
    }

    // Create the physics inArgs.
    MEB::InArgs<Scalar> physInArgs = physics_->createInArgs();
    physInArgs.setArgs(inArgs);
    for (int i(0); i < constraints_.size(); ++i)
    {
      auto parameter = rcp_dynamic_cast<DefaultSpmdVector<Scalar>>(
        createMember(physics_->get_p_space(i)));
      parameter->getRCPtr()[0] = get_ele(*parameters, i);
      physInArgs.set_p(constraints_[i]->parameterIndex(), parameter);
    } // end loop over constraints_
    physInArgs.set_x(x);
    if (isTransient)
      physInArgs.set_x_dot(xDot);

    // Extract the response information that is required on output.
    vector<RCP<VectorBase<Scalar>>> responseValues(outArgs.Ng());
    for (int ri(0); ri < outArgs.Ng(); ++ri)
      responseValues[ri] = outArgs.get_g(ri);

    // Extract the Jacobian information.
    RCP<BlockedLinearOpBase<Scalar>> A =
      rcp_dynamic_cast<BlockedLinearOpBase<Scalar>>(outArgs.get_W_op(), true);
    RCP<ProductVectorBase<Scalar>>   f =
      rcp_dynamic_cast<ProductVectorBase<Scalar>>(outArgs.get_f(), true);
    RCP<LinearOpBase<Scalar>>    dfdx;
    RCP<MultiVectorBase<Scalar>> dfdp;
    RCP<MultiVectorBase<Scalar>> dgdx;
    RCP<MultiVectorBase<Scalar>> dgdp;
    if (not A.is_null())
    {
      // Extract df/dx.
      dfdx = A->getNonconstBlock(0, 0);
      TEUCHOS_ASSERT(not dfdx.is_null())

      // Extract df/dp.
      RCP<LinearOpBase<Scalar>> A01 = A->getNonconstBlock(0, 1);
      dfdp = rcp_dynamic_cast<MultiVectorBase<Scalar>>(A01, true);
      TEUCHOS_ASSERT(not dfdp.is_null())

      // Extract dg/dx.
      RCP<LinearOpBase<Scalar>> A10 = A->getNonconstBlock(1, 0);
      A10 = rcp_dynamic_cast<ScaledAdjointLinearOpBase<Scalar>>(A10,
        true)->getNonconstOrigOp();
      dgdx = rcp_dynamic_cast<MultiVectorBase<Scalar>>(A10, true);
      TEUCHOS_ASSERT(not dgdx.is_null())

      // Extract dg/dp.
      RCP<LinearOpBase<Scalar>> A11 = A->getNonconstBlock(1, 1);
      dgdp = rcp_dynamic_cast<MultiVectorBase<Scalar>>(A11, true);
      TEUCHOS_ASSERT(not dgdp.is_null())
      assign(dgdp.ptr(), 0.0);

      // Loop over the current constraints...
      for (int i(0); i < constraints_.size(); ++i)
      {
        double value(0);
        if (constraints_[i]->isConstantCurrent())
        {
          // If we have a constant current case, then in 3-D our constraint
          // equation is relatively straightforward:  I = Ia, where I is the
          // current at the contact, Ia is the applied current, and both I and
          // Ia are in Amps.  We can rewrite that as g = I - Ia = 0.  Our
          // parameter is the voltage, that is, p = V, so dg/dp = 0.
          //   In 2-D the current we're working with is actually in Amps per
          // centimeter, so we need to modify our constraint equation using the
          // simulation contact length, L, and the device contact area, A,
          // which have units of micrometers and square micrometers,
          // respectively.  The constraint equation is therefore
          // (I * A / L * 1e-4) = Ia, which we can rewrite as
          // g = I - Ia * (L / A * 1e4) = 0, where the factor of 1e4 comes from
          // the conversion between micrometers and centimeters.  Our parameter
          // is still the voltage, so dg/dp = 0 in 2-D as well.
          value = 0;
        }
        else // if (constraints_[i]->isResistorContact())
        {
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
          value = -1.0 / constraints_[i]->resistorValue();
          if (dimension_ == 2)
            value *= constraints_[i]->contactLength() /
              constraints_[i]->contactArea() * 1e4;
        } // end if this is a Constant Current or Resistor Contact constraint

        // Fill in the diagonal entry of dg/dp.
        set_ele(i, value, dgdp->col(i).ptr());
      } // end loop over the current constraints
    } // end if (not A.is_null())

    // Extract the residual.
    RCP<VectorBase<Scalar>> residual;
    RCP<VectorBase<Scalar>> constraintResidual;
    bool responsesNonNull(true);
    for (int i(0); i < constraints_.size(); ++i)
    {
      if (responseValues[constraints_[i]->responseIndex()].is_null())
      {
        responsesNonNull = false;
        break;
      } // end if the i-th response is null
    } // end loop over constraints_
    if (not f.is_null())
    {
      residual           = f->getNonconstVectorBlock(0);
      constraintResidual = f->getNonconstVectorBlock(1);
      TEUCHOS_ASSERT(not residual.is_null()          )
      TEUCHOS_ASSERT(not constraintResidual.is_null())
    }
    else // if (f.is_null())
    {
      constraintResidual = createMember(*constraintVS_);
      if (responsesNonNull)
        for (int i(0); i < constraints_.size(); ++i)
          set_ele(i,
            get_ele(*responseValues[constraints_[i]->responseIndex()], 0),
            constraintResidual.ptr());
    } // end if (f.is_null()) or not

    // Create the physics outArgs.
    MEB::OutArgs<Scalar> physOutArgs = physics_->createOutArgs();
    physOutArgs.set_W_op(dfdx);
    physOutArgs.set_f(residual);
    if (not dgdx.is_null())
      for (int i(0); i < constraints_.size(); ++i)
        physOutArgs.set_DgDx(constraints_[i]->responseIndex(),
          Derivative(dgdx->col(i), Gradient));
    if (not dfdp.is_null())
      for (int i(0); i < constraints_.size(); ++i)
        physOutArgs.set_DfDp(constraints_[i]->parameterIndex(),
          Derivative(dfdp->col(i), Jacobian));

    // Set up the responses.
    for (int ri(0); ri < outArgs.Ng(); ++ri)
      physOutArgs.set_g(ri, responseValues[ri]);

    // Treat the constraints_[i]->responseIndex() ones specially, as we are
    // using them to form the constraintResidual.
    vector<RCP<VectorBase<Scalar>>> constraintResidualVec;
    for (int i(0); i < constraints_.size(); ++i)
    {
      constraintResidualVec.push_back(createMember(physics_->get_g_space(i)));
      set_ele(0, get_ele(*constraintResidual, i),
        constraintResidualVec[i].ptr());
      physOutArgs.set_g(constraints_[i]->responseIndex(),
        constraintResidualVec[i]);
    } // end loop over constraints_

    // Evaluate the physical model.
    physics_->evalModel(physInArgs, physOutArgs);

    // Reconstruct constraintResidual from the pieces in constraintResidualVec.
    for (int i(0); i < constraints_.size(); ++i)
      set_ele(i, get_ele(*constraintResidualVec[i], 0),
        constraintResidual.ptr());

    // Clean up a corner case in which we simply need to copy the
    // constraintResidual into responseValues.
    if ((f.is_null()) and (responsesNonNull))
    {
      for (int i(0); i < constraints_.size(); ++i)
        set_ele(0, get_ele(*constraintResidual, i),
          responseValues[constraints_[i]->responseIndex()].ptr());
    } // end if ((f.is_null()) and (responsesNonNull))

    // Adjust the constraintResidual.
    if (not f.is_null())
    {
      // Loop over the current constraints...
      for (int i(0); i < constraints_.size(); ++i)
      {
        double value(get_ele(*constraintResidual, i)), adjustment(0);
        if (constraints_[i]->isConstantCurrent())
        {
          // In a constant current case, our constraint equation is
          // g = I - Ia = 0 in 3-D or g = I - Ia * (L / A * 1e4) = 0 in 2-D
          // (see above).  At this point the constraintResidual only contains
          // I, so we need to adjust it by subtracting the applied current.
          adjustment = constraints_[i]->currentValue();
          if (dimension_ == 2)
            adjustment *= constraints_[i]->contactLength() /
              constraints_[i]->contactArea() * 1e4;
        }
        else // if (constraints_[i]->isResistorContact())
        {
          // In a resistor contact case, our constraint equation is
          // g = I - (V - Va) / R = 0 in 3-D or
          // g = I - (V - Va) / R * (L / A * 1e4) = 0 in 2-D (see above).  At
          // this point the constraintResidual only contains I, so we need to
          // adjust it by subtracting the current corresponding to the applied
          // voltage and resistor combination.
          vector<Scalar> voltages;
          getVoltages(voltages, xFull);
          adjustment = (voltages[i] - constraints_[i]->appliedVoltage()) /
            constraints_[i]->resistorValue();
          if (dimension_ == 2)
            adjustment *= constraints_[i]->contactLength() /
              constraints_[i]->contactArea() * 1e4;
        }
        value -= adjustment;
        set_ele(i, value, constraintResidual.ptr());
      } // end loop over the current constraints
    } // end if (not f.is_null())
  } // end of evalModelImpl()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initialize()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluator<Scalar>::
  initialize()
  {
    using Teuchos::Comm;
    using Teuchos::DefaultComm;
    using Teuchos::RCP;
    using Thyra::defaultSpmdVectorSpace;

    // Build a locally-replicated vector space.
    int numConstraints(constraints_.size());
    constraintVS_ =
      defaultSpmdVectorSpace<Scalar>(comm_, numConstraints, numConstraints);
    xSpace_ = buildConstrainedVS(physics_->get_x_space(), constraintVS_);
    fSpace_ = buildConstrainedVS(physics_->get_f_space(), constraintVS_);
  } // end of initialize()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  buildConstrainedVS()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  buildConstrainedVS(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& fullVS,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& constrainedVS)
    const
  {
    using Teuchos::RCP;
    using Thyra::VectorSpaceBase;
    using Thyra::productVectorSpace;
    using std::vector;

    // Augment the full vector space with the constrained one.
    vector<RCP<const VectorSpaceBase<Scalar>>> spaces;
    spaces.push_back(fullVS);
    spaces.push_back(constrainedVS);
    return productVectorSpace<Scalar>(spaces);
  } // end of buildConstrainedVS()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  buildXVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  Teuchos::RCP<const Thyra::VectorBase<Scalar>>
  CurrentConstraintModelEvaluator<Scalar>::
  buildXVector(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& inConstraint) const
  {
    using Teuchos::RCP;
    using Thyra::VectorBase;
    using Thyra::createMember;
    using Thyra::defaultProductVector;
    using std::vector;

    // Augment the solution vector with the constraint.
    vector<RCP<const VectorBase<Scalar>>> vec;
    RCP<const VectorBase<Scalar>> constraint = inConstraint;
    if (constraint.is_null())
      constraint = createMember(constraintVS_);
    vec.push_back(x);
    vec.push_back(constraint);
    return defaultProductVector<Scalar>(xSpace_, vec);
  } // end of buildXVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getVoltages()
  //
  /////////////////////////////////////////////////////////////////////////////
  template <typename Scalar>
  void
  CurrentConstraintModelEvaluator<Scalar>::
  getVoltages(
    std::vector<Scalar>&                                 voltages,
    Teuchos::RCP<const Thyra::ProductVectorBase<Scalar>> xFull) const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::ptrFromRef;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::ProductVectorBase;
    using Thyra::SpmdVectorBase;
    using Thyra::VectorBase;

    // Get the voltage from the solution vector.
    RCP<const VectorBase<double>> voltageBlock =
      rcp_dynamic_cast<const ProductVectorBase<double>>(xFull, true)->
      getVectorBlock(1);
    ArrayRCP<const double> voltageData;
    rcp_dynamic_cast<const SpmdVectorBase<double>>(voltageBlock, true)->
      getLocalData(ptrFromRef(voltageData));
    voltages.clear();
    for (int i(0); i < voltageData.size(); ++i)
      voltages.push_back(voltageData[i]);
  } // end of getVoltages()

} // end of namespace charon

#endif // __Charon_CurrentConstraintModelEvaluator_impl_hpp__
