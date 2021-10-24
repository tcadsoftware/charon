
///////////////////////////////////////////////////////////////////////////////
//
//  Charon_CurrentConstraintModelEvaluator.hpp
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __Charon_CurrentConstraintModelEvaluator_hpp__
#define __Charon_CurrentConstraintModelEvaluator_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_CurrentConstraintList.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Forward Declarations
//
///////////////////////////////////////////////////////////////////////////////
namespace Thyra
{
  template<typename> class DefaultProductVectorSpace;
  template<typename> class ProductVectorBase;
}

namespace charon
{
  /**
   *  \brief Current Constraint `ModelEvaluator`.
   *
   *  This `ModelEvaluator` is used to apply a current constraint to a
   *  boundary.  Two different types of constraints are available:
   *  - Constant Current:\n
   *    In this case we're hooking up a current source to a boundary.  We
   *    therefore need to specify that current value.
   *  - Resistor Contact:\n
   *    In this case we're attaching a resistor to a boundary, and then hooking
   *    up a voltage source to that resistor.  We therefore need to specify the
   *    resistor value and the applied voltage.
   */
  template<typename Scalar>
  class CurrentConstraintModelEvaluator
    :
    public Thyra::ModelEvaluatorDelegatorBase<Scalar>
  {
    public:

      /**
       *  \brief Physics Constructor.
       *
       *  This constructor saves the input as member data and then calls
       *  `initialize()`.
       *
       *  \param[in] physics          The physics `ModelEvaluator`.
       *  \param[in] comm             The MPI communicator.
       *  \param[in] constraints      The list of all the current constraints.
       *  \param[in] dimension        The spatial dimension of the simulation.
       */
      CurrentConstraintModelEvaluator(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
        MPI_Comm                                           comm,
        const charon::CurrentConstraintList&               constraints,
        const int&                                         dimension);

      /**
       *  \brief Get the solution space.
       *
       *  Get the vector space that is the product of the solution and
       *  constraint vector spaces.
       *
       *  \returns The product of the solution and constraint vector spaces.
       */
      Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
      get_x_space() const;

      /**
       *  \brief Get the residual space.
       *
       *  Get the vector space that is the product of the residual and
       *  constraint vector spaces.
       *
       *  \returns The product of the residual and constraint vector spaces.
       */
      Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
      get_f_space() const;

      /**
       *  \brief Create the blocked Jacobian.
       *
       *  Given the nonlinear system of equations
       *  \f[
       *    \begin{align}
       *      f(\vec{u}, \vec{v}) &= 0,\quad\text{and}\\
       *      g(\vec{u}, \vec{v}) &= 0,
       *    \end{align}
       *  \f]
       *  where \f$ \vec{u} \f$ is the solution of some partial differential
       *  equation (PDE), \f$ \vec{v} \f$ is a vector of controls, \f$ f \f$ is
       *  the residual for the nonlinear PDE, and \f$ g \f$ is a set of
       *  nonlinear constraints, we can construct the blocked Jacobian of the
       *  system, given by
       *  \f[
       *    J(u, v) = \begin{bmatrix} \frac{\partial f}{\partial u} &         \
              \frac{\partial f}{\partial v} \\                                \
              \frac{\partial g}{\partial u} &                                 \
              \frac{\partial g}{\partial v} \end{bmatrix}.
       *  \f]
       *  Note that \f$ \frac{\partial f}{\partial u} \f$ is the tangent
       *  stiffness matrix for the unconstrained problem.
       *
       *  \returns The blocked Jacobian of the constrained system.
       */
      Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
      create_W_op() const;

      /**
       *  \brief Get the Jacobian factory.
       *
       *  \returns The Jacobian factory from the physics ModelEvaluator.
       */
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>>
      get_W_factory() const;

      /**
       *  \brief Get the nominal values.
       *
       *  Get the set of nominal values or initial guesses for the supported
       *  input arguments.  In most cases, when a supported input argument is
       *  not specified in an `InArgs` object passed to `evalModel()`, the
       *  nominal value is assumed.
       *
       *  \returns The set of nominal values for the supported input arguments.
       */
      Thyra::ModelEvaluatorBase::InArgs<Scalar>
      getNominalValues() const;

      /**
       *  \brief Get the list of all the constraints.
       *
       *  \returns The `CurrentConstraintList` used to construct this object.
       */
      const charon::CurrentConstraintList&
      constraints() const
      {
        return constraints_;
      } // end of constraints()

    private:

      /**
       *  \brief Default Constructor.
       *
       *  The default constructor is disabled.
       */
      CurrentConstraintModelEvaluator();

      /**
       *  \brief Copy Constructor.
       *
       *  The copy constructor is disabled.
       */
      CurrentConstraintModelEvaluator(
        const CurrentConstraintModelEvaluator& src);

      /**
       *  \brief Evaluate the model.
       *
       *  This is the heart of the class.  Here we create the blocked Jacobian,
       *  augmenting the Jacobian of the original system with the constraint
       *  information.  We then solve the system, and then modify the blocked
       *  residual with the constraint information.
       *
       *  \param[in] inArgs  The input arguments, from which we can get the
       *                     solution vector, and possibly it's
       *                     time-derivative, if it's supported.
       *  \param[in] outArgs The output arguments, from which we can get the
       *                     Jacobian, residual, and constraints.
       */
      void
      evalModelImpl(
        const Thyra::ModelEvaluatorBase::InArgs<Scalar>&  inArgs,
        const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;

      /**
       *  \brief Initialize this `ModelEvaluator`.
       *
       *  Setup internal data structures for vector spaces and operators.
       */
      void
      initialize();

      /**
       *  \brief Build the solution vector space.
       *
       *  This builds a vector space for the constrained system, which is the
       *  product of the vector space for the original system and the one for
       *  the constraints.
       *
       *  \param[in] fullVS   The vector space for the original system.
       *  \param[in] constrVS The vector space for the constraints.
       *
       *  \returns The constrained solution vector space.
       */
      Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar>>
      buildConstrainedVS(
        const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& fullVS,
        const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& constrVS)
        const;

      /**
       *  \brief Build the solution vector.
       *
       *  This creates a constrained solution vector, which is the product of
       *  the original solution vector and the constraint vector.
       *
       *  \param[in] x          The solution vector.
       *  \param[in] constraint The constraint vector.
       *
       *  \returns The constrained solution vector.
       */
      Teuchos::RCP<const Thyra::VectorBase<Scalar>>
      buildXVector(
        const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& x,
        const Teuchos::RCP<const Thyra::VectorBase<Scalar>>& constraint)
        const;

      /**
       *  \brief Get the voltages.
       *
       *  This grabs the voltages out of the constrained solution vector.
       *
       *  \param[in,out] voltages The voltages.
       *  \param[in]     xFull    The constrained solution vector.
       */
      void
      getVoltages(
        std::vector<Scalar>&                                 voltages,
        Teuchos::RCP<const Thyra::ProductVectorBase<Scalar>> xFull) const;

      /**
       *  \brief The physics `ModelEvaluator`.
       */
      Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> physics_;

      /**
       *  \brief The MPI communicator to use with this `ModelEvaluator`.
       */
      Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal>> comm_;

      /**
       *  \brief The list of all the current constraints.
       */
      const charon::CurrentConstraintList constraints_;

      /**
       *  \brief The spatial dimension of the simulation (2D or 3D);
       */
      const int dimension_;

      /**
       *  \brief The solution vector space.
       *
       *  The vector space for the constrained solution, which is the product
       *  of the vector space for the solution of the original system, and the
       *  vector space for the constraints.
       */
      Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar>> xSpace_;

      /**
       *  \brief The residual vector space.
       *
       *  The vector space for the constrained residual, which is the product
       *  of the vector space for the residual of the original system, and the
       *  vector space for the constraints.
       */
      Teuchos::RCP<const Thyra::DefaultProductVectorSpace<Scalar>> fSpace_;

      /**
       *  \brief The vector space for the constraints.
       */
      Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> constraintVS_;

  }; // end of class CurrentConstraintModelEvaluator

} // end of namespace charon

#endif // __Charon_CurrentConstraintModelEvaluator_hpp__

// end of Charon_CurrentConstraintModelEvaluator.hpp
