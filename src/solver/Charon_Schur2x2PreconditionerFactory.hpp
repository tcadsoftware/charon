
///////////////////////////////////////////////////////////////////////////////
//
//  Charon_Schur2x2PreconditionerFactory.hpp
//
///////////////////////////////////////////////////////////////////////////////
#ifndef   __Charon_Schur2x2PreconditionerFactory_hpp__
#define   __Charon_Schur2x2PreconditionerFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Teko
#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_InverseFactory.hpp"

// Teuchos
#include "Teuchos_FancyOStream.hpp"

namespace charon
{
  /**
   *  @brief A preconditioner for a blocked \f$ 2 \times 2 \f$ system using the
   *         Schur complement.
   *
   *  Given a block \f$ 2 \times 2 \f$ matrix
   *  \f[
        A = \begin{bmatrix}
          F & U \\
          L & G
        \end{bmatrix},
      \f]
   *  where \f$ F \in \mathcal{R}^{m \times m}, U
      \in \mathcal{R}^{m \times n}, L \in \mathcal{R}^{n \times m},
      \text{ and }G \in \mathcal{R}^{n \times n} \f$, if \f$ F \f$ is
   *  invertible, the block \f$ LDU \f$ decomposition of \f$ A \f$ is given by
   *  \f[
        A = \begin{bmatrix}
          I       & 0 \\
          LF^{-1} & I
        \end{bmatrix}\begin{bmatrix}
          F & 0 \\
          0 & S
        \end{bmatrix}\begin{bmatrix}
          I & F^{-1}U \\
          0 & I
        \end{bmatrix} = \begin{bmatrix}
          F & U            \\
          L & LF^{-1}U + S
        \end{bmatrix},
      \f]
   *  where \f$ S = G - LF^{-1}U \f$ is known as the Schur complement.  If
   *  \f$ S \f$ is invertible, then
   *  \f[
        A^{-1} = \begin{bmatrix}
          I & -F^{-1}U \\
          0 & I
        \end{bmatrix}\begin{bmatrix}
          F^{-1} & 0      \\
          0      & S^{-1}
        \end{bmatrix}\begin{bmatrix}
          I        & 0 \\
          -LF^{-1} & I
        \end{bmatrix},
      \f]
   *  and an effective left-preconditioner for the system is
   *  \f[
        M^{-1} = \begin{bmatrix}
          F^{-1} & 0      \\
          0      & S^{-1}
        \end{bmatrix}\begin{bmatrix}
          I        & 0 \\
          -LF^{-1} & I
        \end{bmatrix}.
      \f]
   *  This class builds that preconditioner.
   */
  class Schur2x2PreconditionerFactory
    :
    public Teko::BlockPreconditionerFactory
  {
    public:
      /**
       *  @brief Default Constructor.
       *
       *  Does absolutely nothing.
       */
      Schur2x2PreconditionerFactory();

      /**
       *  @brief Build the preconditioner.
       *
       *  Construct the left-preconditioner that comes from the block \f$ LDU
          \f$ decomposition of out input blocked linear operator.
       *
       *  @param[in] op    A blocked linear operator of the form \f$
                           \begin{bmatrix} F & U \\ L & G \end{bmatrix} \f$.
       *  @param[in] state A state object for the preconditioner, from which we
       *                   can `getModifiableOp()`.
       *
       *  @return The preconditioner, given by \f$ \begin{bmatrix} F^{-1} & 0
                  \\ 0 & S^{-1} \end{bmatrix}\begin{bmatrix} I & 0 \\
                  -LF^{-1} & I \end{bmatrix} \f$, where \f$ S = G - LF^{-1}U
                  \f$ is the Schur complement.
       */
      Teko::LinearOp
      buildPreconditionerOperator(
        Teko::BlockedLinearOp&          op,
        Teko::BlockPreconditionerState& state) const;

      /**
       *  @brief Initialize from a `ParameterList`.
       *
       *  Determines the inverse types to use for the inverse of \f$ F \f$ and
       *  for the inverse of \f$ F \f$ used in constructing the Schur
       *  complement, which aren't necessarily the same.  Also builds the
       *  corresponding inverse factories.
       *
       *  @param[in] pl A ParameterList possibly specifying "Inverse Type" and
       *                "Inverse Type for SC".
       */
      virtual void
      initializeFromParameterList(
        const Teuchos::ParameterList& pl);

      /**
       *  @brief Print out an operator.
       *
       *  Prints out the contents of a `Teko::LinearOp` for debugging purposes.
       *
       *  @param[in]     op   The operator you'd like to print out.
       *  @param[in,out] os   The `FancyOStream` to which the operator should
       *                      be printed.
       *  @param[in]     name The name of the operator, if you'd like that
       *                      included in the output.
       */
      virtual void
      printOperator(
        const Teko::LinearOp&  op,
        Teuchos::FancyOStream& os,
        const std::string      name = "op") const;

    protected:

      /**
       *  @brief Inverse Factory.
       *
       *  The factory used to build the inverse of \f$ F \f$.
       */
      Teuchos::RCP<Teko::InverseFactory> invFactory_;

      /**
       *  @brief Inverse Factory for the Schur Complement.
       *
       *  The factory used to build the inverse of \f$ F \f$ used in
       *  constructing the Schur complement.  This is not necessarily the same
       *  as `invFactory_`, though it can be.
       */
      Teuchos::RCP<Teko::InverseFactory> invFactorySC_;

  }; // end of class Schur2x2PreconditionerFactory

} // end of namespace charon

#endif //__Charon_Schur2x2PreconditionerFactory_hpp__

// end of Charon_Schur2x2PreconditionerFactory.hpp
