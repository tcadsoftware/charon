
///////////////////////////////////////////////////////////////////////////////
//
//  Charon_Schur2x2PreconditionerFactory.cpp
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <vector>

// Charon
#include "Charon_Schur2x2PreconditionerFactory.hpp"

// NOX
#include "NOX_Exceptions.H"

// Teko
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_Utilities.hpp"

// Teuchos
#include "Teuchos_Time.hpp"

// Thyra
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

namespace charon
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Default Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  Schur2x2PreconditionerFactory::
  Schur2x2PreconditionerFactory()
  {
  } // end of Default Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  buildPreconditionerOperator()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teko::LinearOp
  Schur2x2PreconditionerFactory::
  buildPreconditionerOperator(
    Teko::BlockedLinearOp&          op,
    Teko::BlockPreconditionerState& state) const
  {
    using NOX::Exceptions::SolverFailure;
    using std::invalid_argument;
    using std::vector;
    using Teko::blockColCount;
    using Teko::BlockedLinearOp;
    using Teko::blockRowCount;
    using Teko::buildInverse;
    using Teko::createBlockLowerTriInverseOp;
    using Teko::domainSpace;
    using Teko::endBlockFill;
    using Teko::getBlock;
    using Teko::isPhysicallyBlockedLinearOp;
    using Teko::LinearOp;
    using Teko::ModifiableLinearOp;
    using Teko::rangeSpace;
    using Teko::rebuildInverse;
    using Teko::setBlock;
    using Teko::zeroBlockedOp;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::apply;
    using Thyra::createMembers;
    using Thyra::defaultSerialDenseLinearOpWithSolveFactory;
    using Thyra::DefaultZeroLinearOp;
    using Thyra::linearOpWithSolve;
    using Thyra::LinearOpWithSolveBase;
    using Thyra::LinearOpWithSolveFactoryBase;
    using Thyra::MultiVectorBase;
    using Thyra::nonconstInverse;
    using Thyra::scale;
    using Thyra::update;

    // Ensure we're dealing with a 2x2 blocked system, which we'll represent by
    //      ⎡F U⎤
    // op = ⎢   ⎥.  The block preconditioner we'll be building is given by
    //      ⎣L G⎦
    //        ⎡invF 0   ⎤⎡I       0⎤
    // invM = ⎢         ⎥⎢         ⎥.
    //        ⎣0    invS⎦⎣-L*invF I⎦
    TEUCHOS_TEST_FOR_EXCEPTION(not isPhysicallyBlockedLinearOp(op),
      invalid_argument,
      "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "       \
      "Error:  This preconditioner is intended to be used with a blocked "    \
      "2x2 system, but the system given is not blocked.")
    int rows(blockRowCount(op)), cols(blockColCount(op));
    TEUCHOS_TEST_FOR_EXCEPTION((rows != 2) or (cols != 2), invalid_argument,
      "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "       \
      "Error:  This preconditioner is intended to be used with a blocked "    \
      "2x2 system, but the system given is " << rows << "x" << cols << ".")

    // Extract F, and build its inverse.
    const LinearOp F = getBlock(0, 0, op);
    ModifiableLinearOp& invF   = state.getModifiableOp("invF"  );
    ModifiableLinearOp& invFSC = state.getModifiableOp("invFSC");
    if (invF.is_null())
    {
      try
      {
        invF = buildInverse(*invFactory_, F);
      }
      catch (...)
        TEUCHOS_TEST_FOR_EXCEPTION(true, SolverFailure,
          "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "   \
          "I'm afraid it looks like F is not invertible.")
    }
    else // if invF isn't null
    {
      try
      {
        rebuildInverse(*invFactory_, F, invF);
      }
      catch (...)
        TEUCHOS_TEST_FOR_EXCEPTION(true, SolverFailure,
          "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "   \
          "I'm afraid it looks like F is not invertible.")
    } // end if invF is null or not

    // Set up the inverse of F used when building the Schur complement, which
    // may be different from the inverse of F above.
    if (invFactory_ == invFactorySC_)
      invFSC = invF;
    else if (invFSC.is_null())
    {
      try
      {
        invFSC = buildInverse(*invFactorySC_, F);
      }
      catch (...)
        TEUCHOS_TEST_FOR_EXCEPTION(true, SolverFailure,
          "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "   \
          "I'm afraid it looks like F is not invertible.")
    }
    else // if invFSC isn't null
    {
      try
      {
        rebuildInverse(*invFactorySC_, F, invFSC);
      }
      catch (...)
        TEUCHOS_TEST_FOR_EXCEPTION(true, SolverFailure,
          "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  "   \
          "I'm afraid it looks like F is not invertible.")
    } // end if invFSC is null or not

    // Get the upper- and lower-triangular blocks of our operator.
    const LinearOp U = getBlock(0, 1, op);
    const LinearOp L = getBlock(1, 0, op);
    auto mv_U = rcp_dynamic_cast<const MultiVectorBase<double>>(U, true);

    // Form the Schur complement:  S = G - L*invF*U.
    RCP<MultiVectorBase<double>> iFU =
      createMembers(rangeSpace(invFSC), domainSpace(U));
    RCP<MultiVectorBase<double>> S   =
      createMembers(rangeSpace(L), domainSpace(U));
    apply(*invFSC, Thyra::NOTRANS, *mv_U, iFU.ptr());
    apply(*L,      Thyra::NOTRANS, *iFU,  S.ptr()  );
    scale(-1.0, S.ptr());
    const LinearOp G = getBlock(1, 1, op);
    auto zeroG = rcp_dynamic_cast<const DefaultZeroLinearOp<double>>(G);
    if (zeroG.is_null())
    {
      auto mv_G = rcp_dynamic_cast<const MultiVectorBase<double>>(G, true);
      update(1.0, *mv_G, S.ptr());
    } // end if (zeroG.is_null())

    // Create the inverse of S.
    ModifiableLinearOp& invS = state.getModifiableOp("invS");
    try
    {
      RCP<LinearOpWithSolveFactoryBase<double>> dnsInvFact =
        defaultSerialDenseLinearOpWithSolveFactory<double>();
      RCP<LinearOpWithSolveBase<double>> S_withInv =
        linearOpWithSolve<double>(*dnsInvFact, S);
      invS = nonconstInverse<double>(S_withInv);
    } // end of creating the inverse of S
    catch (...)
      TEUCHOS_TEST_FOR_EXCEPTION(true, SolverFailure,
        "Schur2x2PreconditionerFactory::buildPreconditionerOperator():  I'm " \
        "afraid it looks like S is not invertible.")

    // Create and return the preconditioner, invM, described above.
    BlockedLinearOp blockL = zeroBlockedOp(op);
    setBlock(1, 0, blockL, L);
    endBlockFill(blockL);
    vector<LinearOp> invDiag{invF, invS};
    LinearOp invM = createBlockLowerTriInverseOp(blockL, invDiag);
    return invM;
  } // end of buildPreconditionerOperator()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initializeFromParameterList()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  Schur2x2PreconditionerFactory::
  initializeFromParameterList(
    const Teuchos::ParameterList& pl)
  {
    using std::string;
    using Teko::InverseLibrary;
    using Teko::PreconditionerFactory;
    using Teuchos::RCP;

    // Determine the inverse type to use for both inverses of F.
    string invStr("Amesos"), invStrSC("Amesos");
    if (pl.isParameter("Inverse Type"))
      invStr = pl.get<string>("Inverse Type");
    if (pl.isParameter("Inverse Type for SC"))
      invStrSC = pl.get<string>("Inverse Type for SC");

    // Build the inverse factory.
    RCP<const InverseLibrary> invLib =
      PreconditionerFactory::getInverseLibrary();
    invFactory_ = invLib->getInverseFactory(invStr);
    if (invStr != invStrSC)
      invFactorySC_ = invLib->getInverseFactory(invStrSC);
    else
      invFactorySC_ = invFactory_;
  } // end of initializeFromParameterList()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  printOperator()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  Schur2x2PreconditionerFactory::
  printOperator(
    const Teko::LinearOp&  op,
    Teuchos::FancyOStream& os,
    const std::string      name /* = "op" */) const
  {
    using std::endl;
    using Teuchos::Ordinal;
    using Teuchos::RCP;
    using Thyra::apply;
    using Thyra::assign;
    using Thyra::ConstDetachedMultiVectorView;
    using Thyra::createMembers;
    using Thyra::MultiVectorBase;
    using Thyra::set_ele;

    // Get the number of columns and rows in the operator.
    const Ordinal cols = op->domain()->dim();
    const Ordinal rows = op->range()->dim();
    os << name << " has " << cols << " column" << (cols > 1 ? "s" : "")
       << " and " << rows << " row" << (rows > 1 ? "s" : "") << "."
       << endl << "Note:  Only non-zero entries will be printed." << endl;

    // Construct an identity multi-vector.
    const RCP<MultiVectorBase<double>> X = createMembers(op->domain(), cols);
    assign<double>(X.ptr(), 0.0);
    for (Ordinal j(0); j < cols; ++j)
      set_ele(j, 1.0, X->col(j).ptr());

    // Apply the operator to the identity multi-vector to get a multi-vector
    // with the contents of the operator in it.
    const RCP<MultiVectorBase<double>> Y = createMembers(op->range(), cols);
    apply<double>(*op, Thyra::NOTRANS, *X, Y.ptr());

    // Print out the results.
    ConstDetachedMultiVectorView<double> Ydmv(Y);
    for (Ordinal i(0); i < rows; ++i)
      for (Ordinal j(0); j < cols; ++j)
        if (Ydmv(i, j) != 0)
          os << name << "(" << i << ", " << j << ") = " << Ydmv(i, j) << endl;
  } // end of printOperator()

} // end of namespace charon

// end of Charon_Schur2x2PreconditionerFactory.cpp
