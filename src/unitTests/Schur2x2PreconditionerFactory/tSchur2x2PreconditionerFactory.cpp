
///////////////////////////////////////////////////////////////////////////////
//
//  tSchur2x2PreconditionerFactory.cpp
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <random>

// Charon
#include "Charon_Schur2x2PreconditionerFactory.hpp"

// Epetra
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

// Teko
#include "Teko_Utilities.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"

// Teuchos
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// Thyra
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSpmdVectorSpace_decl.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_SpmdMultiVectorBase.hpp"

namespace charon
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  createEpetraCrsMatrix()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Epetra_CrsMatrix> createEpetraCrsMatrix(
    std::vector<std::vector<double>> input)
  {
    //
    // This function creates an Epetra_CrsMatrix out of a vector of vectors of
    // doubles.
    //

    // Ensure that the input has the right shape.
    int numRows(input.size()), numCols(input[0].size());
    for (int row(0); row < numRows; ++row)
      TEUCHOS_ASSERT((int)(input[row].size()) == numCols);

    // Create row and column maps for the Epetra_CrsMatrix, and then create the
    // matrix itself.
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    int indexBase(0);
    Teuchos::RCP<Epetra_Map> rowMap =
      Teuchos::rcp(new Epetra_Map(numRows, indexBase, comm));
    Teuchos::RCP<Epetra_Map> colMap =
      Teuchos::rcp(new Epetra_Map(numCols, indexBase, comm));
    Teuchos::RCP<Epetra_CrsMatrix> mat =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, numCols));

    // Determine the number of nonzeros in each row.
    std::vector<int> numNonzeros(numRows, 0);
    int maxNumNonzeros(0);
    for (int row(0); row < numRows; ++row)
    {
      for (int col(0); col < numCols; ++col)
      {
        if (input[row][col] != 0)
          ++numNonzeros[row];
      }
      maxNumNonzeros = std::max(maxNumNonzeros, numNonzeros[row]);
    }

    // Loop over the rows of the input...
    double* values(new double[maxNumNonzeros]);
    int*    indices(new int[maxNumNonzeros]);
    for (int row(0); row < numRows; ++row)
    {
      // If this row is non-empty...
      if (numNonzeros[row] > 0)
      {
        // Loop over the columns of this row...
        for (int col(0), k(0); col < numCols; ++col)
        {
          // ...filling arrays with the nonzero values and their corresponding
          // indices.
          if (input[row][col] != 0)
          {
            values[k]  = input[row][col];
            indices[k] = col;
            ++k;
          }
        }

        // Insert these values into the Epetra_CrsMatrix.
        int result(mat->InsertGlobalValues(row, numNonzeros[row], &values[0],
          &indices[0]));

        // If something goes wrong with the insertion, print an error message.
        if (result != 0)
        {
          Teuchos::FancyOStream ferr(Teuchos::rcpFromRef(std::cerr));
          ferr.setShowProcRank(true);
          ferr.setOutputToRootOnly(-1);
          ferr << "Error in mat->InsertGlobalValues(row, numNonzeros[row], "
               << "&values[0],"                                  << std::endl
               << "&indices[0]) in createEpetraCrsMatrix() with" << std::endl
               << "  row              = " << row                 << std::endl
               << "  numNonzeros[row] = " << numNonzeros[row]    << std::endl
               << "  values           = {";
          for (int i(0); i < numNonzeros[row]; ++i)
          {
            ferr << values[i];
            if (i < numNonzeros[row] - 1)
              ferr << ", ";
          }
          ferr << "}" << std::endl << "  indices          = {";
          for (int i(0); i < numNonzeros[row]; ++i)
          {
            ferr << indices[i];
            if (i < numNonzeros[row] - 1)
              ferr << ", ";
          }
          ferr << "}" << std::endl;
        }
      }
    }
    delete[] values;
    delete[] indices;

    // Specify both the row and column maps in FillComplete() to allow for non-
    // square matrices.
    mat->FillComplete(*colMap, *rowMap);
    return mat;
  } // end of createEpetraCrsMatrix()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  createEpetraMultiVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Epetra_MultiVector> createEpetraMultiVector(
    std::vector<std::vector<double>> input)
  {
    //
    // This function creates an Epetra_MultiVector out of a vector of vectors
    // of doubles.
    //

    // Ensure that the input has the right shape.
    int numRows(input.size()), numCols(input[0].size());
    for (int row(0); row < numRows; ++row)
      TEUCHOS_ASSERT((int)(input[row].size()) == numCols);

    // Create the Epetra_MultiVector.
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    int indexBase(0);
    Teuchos::RCP<Epetra_Map> rowMap =
      Teuchos::rcp(new Epetra_Map(numRows, indexBase, comm));
    Teuchos::RCP<Epetra_MultiVector> vec =
      Teuchos::rcp(new Epetra_MultiVector(*rowMap, numCols));

    // Loop over the rows and columns of the input, placing those values in the
    // corresponding locations in the Epetra_MultiVector.
    for (int row(0); row < numRows; ++row)
    {
      for (int col(0); col < numCols; ++col)
        vec->ReplaceGlobalValue(row, col, input[row][col]);
    }
    return vec;
  } // end of createEpetraMultiVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  createThyraMultiVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::MultiVectorBase<double>>
  createThyraMultiVector(std::vector<std::vector<double>> input)
  {
    //
    // This function creates a Thyra::MultiVectorBase<double> out of a vector
    // of vectors of doubles.
    //
    // It's awfully convoluted--why on earth should I create an
    // Epetra_MultiVector, an Epetra_CrsMatrix, and a Teko::LinearOp before
    // creating a Thyra::MultiVectorBase?--but I know of no other way to create
    // a Thyra::MultiVectorBase without first having an operator from which you
    // can grab a Thyra::VectorSpaceBase.  I'm sure this function could be
    // rewritten in a better way, but that'll have to wait till I understand
    // Thyra better.
    //

    // Ensure that the input has the right shape.
    std::size_t numRows(input.size()), numCols(input[0].size());
    for (std::size_t row(0); row < numRows; ++row)
      TEUCHOS_ASSERT(input[row].size() == numCols);

    // Create an Epetra_MultiVector out of the input, which we'll copy into the
    // Thyra::MultiVectorBase<double> later.
    Teuchos::RCP<Epetra_MultiVector> epetraA = createEpetraMultiVector(input);

    // Create a Teko::LinearOp out of the input, just to get at its range space
    // for the sake of creating the Thyra::MultiVectorBase<double>.
    Teko::LinearOp A = Thyra::epetraLinearOp(createEpetraCrsMatrix(input));
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> vs = A->range();

    // Create the Thyra::MultiVectorBase<double>.
    int numMembers(numCols);
    Teuchos::RCP<Thyra::MultiVectorBase<double>> thyraA =
      Thyra::createMembers(vs, numMembers);

    // Copy the Epetra_MultiVector into the Thyra::MultiVectorBase<double>.
    const Teuchos::RCP<Teko::Epetra::EpetraOperatorWrapper> eow =
      Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(A));
    Teuchos::RCP<const Teko::Epetra::MappingStrategy> ms =
      eow->getMapStrategy();
    ms->copyEpetraIntoThyra(*epetraA, thyraA.ptr());
    return thyraA;
  } // end of createThyraMultiVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  create2x2()
  //
  /////////////////////////////////////////////////////////////////////////////
  const Teko::BlockedLinearOp create2x2(std::vector<std::vector<double>> a,
    std::vector<std::vector<double>> b, std::vector<std::vector<double>> c,
    std::vector<std::vector<double>> d)
  {
    //
    // This function creates a 2x2 blocked linear operator out of some input
    // vectors of vectors of doubles.
    //

    // Create the top-left block as a Thyra::EpetraLinearOp because we'll need
    // to be able to pass it to Teko::buildInverse() later on in
    // Schur2x2PreconditionerFactory::buildPreconditionerOperator().
    Teko::LinearOp A00 = Thyra::epetraLinearOp(createEpetraCrsMatrix(a));

    // Create the rest of the blocks as Thyra::MultiVectorBase<double>s,
    // because that'll need to be their concrete type later on in
    // Schur2x2PreconditionerFactory::buildPreconditionerOperator().
    Teko::LinearOp A01 = createThyraMultiVector(b);
    Teko::LinearOp A10 = createThyraMultiVector(c);
    Teko::LinearOp A11 = createThyraMultiVector(d);

    // Create the Teko::BlockedLinearOp.
    return Teko::toBlockedLinearOp(Thyra::block2x2(A00, A01, A10, A11));
  } // end of create2x2()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  create2x2()
  //
  /////////////////////////////////////////////////////////////////////////////
  const Teko::BlockedLinearOp create2x2(double a, double b, double c, double d)
  {
    //
    // This function creates a 2x2 blocked linear operator out of some doubles.
    //

    // Create vectors of vectors out of the input doubles and pass them on to
    // the method above.
    std::vector<std::vector<double>> aVec = {{a}}, bVec = {{b}}, cVec = {{c}},
      dVec = {{d}};
    return create2x2(aVec, bVec, cVec, dVec);
  } // end of create2x2()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  buildPreconditioner()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teko::LinearOp buildPreconditioner(std::vector<std::vector<double>> a,
    std::vector<std::vector<double>> b, std::vector<std::vector<double>> c,
    std::vector<std::vector<double>> d)
  {
    //
    // This function builds the preconditioner for the blocked linear operator
    // given by the input vectors of vectors of doubles.
    //

    // Create the blocked linear operator from the input.
    Teko::BlockedLinearOp A = create2x2(a, b, c, d);

    // Create an initialize a Schur2x2PreconditionerFactory.
    Teuchos::RCP<Teko::BlockPreconditionerState> state =
      Teuchos::rcp(new Teko::BlockPreconditionerState());
    Schur2x2PreconditionerFactory factory;
    Teuchos::ParameterList pl;
    factory.initializeFromParameterList(pl);

    // Build the preconditioner for the blocked linear operator.
    return factory.buildPreconditionerOperator(A, *state);
  } // end of buildPreconditioner()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  buildPreconditioner()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teko::LinearOp buildPreconditioner(double a, double b, double c, double d)
  {
    //
    // This function builds the preconditioner for the blocked linear operator
    // given by the input doubles.
    //

    // Create vectors of vectors out of the input doubles and pass them on to
    // the method above.
    std::vector<std::vector<double>> aVec = {{a}}, bVec = {{b}}, cVec = {{c}},
      dVec = {{d}};
    return buildPreconditioner(aVec, bVec, cVec, dVec);
  } // end of buildPreconditioner()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  checkPreconditioner()
  //
  /////////////////////////////////////////////////////////////////////////////
  bool checkPreconditioner(Teko::LinearOp invM,
    Teko::BlockedLinearOp invMTruth, Teuchos::FancyOStream& out)
  {
    //
    // This function checks that the preconditioner built by the
    // Schur2x2PreconditionerFactory matches the one built by hand.
    //

    // Set up a Thyra::LinearOpTester.
    Thyra::LinearOpTester<double> tester;
    tester.dump_all(true);
    tester.show_all_tests(true);

    // Compare the two operators.
    const bool result(tester.compare(*invM, *invMTruth,
      Teuchos::ptrFromRef(out)));

    // Print out the results.
    out << "Preconditioner ";
    if (!result)
      out << "not ";
    out << "constructed correctly." << std::endl;
    return result;
  } // end of checkPreconditioner()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, random2x2)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, random2x2)
  {
    //
    // Given a matrix of the form
    //     ⎡a │ b⎤
    // A = ⎢──┼──⎥,
    //     ⎣c │ d⎦
    // the preconditioner should have the form
    //        ⎡   1/a │     0⎤
    // invM = ⎢───────┼──────⎥, where |A| = ad - bc.
    //        ⎣-c/|A| │ a/|A|⎦
    //

    // Set up everything for generating random numbers between 0 and 1.
    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    // Set up the matrix A and build the preconditioner invM.
    double a = dist(e2); double b = dist(e2);
    double c = dist(e2); double d = dist(e2);
    Teko::LinearOp invM = buildPreconditioner(a, b, c, d);

    // Construct the preconditioner by hand.
    double detA = a * d - b * c;
    Teko::BlockedLinearOp invMTruth =
      create2x2(1.0 / a, 0, -c / detA, a / detA);

    // Check to see that invM from the preconditioner routine is the same as
    // the one built by hand.
    TEST_ASSERT(checkPreconditioner(invM, invMTruth, out));
  } // end of TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, random2x2)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, 4x4)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, 4x4)
  {
    //
    // Given the matrix
    //     ⎡ 1  2 │  5   6⎤
    //     ⎢      │       ⎥
    //     ⎢ 3  4 │  7   8⎥
    // A = ⎢──────┼───────⎥,
    //     ⎢ 9 10 │-13 -14⎥
    //     ⎢      │       ⎥
    //     ⎣11 12 │-15 -16⎦
    // the preconditioner should be
    //        ⎡  -2     1 │    0    0⎤
    //        ⎢           │          ⎥
    //        ⎢ 1.5  -0.5 │    0    0⎥
    // invM = ⎢───────────┼──────────⎥.
    //        ⎢  -2   1.5 │    4 -3.5⎥
    //        ⎢           │          ⎥
    //        ⎣1.75 -1.25 │-3.75 3.25⎦
    //

    // Set up the matrix A and build the preconditioner invM.
    std::vector<std::vector<double>> a = {{  1,   2}, {  3,   4}};
    std::vector<std::vector<double>> b = {{  5,   6}, {  7,   8}};
    std::vector<std::vector<double>> c = {{  9,  10}, { 11,  12}};
    std::vector<std::vector<double>> d = {{-13, -14}, {-15, -16}};
    Teko::LinearOp invM = buildPreconditioner(a, b, c, d);

    // Construct the preconditioner by hand.
    a = {{-2,    1}, {  1.5,  -0.5}};
    b = {{ 0,    0}, {    0,     0}};
    c = {{-2,  1.5}, { 1.75, -1.25}};
    d = {{ 4, -3.5}, {-3.75,  3.25}};
    Teko::BlockedLinearOp invMTruth = create2x2(a, b, c, d);

    // Check to see that invM from the preconditioner routine is the same as
    // the one built by hand.
    TEST_ASSERT(checkPreconditioner(invM, invMTruth, out));
  } // end of TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, 4x4)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, columnVector)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, columnVector)
  {
    //
    // Given the matrix
    //     ⎡1 2 │ 5⎤
    //     ⎢    │  ⎥
    // A = ⎢3 4 │ 6⎥,
    //     ⎢────┼──⎥
    //     ⎣7 8 │ 9⎦
    // the preconditioner should be
    //        ⎡ -2    1 │ 0⎤
    //        ⎢         │  ⎥
    // invM = ⎢1.5 -0.5 │ 0⎥.
    //        ⎢─────────┼──⎥
    //        ⎣  2   -3 │ 1⎦
    //

    // Set up the matrix A and build the preconditioner invM.
    std::vector<std::vector<double>> a = {{1, 2}, {3, 4}};
    std::vector<std::vector<double>> b = {{5}, {6}};
    std::vector<std::vector<double>> c = {{7, 8}};
    std::vector<std::vector<double>> d = {{9}};
    Teko::LinearOp invM = buildPreconditioner(a, b, c, d);

    // Construct the preconditioner by hand.
    a = {{-2, 1}, {1.5, -0.5}};
    b = {{0}, {0}};
    c = {{2, -3}};
    d = {{1}};
    Teko::BlockedLinearOp invMTruth = create2x2(a, b, c, d);

    // Check to see that invM from the preconditioner routine is the same as
    // the one built by hand.
    TEST_ASSERT(checkPreconditioner(invM, invMTruth, out));
  } // end of TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, columnVector)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, rectangle)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, rectangle)
  {
    //
    // Given the matrix
    //     ⎡-1  2  3 │ 10 11⎤
    //     ⎢         │      ⎥
    //     ⎢ 4  5  6 │ 12 13⎥
    //     ⎢         │      ⎥
    // A = ⎢ 7  8  9 │ 14 15⎥,
    //     ⎢─────────┼──────⎥
    //     ⎢16 17 18 │ 21 23⎥
    //     ⎢         │      ⎥
    //     ⎣19 20 21 │ 25 27⎦
    // the preconditioner should be
    //        ⎡-1/2    1  -1/2 │   0    0⎤
    //        ⎢                │         ⎥
    //        ⎢   1   -5     3 │   0    0⎥
    //        ⎢                │         ⎥
    // invM = ⎢-1/2 11/3 -13/6 │   0    0⎥.
    //        ⎢────────────────┼─────────⎥
    //        ⎢   0   -2     3 │  -2    1⎥
    //        ⎢                │         ⎥
    //        ⎣   0  5/2  -7/2 │ 3/2 -1/2⎦
    //

    // Set up the matrix A and build the preconditioner invM.
    std::vector<std::vector<double>> a = {{-1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::vector<std::vector<double>> b = {{10, 11}, {12, 13}, {14, 15}};
    std::vector<std::vector<double>> c = {{16, 17, 18}, {19, 20, 21}};
    std::vector<std::vector<double>> d = {{21, 23}, {25, 27}};
    Teko::LinearOp invM = buildPreconditioner(a, b, c, d);

    // Construct the preconditioner by hand.
    a = {{-0.5, 1, -0.5}, {1, -5, 3}, {-0.5, 11.0 / 3.0, -13.0 / 6.0}};
    b = {{0, 0}, {0, 0}, {0, 0}};
    c = {{0, -2, 3}, {0, 2.5, -3.5}};
    d = {{-2, 1}, {1.5, -0.5}};
    Teko::BlockedLinearOp invMTruth = create2x2(a, b, c, d);

    // Check to see that invM from the preconditioner routine is the same as
    // the one built by hand.
    TEST_ASSERT(checkPreconditioner(invM, invMTruth, out));
  } // end of TEUCHOS_UNIT_TEST(Schur2x2PreconditionerFactory, rectangle)

} // end of namespace charon

// end of tSchur2x2PreconditionerFactory.cpp
