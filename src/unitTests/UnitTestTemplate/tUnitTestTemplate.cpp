
///////////////////////////////////////////////////////////////////////////////
//
//  tUnitTestTemplate.cpp
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Include the header file for the class under testing below.  Also include any
// other Trilinos header files that are needed.

// Charon
//#include "Charon_ClassBeingTested.hpp"

// Teuchos
#include "Teuchos_UnitTestHarness.hpp"

// A set of unit tests for a class has the following general form.  Feel free
// to delete this comment.
namespace charon
{
  // If you need to create some functions to be used in the bodies of the unit
  // tests, insert them here.  When done, feel free to delete this comment.

  // This is what a single unit test looks like.  The TEUCHOS_UNIT_TEST macro
  // takes two arguments, the general and specific test names.  The general
  // test name should be the same throughout the file and corresponds to the
  // class under testing.  The specific test name should indicate what exactly
  // is being tested.  Feel free to delete this comment.

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(UnitTestTemplate, specificTestName)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(UnitTestTemplate, specificTestName)
  {
    //
    // Insert some documentation on what this particular unit test is testing.
    // It should be detailed enough that someone else could understand
    // the gist of what's happening.
    //

    // Now insert whatever you need to build the unit test.  See
    // https://cee-gitlab.sandia.gov/tcad-charon/wikis/Creating-a-Unit-Test
    // for details on what Teuchos macros you can use here.
    TEST_ASSERT(true);
  } // end of TEUCHOS_UNIT_TEST(UnitTestTemplate, specificTestName)

  // You can create as many unit tests in this collection as you like.  For
  // any functionality of the class that needs to be tested, create another
  // unit test.  Feel free to delete this comment.

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(UnitTestTemplate, anotherSpecificTestName)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(UnitTestTemplate, anotherSpecificTestName)
  {
    //
    // Insert some documentation on what this additional unit test is testing.
    // Be sure to explain how this test differs from the one above.
    //

    // Again, do whatever you need here to test whatever it is you're testing.
    TEST_EQUALITY(2 + 2, 4);
  } // end of TEUCHOS_UNIT_TEST(UnitTestTemplate, anotherSpecificTestName)

} // end of namespace charon

// end of tUnitTestTemplate.cpp
