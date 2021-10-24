
///////////////////////////////////////////////////////////////////////////////
//
//  tCurrentConstraintList.cpp
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_CurrentConstraintList.hpp"

// Teuchos
#include "Teuchos_UnitTestHarness.hpp"

namespace charon
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  testEmpty()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  testEmpty(
    charon::CurrentConstraintList& in,
    std::ostream&                  out,
    bool&                          success)
  {
    //
    // Test to see whether or not the given CurrentConstraintList is actually
    // empty.
    //
    using std::endl;
    using std::out_of_range;
    using std::string;
    using std::stringstream;

    // Make sure there are no constraints.
    TEST_ASSERT(in.hasConstantCurrent() == false)
    TEST_ASSERT(in.hasResistorContact() == false)
    TEST_ASSERT(in.empty() == true)
    TEST_EQUALITY_CONST(in.numConstantCurrents(), 0)
    TEST_EQUALITY_CONST(in.numResistorContacts(), 0)
    TEST_EQUALITY_CONST(in.size(), 0)
    TEST_ITER_EQUALITY(in.begin(), in.end())
    TEST_THROW(in[0], out_of_range)

    // Test the print() function.
    stringstream os, expected;
    in.print(os, "TEST");
    expected << "TESTCurrentConstraintList:"            << endl
             << "TEST  Summary:"                        << endl
             << "TEST    hasConstantCurrent()  = false" << endl
             << "TEST    hasResistorContact()  = false" << endl
             << "TEST    empty()               = true"  << endl
             << "TEST    numConstantCurrents() = 0"     << endl
             << "TEST    numResistorContacts() = 0"     << endl
             << "TEST    size()                = 0"     << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Test the output streaming operator.
    os.str(string());
    expected.str(string());
    expected << "CurrentConstraintList:"            << endl
             << "  Summary:"                        << endl
             << "    hasConstantCurrent()  = false" << endl
             << "    hasResistorContact()  = false" << endl
             << "    empty()               = true"  << endl
             << "    numConstantCurrents() = 0"     << endl
             << "    numResistorContacts() = 0"     << endl
             << "    size()                = 0"     << endl;
    os << in;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;
  } // end of testEmpty()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  testConstantCurrentExists()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  testConstantCurrentExists(
    charon::CurrentConstraintList& in,
    std::ostream&                  out,
    bool&                          success)
  {
    //
    // Test that the given CurrentConstraintList has a Constant Current
    // constraint on it.
    //
    using std::endl;
    using std::out_of_range;

    // Make sure there are no constraints.
    TEST_ASSERT(in.hasConstantCurrent() == true)
    TEST_ASSERT(in.empty() == false)
    TEST_EQUALITY_CONST(in.numConstantCurrents(), 1)
    TEST_COMPARE_CONST(in.size(), >, 0)
    TEST_ITER_INEQUALITY(in.begin(), in.end())
    TEST_NOTHROW(in[0])
  } // end of testConstantCurrentExists()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  testResistorContactExists()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  testResistorContactExists(
    charon::CurrentConstraintList& in,
    std::ostream&                  out,
    bool&                          success)
  {
    //
    // Test that the given CurrentConstraintList has a Resistor Contact
    // constraint on it.
    //
    using std::endl;
    using std::out_of_range;

    // Make sure there are no constraints.
    TEST_ASSERT(in.hasResistorContact() == true)
    TEST_ASSERT(in.empty() == false)
    TEST_COMPARE_CONST(in.numResistorContacts(), >, 0)
    TEST_COMPARE_CONST(in.size(), >, 0)
    TEST_ITER_INEQUALITY(in.begin(), in.end())
    TEST_NOTHROW(in[0])
  } // end of testResistorContactExists()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList, ConstructEmptyList)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList, ConstructEmptyList)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor
    // and then make sure it really is empty.
    //
    using charon::CurrentConstraintList;
    CurrentConstraintList ccList;
    testEmpty(ccList, out, success);
  } // end of TEUCHOS_UNIT_TEST(CurrentConstraintList, ConstructEmptyList)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList, AddConstantCurrent)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList, AddConstantCurrent)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor,
    // then add a Constant Current constraint to it.  Make sure that constraint
    // was added correctly, then clear the list and make sure it's empty.
    //
    using charon::CurrentConstraintList;
    using std::endl;
    using std::string;
    using std::stringstream;
    CurrentConstraintList ccList;
    ccList.addConstantCurrent(42, "MySideset");
    testConstantCurrentExists(ccList, out, success);

    // Test the print() function.
    stringstream os, expected;
    ccList.print(os, "TEST");
    expected << "TESTCurrentConstraintList:"                           << endl
             << "TEST  Summary:"                                       << endl
             << "TEST    hasConstantCurrent()  = true"                 << endl
             << "TEST    hasResistorContact()  = false"                << endl
             << "TEST    empty()               = false"                << endl
             << "TEST    numConstantCurrents() = 1"                    << endl
             << "TEST    numResistorContacts() = 0"                    << endl
             << "TEST    size()                = 1"                    << endl
             << "TEST  Constraint 1:"                                  << endl
             << "TEST    Type:                       Constant Current" << endl
             << "TEST    Current Value:              42 A"             << endl
             << "TEST    Sideset ID:                 MySideset"        << endl
             << "TEST    Simulation Contact Length:  1 μm"             << endl
             << "TEST    Device Contact Area:        1 μm^2"           << endl
             << "TEST    Initial Voltage:            0 V"              << endl
             << "TEST    Element Block ID:           "                 << endl
             << "TEST    Response Index:             -1"               << endl
             << "TEST    Parameter Index:            -1"               << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Test the output streaming operator.
    os.str(string());
    expected.str(string());
    os << ccList;
    expected << "CurrentConstraintList:"                           << endl
             << "  Summary:"                                       << endl
             << "    hasConstantCurrent()  = true"                 << endl
             << "    hasResistorContact()  = false"                << endl
             << "    empty()               = false"                << endl
             << "    numConstantCurrents() = 1"                    << endl
             << "    numResistorContacts() = 0"                    << endl
             << "    size()                = 1"                    << endl
             << "  Constraint 1:"                                  << endl
             << "    Type:                       Constant Current" << endl
             << "    Current Value:              42 A"             << endl
             << "    Sideset ID:                 MySideset"        << endl
             << "    Simulation Contact Length:  1 μm"             << endl
             << "    Device Contact Area:        1 μm^2"           << endl
             << "    Initial Voltage:            0 V"              << endl
             << "    Element Block ID:           "                 << endl
             << "    Response Index:             -1"               << endl
             << "    Parameter Index:            -1"               << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Clear the list and make sure it's empty.
    ccList.clear();
    testEmpty(ccList, out, success);
  } // end of TEUCHOS_UNIT_TEST(CurrentConstraintList, AddConstantCurrent)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList, AddResistorContact)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList, AddResistorContact)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor,
    // then add a Resistor Contact constraint to it.  Make sure that constraint
    // was added correctly, then clear the list and make sure it's empty.
    //
    using charon::CurrentConstraintList;
    using std::endl;
    using std::string;
    using std::stringstream;
    CurrentConstraintList ccList;
    ccList.addResistorContact(24, 3.14159, "MySideset");
    testResistorContactExists(ccList, out, success);

    // Test the print() function.
    stringstream os, expected;
    ccList.print(os, "TEST");
    expected << "TESTCurrentConstraintList:"                           << endl
             << "TEST  Summary:"                                       << endl
             << "TEST    hasConstantCurrent()  = false"                << endl
             << "TEST    hasResistorContact()  = true"                 << endl
             << "TEST    empty()               = false"                << endl
             << "TEST    numConstantCurrents() = 0"                    << endl
             << "TEST    numResistorContacts() = 1"                    << endl
             << "TEST    size()                = 1"                    << endl
             << "TEST  Constraint 1:"                                  << endl
             << "TEST    Type:                       Resistor Contact" << endl
             << "TEST    Resistor Value:             24 Ω"             << endl
             << "TEST    Applied Voltage:            3.14159 V"        << endl
             << "TEST    Sideset ID:                 MySideset"        << endl
             << "TEST    Simulation Contact Length:  1 μm"             << endl
             << "TEST    Device Contact Area:        1 μm^2"           << endl
             << "TEST    Initial Voltage:            0 V"              << endl
             << "TEST    Element Block ID:           "                 << endl
             << "TEST    Response Index:             -1"               << endl
             << "TEST    Parameter Index:            -1"               << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Test the output streaming operator.
    os.str(string());
    expected.str(string());
    os << ccList;
    expected << "CurrentConstraintList:"                           << endl
             << "  Summary:"                                       << endl
             << "    hasConstantCurrent()  = false"                << endl
             << "    hasResistorContact()  = true"                 << endl
             << "    empty()               = false"                << endl
             << "    numConstantCurrents() = 0"                    << endl
             << "    numResistorContacts() = 1"                    << endl
             << "    size()                = 1"                    << endl
             << "  Constraint 1:"                                  << endl
             << "    Type:                       Resistor Contact" << endl
             << "    Resistor Value:             24 Ω"             << endl
             << "    Applied Voltage:            3.14159 V"        << endl
             << "    Sideset ID:                 MySideset"        << endl
             << "    Simulation Contact Length:  1 μm"             << endl
             << "    Device Contact Area:        1 μm^2"           << endl
             << "    Initial Voltage:            0 V"              << endl
             << "    Element Block ID:           "                 << endl
             << "    Response Index:             -1"               << endl
             << "    Parameter Index:            -1"               << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Clear the list and make sure it's empty.
    ccList.clear();
    testEmpty(ccList, out, success);
  } // end of TEUCHOS_UNIT_TEST(CurrentConstraintList, AddResistorContact)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList,
  //    MultipleConstraintsCopyConstructorOperatorAndIndices)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList,
    MultipleConstraintsCopyConstructorOperatorAndIndices)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor.
    // Add a Constant Current and two Resistor Contact constraints to it.  Set
    // their response and parameter indices, then use the Copy Constructor to
    // duplicate it.  Test that the copy has all the right stuff in it.
    //
    using charon::CurrentConstraintList;
    using std::endl;
    using std::string;
    using std::stringstream;
    CurrentConstraintList ccList;
    ccList.addConstantCurrent(42, "Sideset1", 5, 6, 7, "ElementBlock1");
    ccList.addResistorContact(24, 3.14159, "Sideset2", 8, 9, 10,
      "ElementBlock1");
    ccList.addResistorContact(6, 1.618, "Sideset3", 11, 12, 13,
      "ElementBlock1");
    for (int i(0); i < ccList.size(); ++i)
    {
      ccList[i]->responseIndex(i + 1);
      ccList[i]->parameterIndex(i);
    }
    CurrentConstraintList copyList(ccList);
    testConstantCurrentExists(copyList, out, success);
    testResistorContactExists(copyList, out, success);

    // Test the print() function.
    stringstream os, expected;
    copyList.print(os, "TEST");
    expected << "TESTCurrentConstraintList:"                           << endl
             << "TEST  Summary:"                                       << endl
             << "TEST    hasConstantCurrent()  = true"                 << endl
             << "TEST    hasResistorContact()  = true"                 << endl
             << "TEST    empty()               = false"                << endl
             << "TEST    numConstantCurrents() = 1"                    << endl
             << "TEST    numResistorContacts() = 2"                    << endl
             << "TEST    size()                = 3"                    << endl
             << "TEST  Constraint 1:"                                  << endl
             << "TEST    Type:                       Constant Current" << endl
             << "TEST    Current Value:              42 A"             << endl
             << "TEST    Sideset ID:                 Sideset1"         << endl
             << "TEST    Simulation Contact Length:  5 μm"             << endl
             << "TEST    Device Contact Area:        6 μm^2"           << endl
             << "TEST    Initial Voltage:            7 V"              << endl
             << "TEST    Element Block ID:           ElementBlock1"    << endl
             << "TEST    Response Index:             1"                << endl
             << "TEST    Parameter Index:            0"                << endl
             << "TEST  Constraint 2:"                                  << endl
             << "TEST    Type:                       Resistor Contact" << endl
             << "TEST    Resistor Value:             24 Ω"             << endl
             << "TEST    Applied Voltage:            3.14159 V"        << endl
             << "TEST    Sideset ID:                 Sideset2"         << endl
             << "TEST    Simulation Contact Length:  8 μm"             << endl
             << "TEST    Device Contact Area:        9 μm^2"           << endl
             << "TEST    Initial Voltage:            10 V"             << endl
             << "TEST    Element Block ID:           ElementBlock1"    << endl
             << "TEST    Response Index:             2"                << endl
             << "TEST    Parameter Index:            1"                << endl
             << "TEST  Constraint 3:"                                  << endl
             << "TEST    Type:                       Resistor Contact" << endl
             << "TEST    Resistor Value:             6 Ω"              << endl
             << "TEST    Applied Voltage:            1.618 V"          << endl
             << "TEST    Sideset ID:                 Sideset3"         << endl
             << "TEST    Simulation Contact Length:  11 μm"            << endl
             << "TEST    Device Contact Area:        12 μm^2"          << endl
             << "TEST    Initial Voltage:            13 V"             << endl
             << "TEST    Element Block ID:           ElementBlock1"    << endl
             << "TEST    Response Index:             3"                << endl
             << "TEST    Parameter Index:            2"                << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Test the output streaming operator.
    os.str(string());
    expected.str(string());
    os << copyList;
    expected << "CurrentConstraintList:"                           << endl
             << "  Summary:"                                       << endl
             << "    hasConstantCurrent()  = true"                 << endl
             << "    hasResistorContact()  = true"                 << endl
             << "    empty()               = false"                << endl
             << "    numConstantCurrents() = 1"                    << endl
             << "    numResistorContacts() = 2"                    << endl
             << "    size()                = 3"                    << endl
             << "  Constraint 1:"                                  << endl
             << "    Type:                       Constant Current" << endl
             << "    Current Value:              42 A"             << endl
             << "    Sideset ID:                 Sideset1"         << endl
             << "    Simulation Contact Length:  5 μm"             << endl
             << "    Device Contact Area:        6 μm^2"           << endl
             << "    Initial Voltage:            7 V"              << endl
             << "    Element Block ID:           ElementBlock1"    << endl
             << "    Response Index:             1"                << endl
             << "    Parameter Index:            0"                << endl
             << "  Constraint 2:"                                  << endl
             << "    Type:                       Resistor Contact" << endl
             << "    Resistor Value:             24 Ω"             << endl
             << "    Applied Voltage:            3.14159 V"        << endl
             << "    Sideset ID:                 Sideset2"         << endl
             << "    Simulation Contact Length:  8 μm"             << endl
             << "    Device Contact Area:        9 μm^2"           << endl
             << "    Initial Voltage:            10 V"             << endl
             << "    Element Block ID:           ElementBlock1"    << endl
             << "    Response Index:             2"                << endl
             << "    Parameter Index:            1"                << endl
             << "  Constraint 3:"                                  << endl
             << "    Type:                       Resistor Contact" << endl
             << "    Resistor Value:             6 Ω"              << endl
             << "    Applied Voltage:            1.618 V"          << endl
             << "    Sideset ID:                 Sideset3"         << endl
             << "    Simulation Contact Length:  11 μm"            << endl
             << "    Device Contact Area:        12 μm^2"          << endl
             << "    Initial Voltage:            13 V"             << endl
             << "    Element Block ID:           ElementBlock1"    << endl
             << "    Response Index:             3"                << endl
             << "    Parameter Index:            2"                << endl;
    TEST_ASSERT(os && os.str() == expected.str())
    if (os.str() != expected.str())
      out << "Output Expected:" << endl
          << expected.str()     << endl
          << "Output Received:" << endl
          << os.str()           << endl;

    // Clear the list and make sure it's empty.
    copyList.clear();
    testEmpty(copyList, out, success);
  } // end of TEUCHOS_UNIT_TEST(CurrentConstraintList,
    //   MultipleConstraintsCopyConstructorOperatorAndIndices)

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList, OnlySingleConstantCurrent)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList, OnlySingleConstantCurrent)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor,
    // then add a Constant Current constraint to it.  Try to add a second and
    // make sure that throws an exception.
    //
    using charon::CurrentConstraintList;
    using std::logic_error;
    CurrentConstraintList ccList;
    ccList.addConstantCurrent(42, "Sideset1");
    TEST_THROW(ccList.addConstantCurrent(43, "Sideset2"), logic_error)
  } // end of OnlySingleConstantCurrent

  /////////////////////////////////////////////////////////////////////////////
  //
  //  TEUCHOS_UNIT_TEST(CurrentConstraintList, OneConstraintPerTerminal)
  //
  /////////////////////////////////////////////////////////////////////////////
  TEUCHOS_UNIT_TEST(CurrentConstraintList, OneConstraintPerTerminal)
  {
    //
    // Construct an empty CurrentConstraintList with the Default Constructor,
    // then add a Constant Current constraint to it.  Try to add a Resistor
    // Contact constraint to the same terminal of the device, and ensure that
    // throws an exception.  Also test that adding two Resistor Contact
    // constraints to the same terminal also throws.
    //
    using charon::CurrentConstraintList;
    using std::logic_error;
    CurrentConstraintList ccList;
    ccList.addConstantCurrent(42, "Sideset1");
    TEST_THROW(ccList.addResistorContact(43, 6, "Sideset1"), logic_error)
    ccList.addResistorContact(43, 6, "Sideset2");
    TEST_THROW(ccList.addResistorContact(44, 7, "Sideset2"), logic_error)
  } // end of OneConstraintPerTerminal

} // end of namespace charon

// end of tCurrentConstraintList.cpp
