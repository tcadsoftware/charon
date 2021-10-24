
///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <ios>
#include <stdexcept>

// Charon
#include "Charon_CurrentConstraintList.hpp"

// Teuchos
#include "Teuchos_TestForException.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
charon::CurrentConstraintList::ConstraintBase::
ConstraintBase(
  const std::string& sidesetId,      /* = "" */
  const double&      contactLength,  /* = 1  */
  const double&      contactArea,    /* = 1  */
  const double&      initialVoltage, /* = 0  */
  const std::string& elementBlockId  /* = "" */)
  :
  sidesetId_(sidesetId),
  contactLength_(contactLength),
  contactArea_(contactArea),
  initialVoltage_(initialVoltage),
  elementBlockId_(elementBlockId),
  responseIndex_(-1),
  parameterIndex_(-1)
{
} // end of ConstraintBase Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::sidesetId()
//
///////////////////////////////////////////////////////////////////////////////
std::string
charon::CurrentConstraintList::ConstraintBase::
sidesetId() const
{
  return sidesetId_;
} // end of ConstraintBase::sidesetId()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::contactLength()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
contactLength() const
{
  return contactLength_;
} // end of ConstraintBase::contactLength()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::contactArea()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
contactArea() const
{
  return contactArea_;
} // end of ConstraintBase::contactArea()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::initialVoltage()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
initialVoltage() const
{
  return initialVoltage_;
} // end of ConstraintBase::initialVoltage()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::elementBlockId()
//
///////////////////////////////////////////////////////////////////////////////
std::string
charon::CurrentConstraintList::ConstraintBase::
elementBlockId() const
{
  return elementBlockId_;
} // end of ConstraintBase::elementBlockId()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::responseIndex() (Set)
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::ConstraintBase::
responseIndex(const int& i)
{
  responseIndex_ = i;
} // end of ConstraintBase::responseIndex() (Set)

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::responseIndex() (Get)
//
///////////////////////////////////////////////////////////////////////////////
int
charon::CurrentConstraintList::ConstraintBase::
responseIndex() const
{
  return responseIndex_;
} // end of ConstraintBase::responseIndex() (Get)

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::parameterIndex() (Set)
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::ConstraintBase::
parameterIndex(const int& i)
{
  parameterIndex_ = i;
} // end of ConstraintBase::parameterIndex() (Set)

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::parameterIndex() (Get)
//
///////////////////////////////////////////////////////////////////////////////
int
charon::CurrentConstraintList::ConstraintBase::
parameterIndex() const
{
  return parameterIndex_;
} // end of ConstraintBase::parameterIndex() (Get)

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::typeImpl()
//
///////////////////////////////////////////////////////////////////////////////
std::string
charon::CurrentConstraintList::ConstraintBase::
typeImpl() const
{
  return "Undefined";
} // end of ConstraintBase::typeImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::isConstantCurrentImpl()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::ConstraintBase::
isConstantCurrentImpl() const
{
  return false;
} // end of ConstraintBase::isConstantCurrentImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::currentValueImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
currentValueImpl() const
{
  using std::logic_error;
  using std::numeric_limits;
  TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Trying to get the "  \
    "currentValue() out of a non-Constant Current constraint.")
  return numeric_limits<double>::quiet_NaN();
} // end of ConstraintBase::currentValueImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::isResistorContactImpl()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::ConstraintBase::
isResistorContactImpl() const
{
  return false;
} // end of ConstraintBase::isResistorContactImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::resistorValueImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
resistorValueImpl() const
{
  using std::logic_error;
  using std::numeric_limits;
  TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Trying to get the "  \
    "resistorValue() out of a non-Constant Current constraint.")
  return numeric_limits<double>::quiet_NaN();
} // end of ConstraintBase::resistorValueImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::appliedVoltageImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstraintBase::
appliedVoltageImpl() const
{
  using std::logic_error;
  using std::numeric_limits;
  TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Trying to get the "  \
    "appliedVoltage() out of a non-Constant Current constraint.")
  return numeric_limits<double>::quiet_NaN();
} // end of ConstraintBase::appliedVoltageImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstraintBase::printImpl()
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::ConstraintBase::
printImpl(
  std::ostream&      os,
  const std::string& tab /* = "" */) const
{
  using std::endl;
  os << tab << "Sideset ID:                 " << sidesetId_      << endl
     << tab << "Simulation Contact Length:  " << contactLength_  << " μm"
     << endl
     << tab << "Device Contact Area:        " << contactArea_    << " μm^2"
     << endl
     << tab << "Initial Voltage:            " << initialVoltage_ << " V"
     << endl
     << tab << "Element Block ID:           " << elementBlockId_ << endl
     << tab << "Response Index:             " << responseIndex_  << endl
     << tab << "Parameter Index:            " << parameterIndex_ << endl;
} // end of ConstraintBase::printImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstantCurrent Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
charon::CurrentConstraintList::ConstantCurrent::
ConstantCurrent(
  const double&      currentValue,
  const std::string& sidesetId,
  const double&      contactLength,  /* = 1  */
  const double&      contactArea,    /* = 1  */
  const double&      initialVoltage, /* = 0  */
  const std::string& elementBlockId  /* = "" */)
  :
  ConstraintBase(sidesetId, contactLength, contactArea, initialVoltage,
    elementBlockId),
  currentValue_(currentValue)
{
} // end of ConstantCurrent Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  ConstantCurrent::typeImpl()
//
///////////////////////////////////////////////////////////////////////////////
std::string
charon::CurrentConstraintList::ConstantCurrent::
typeImpl() const
{
  return "ConstantCurrent";
} // end of ConstantCurrent::typeImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstantCurrent::isConstantCurrentImpl()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::ConstantCurrent::
isConstantCurrentImpl() const
{
  return true;
} // end of ConstantCurrent::isConstantCurrentImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstantCurrent::currentValueImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ConstantCurrent::
currentValueImpl() const
{
  return currentValue_;
} // end of ConstantCurrent::currentValueImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ConstantCurrent::printImpl()
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::ConstantCurrent::
printImpl(
  std::ostream&      os,
  const std::string& tab /* = "" */) const
{
  using std::endl;
  os << tab << "Type:                       Constant Current"          << endl
     << tab << "Current Value:              " << currentValue_ << " A" << endl;
  ConstraintBase::printImpl(os, tab);
} // end of ConstantCurrent::printImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
charon::CurrentConstraintList::ResistorContact::
ResistorContact(
  const double&      resistorValue,
  const double&      appliedVoltage,
  const std::string& sidesetId,
  const double&      contactLength,  /* = 1  */
  const double&      contactArea,    /* = 1  */
  const double&      initialVoltage, /* = 0  */
  const std::string& elementBlockId  /* = "" */)
  :
  ConstraintBase(sidesetId, contactLength, contactArea, initialVoltage,
    elementBlockId),
  resistorValue_(resistorValue),
  appliedVoltage_(appliedVoltage)
{
} // end of ResistorContact Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact::typeImpl()
//
///////////////////////////////////////////////////////////////////////////////
std::string
charon::CurrentConstraintList::ResistorContact::
typeImpl() const
{
  return "ResistorContact";
} // end of ResistorContact::typeImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact::isResistorContactImpl()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::ResistorContact::
isResistorContactImpl() const
{
  return true;
} // end of ResistorContact::isResistorContactImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact::resistorValueImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ResistorContact::
resistorValueImpl() const
{
  return resistorValue_;
} // end of ResistorContact::resistorValueImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact::appliedVoltageImpl()
//
///////////////////////////////////////////////////////////////////////////////
double
charon::CurrentConstraintList::ResistorContact::
appliedVoltageImpl() const
{
  return appliedVoltage_;
} // end of ResistorContact::appliedVoltageImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  ResistorContact::printImpl()
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::ResistorContact::
printImpl(
  std::ostream&      os,
  const std::string& tab /* = "" */) const
{
  using std::endl;
  os << tab << "Type:                       Resistor Contact"    << endl
     << tab << "Resistor Value:             " << resistorValue_  << " Ω"
     << endl
     << tab << "Applied Voltage:            " << appliedVoltage_ << " V"
     << endl;
  ConstraintBase::printImpl(os, tab);
} // end of ResistorContact::printImpl()

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
charon::CurrentConstraintList::
CurrentConstraintList()
  :
  numConstantCurrents_(0),
  numResistorContacts_(0)
{
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  Copy Constructor
//
///////////////////////////////////////////////////////////////////////////////
charon::CurrentConstraintList::
CurrentConstraintList(
  const CurrentConstraintList& original)
  :
  constraints_(original.constraints_),
  numConstantCurrents_(original.numConstantCurrents_),
  numResistorContacts_(original.numResistorContacts_)
{
} // end of Copy Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  hasConstantCurrent()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
hasConstantCurrent() const
{
  return (numConstantCurrents_ == 1);
} // end of hasConstantCurrent()

///////////////////////////////////////////////////////////////////////////////
//
//  hasResistorContact()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
hasResistorContact() const
{
  return (numResistorContacts_ > 0);
} // end of hasResistorContact()

///////////////////////////////////////////////////////////////////////////////
//
//  empty()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
empty() const
{
  return (constraints_.empty());
} // end of empty()

///////////////////////////////////////////////////////////////////////////////
//
//  numConstantCurrents()
//
///////////////////////////////////////////////////////////////////////////////
int
charon::CurrentConstraintList::
numConstantCurrents() const
{
  return numConstantCurrents_;
} // end of numConstantCurrents()

///////////////////////////////////////////////////////////////////////////////
//
//  numResistorContacts()
//
///////////////////////////////////////////////////////////////////////////////
int
charon::CurrentConstraintList::
numResistorContacts() const
{
  return numResistorContacts_;
} // end of numResistorContacts()

///////////////////////////////////////////////////////////////////////////////
//
//  size()
//
///////////////////////////////////////////////////////////////////////////////
int
charon::CurrentConstraintList::
size() const
{
  return (constraints_.size());
} // end of size()

///////////////////////////////////////////////////////////////////////////////
//
//  begin()
//
///////////////////////////////////////////////////////////////////////////////
std::vector<Teuchos::RCP<charon::CurrentConstraintList::ConstraintBase>>::
  const_iterator
charon::CurrentConstraintList::
begin() const
{
  return constraints_.begin();
} // end of begin()

///////////////////////////////////////////////////////////////////////////////
//
//  end()
//
///////////////////////////////////////////////////////////////////////////////
std::vector<Teuchos::RCP<charon::CurrentConstraintList::ConstraintBase>>::
  const_iterator
charon::CurrentConstraintList::
end() const
{
  return constraints_.end();
} // end of end()

///////////////////////////////////////////////////////////////////////////////
//
//  operator[]() (Const)
//
///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<const charon::CurrentConstraintList::ConstraintBase>
charon::CurrentConstraintList::
operator[](
  int i) const
{
  using std::out_of_range;
  using std::stringstream;
  stringstream message;
  message << "Error:  Attempted to access element " << i << " of the "
          << "CurrentConstrinatList.  ";
  if (constraints_.size() > 0)
    message << "The index must be between 0 and " << constraints_.size() - 1
            << ".";
  else
    message << "The list is empty.";
  TEUCHOS_TEST_FOR_EXCEPTION((i < 0) or (i >= size()), out_of_range,
    message.str())
  return constraints_[i];
} // end of operator[]()

///////////////////////////////////////////////////////////////////////////////
//
//  operator[]() (Non-const)
//
///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<charon::CurrentConstraintList::ConstraintBase>
charon::CurrentConstraintList::
operator[](
  int i)
{
  using std::out_of_range;
  using std::stringstream;
  stringstream message;
  message << "Error:  Attempted to access element " << i << " of the "
          << "CurrentConstrinatList.  ";
  if (constraints_.size() > 0)
    message << "The index must be between 0 and " << constraints_.size() - 1
            << ".";
  else
    message << "The list is empty.";
  TEUCHOS_TEST_FOR_EXCEPTION((i < 0) or (i >= size()), out_of_range,
    message.str())
  return constraints_[i];
} // end of operator[]()

///////////////////////////////////////////////////////////////////////////////
//
//  print()
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::
print(
  std::ostream&      os,
  const std::string& tab /* = "" */) const
{
  using std::endl;
  using std::ios_base;
  ios_base::fmtflags flags = os.setf(ios_base::boolalpha);
  os << tab << "CurrentConstraintList:"                                << endl
     << tab << "  Summary:"                                            << endl
     << tab << "    hasConstantCurrent()  = " << hasConstantCurrent()  << endl
     << tab << "    hasResistorContact()  = " << hasResistorContact()  << endl
     << tab << "    empty()               = " << empty()               << endl
     << tab << "    numConstantCurrents() = " << numConstantCurrents() << endl
     << tab << "    numResistorContacts() = " << numResistorContacts() << endl
     << tab << "    size()                = " << size()                << endl;
  for (int i(0); i < size(); ++i)
  {
    os << tab << "  Constraint " << i + 1 << ":" << endl;
    constraints_[i]->print(os, tab + "    ");
  } // end loop over the constraints
  os.flags(flags);
} // end of print()

///////////////////////////////////////////////////////////////////////////////
//
//  clear()
//
///////////////////////////////////////////////////////////////////////////////
void
charon::CurrentConstraintList::
clear()
{
  constraints_.clear();
  numConstantCurrents_ = 0;
  numResistorContacts_ = 0;
} // end of clear()

///////////////////////////////////////////////////////////////////////////////
//
//  addConstantCurrent()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
addConstantCurrent(
  const double&      currentValue,
  const std::string& sidesetId,
  const double&      contactLength,  /* = 1  */
  const double&      contactArea,    /* = 1  */
  const double&      initialVoltage, /* = 0  */
  const std::string& elementBlockId  /* = "" */)
{
  using std::logic_error;
  using Teuchos::rcp;
  if (hasConstantCurrent())
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Attempting to "    \
      "add a second Constant Current constraint.  Only one Constant Current " \
      "constraint per device is supported.")
    return false;
  }
  if (constraintOnContact(sidesetId))
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Attempting to "    \
      "add a second constraint to the \"" + sidesetId + "\".  Only one "      \
      "constraint per device terminal is supported.")
    return false;
  }
  constraints_.push_back(rcp(new ConstantCurrent(currentValue, sidesetId,
    contactLength, contactArea, initialVoltage, elementBlockId)));
  ++numConstantCurrents_;
  return true;
} // end of addConstantCurrent()

///////////////////////////////////////////////////////////////////////////////
//
//  addResistorContact()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
addResistorContact(
  const double&      resistorValue,
  const double&      appliedVoltage,
  const std::string& sidesetId,
  const double&      contactLength,  /* = 1  */
  const double&      contactArea,    /* = 1  */
  const double&      initialVoltage, /* = 0  */
  const std::string& elementBlockId  /* = "" */)
{
  using std::logic_error;
  using Teuchos::rcp;
  if (constraintOnContact(sidesetId))
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "Error:  Attempting to "    \
      "add a second constraint to the \"" + sidesetId + "\".  Only one "      \
      "constraint per device terminal is supported.")
    return false;
  }
  constraints_.push_back(rcp(new ResistorContact(resistorValue, appliedVoltage,
    sidesetId, contactLength, contactArea, initialVoltage, elementBlockId)));
  ++numResistorContacts_;
  return true;
} // end of addResistorContact()

///////////////////////////////////////////////////////////////////////////////
//
//  constraintOnContact()
//
///////////////////////////////////////////////////////////////////////////////
bool
charon::CurrentConstraintList::
constraintOnContact(
  const std::string& sidesetId)
{
  using std::size_t;
  for (size_t i(0); i < constraints_.size(); ++i)
    if (constraints_[i]->sidesetId() == sidesetId)
      return true;
  return false;
} // end of constraintOnContact()

// end of Charon_CurrentConstraintList.cpp
