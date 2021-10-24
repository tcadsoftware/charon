
#ifndef   __Charon_CurrentConstraintList_hpp__
#define   __Charon_CurrentConstraintList_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <ostream>
#include <string>
#include <vector>

// Teuchos
#include "Teuchos_RCP.hpp"

namespace charon
{
  /**
   *  \brief Current Constraint List.
   *
   *  This class is a "list" of all the current constraints we have added to
   *  our device.  Two types of constraints are supported:
   *  - Constant Current:  This constraint consists of:
   *    + the value of the current being applied.
   *    .
   *    There can be at most one Constant Current constraint on the device.
   *  - Resistor Contact:  This constraint consists of:
   *    + the resistance value of the resistor being hooked up; and
   *    + the value of the voltage being applied.
   *  .
   *  All constraints have the following common attributes:
   *  - the sideset ID denoting the terminal to which it's applied;
   *  - the length of the contact in the simulation in \f$ \mu\text{m} \f$;
   *  - the actual device contact area in \f$ \mu\text{m}^2 \f$;
   *  - the initial value given to the voltage parameter corresponding to this
   *    constraint;
   *  - the element block ID corresponding to the sideset ID;
   *  - an index indicating the response to which this constraint corresponds;
   *    and
   *  - an index indicating the "Active Parameter" to which this constraint
   *    corresponds.
   *  .
   *  There can be at most one constraint on any given terminal of the device.
   */
  class
  CurrentConstraintList
  {
    private:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Private Enums, Typedefs, Structs, Classes, etc.              */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Data common to all constraints.
       *
       *  Any constraint we support will have the following associated data:
       *  - `sidesetId_`, which is the sideset ID denoting the terminal to
       *    which we're attaching the constraint;
       *  - `contactLength_`, which is the length of the contact in the
       *    simulation in \f$ \mu\text{m} \f$;
       *  - `contactArea_`, which is the actual area of the device contact in
       *    \f$ \mu\text{m}^2 \f$;
       *  - `initialVoltage_`, which is the initial value (in V) given to the
       *    voltage parameter corresponding to this constraint;
       *  - `elementBlockId_`, which is the element block ID corresponding to
       *    the sideset ID;
       *  - `responseIndex_`, which indicates the response to which this
       *    constraint corresponds; and
       *  - `parameterIndex_`, which indicates the "Active Parameter" to which
       *    this constraint corresponds.
       */
      class ConstraintBase
      {
        public:

          /**
           *  \brief Default Constructor.
           *
           *  Creates a `ConstraintBase` object, storing the three inputs.
           *
           *  \param[in] sidesetId      The sideset ID denoting the terminal to
           *                            which we're attaching the constraint.
           *                            If omitted, this defaults to the empty
           *                            string.
           *  \param[in] contactLength  The length of the contact in the
           *                            simulation (in \f$ \mu\text{m} \f$).
           *                            If omitted, this defaults to 1.
           *  \param[in] contactArea    The actual area of the device contact
           *                            (in \f$ \mu\text{m}^2 \f$).  If
           *                            omitted, this defaults to 1.
           *  \param[in] initialVoltage The initial value (in V) to be given to
           *                            the voltage parameter corresponding to
           *                            this constraint.  If omitted, this
           *                            defaults to 0 V.
           *  \param[in] elementBlockId The element block ID corresponding to
           *                            the sideset ID.  If omitted, this
           *                            defaults to the empty string.
           */
          ConstraintBase(
            const std::string& sidesetId      = "",
            const double&      contactLength  = 1,
            const double&      contactArea    = 1,
            const double&      initialVoltage = 0,
            const std::string& elementBlockId = "");

          /**
           *  \brief Destructor.
           */
          virtual
          ~ConstraintBase()
          {
          } // end of Destructor

          /**
           *  \brief Get the sideset ID.
           */
          std::string
          sidesetId() const;

          /**
           *  \brief Get the contact length in the simulation in
           *         \f$ \mu\text{m} \f$.
           */
          double
          contactLength() const;

          /**
           *  \brief Get the actual device contact area in
           *         \f$ \mu\text{m}^2 \f$.
           */
          double
          contactArea() const;

          /**
           *  \brief Get the initial value of the voltage parameter in Volts.
           */
          double
          initialVoltage() const;

          /**
           *  \brief Get the element block ID.
           */
          std::string
          elementBlockId() const;

          /**
           *  \brief Set the index of the response to which this constraint
           *         corresponds.
           */
          void
          responseIndex(const int& i);

          /**
           *  \brief Get the index of the response to which this constraint
           *         corresponds.
           */
          int
          responseIndex() const;

          /**
           *  \brief Set the index of the "Active Parameter" to which this
           *         constraint corresponds.
           */
          void
          parameterIndex(const int& i);

          /**
           *  \brief Get the index of the "Active Parameter" to which this
           *         constraint corresponds.
           */
          int
          parameterIndex() const;

          /**
           *  \brief Get the type of constraint.
           *
           *  \returns "Undefined" for this base class.
           */
          std::string
          type() const
          {
            return typeImpl();
          } // end of type()

          /**
           *  \brief Is this a Constant Current constraint?
           *
           *  \returns False for this base class.
           */
          bool
          isConstantCurrent() const
          {
            return isConstantCurrentImpl();
          } // end of isConstantCurrent()

          /**
           *  \brief Get the current value (in A).
           *
           *  \throws std::logic_error Because this base class has no current
           *                           value.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          double
          currentValue() const
          {
            return currentValueImpl();
          } // end of currentValue()

          /**
           *  \brief Is this a Resistor Contact constraint?
           *
           *  \returns False for this base class.
           */
          bool
          isResistorContact() const
          {
            return isResistorContactImpl();
          } // end of isResistorContact()

          /**
           *  \brief Get the resistor value (in \f$ \Omega \f$).
           *
           *  \throws std::logic_error Because this base class has no resistor
           *                           value.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          double
          resistorValue() const
          {
            return resistorValueImpl();
          } // end of resistorValue()

          /**
           *  \brief Get the applied voltage (in V).
           *
           *  \throws std::logic_error Because this base class has no applied
           *                           voltage.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          double
          appliedVoltage() const
          {
            return appliedVoltageImpl();
          } // end of appliedVoltage()

          /**
           *  \brief Print out this object.
           *
           *  \param[in,out] os  The output stream to which you'd like to
           *                     print.
           *  \param[in]     tab An optional string with which to prefix each
           *                     line of the ouput.  For instance, you could
           *                     use four spaces within quotes to indent each
           *                     line, or you could use something like "XYZ::"
           *                     for the sake of more easily grepping through
           *                     the output.
           */
          void
          print(
            std::ostream&      os,
            const std::string& tab = "") const
          {
            printImpl(os, tab);
          } // end of print()

          /**
           *  \brief Output streaming operator.
           *
           *  \param[in,out] os The output stream to which you'd like to print.
           *  \param[in]     in The `ConstraintBase` object that you'd like to
           *                    print.
           */
          friend std::ostream& operator<<(
            std::ostream& os,
            const ConstraintBase& in)
          {
            in.printImpl(os);
            return os;
          } // end of operator<<()

        protected:

          /**
           *  \brief Print out this object.
           *
           *  \param[in,out] os  The output stream to which you'd like to
           *                     print.
           *  \param[in]     tab An optional string with which to prefix each
           *                     line of the ouput.  For instance, you could
           *                     use four spaces within quotes to indent each
           *                     line, or you could use something like "XYZ::"
           *                     for the sake of more easily grepping through
           *                     the output.
           */
          virtual void
          printImpl(
            std::ostream&      os,
            const std::string& tab = "") const;

        private:

          /**
           *  \brief Get the type of constraint.
           *
           *  \returns "Undefined" for this base class.
           */
          virtual std::string
          typeImpl() const;

          /**
           *  \brief Is this a Constant Current constraint?
           *
           *  \returns False for this base class.
           */
          virtual bool
          isConstantCurrentImpl() const;

          /**
           *  \brief Get the current value (in A).
           *
           *  \throws std::logic_error Because this base class has no current
           *                           value.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          virtual double
          currentValueImpl() const;

          /**
           *  \brief Is this a Resistor Contact constraint?
           *
           *  \returns False for this base class.
           */
          virtual bool
          isResistorContactImpl() const;

          /**
           *  \brief Get the resistor value (in \f$ \Omega \f$).
           *
           *  \throws std::logic_error Because this base class has no resistor
           *                           value.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          virtual double
          resistorValueImpl() const;

          /**
           *  \brief Get the applied voltage (in V).
           *
           *  \throws std::logic_error Because this base class has no applied
           *                           voltage.
           *
           *  \returns std::numeric_limits<double>::quiet_NaN().
           */
          virtual double
          appliedVoltageImpl() const;

          /**
           *  \brief The sideset ID.
           */
          std::string sidesetId_;

          /**
           *  \brief The contact length in the simulation (in
           *         \f$ \mu\text{m} \f$).
           */
          double contactLength_;

          /**
           *  \brief The actual device contact area (in \f$ \mu\text{m}^2 \f$).
           */
          double contactArea_;

          /**
           *  \brief The initial value (in V) given to the voltage parameter
           *         associated with this constraint.
           */
          double initialVoltage_;

          /**
           *  \brief The element block ID corresponding to the sideset ID.
           */
          std::string elementBlockId_;

          /**
           *  \brief The index of the response to which this constraint
           *         corresponds.
           */
          int responseIndex_;

          /**
           *  \brief The index of the "Active Parameter" to which this
           *         constraint corresponds.
           */
          int parameterIndex_;
      }; // end of class ConstraintBase

      /**
       *  \brief Constant Current Constraint.
       *
       *  A Constant Current constraint consists of:
       *  - `currentValue_`, which is the current value (in A) for the constant
       *    current source.
       */
      class ConstantCurrent
        :
        public ConstraintBase
      {
        public:

          /**
           *  \brief Default Constructor.
           *
           *  Creates a `ConstantCurrent` object, storing the inputs.
           *
           *  \param[in] currentValue   The value (in A) of the applied
           *                            constant current source.
           *  \param[in] sidesetId      The sideset ID denoting the terminal to
           *                            which we're attaching the constraint.
           *  \param[in] contactLength  The length of the contact in the
           *                            simulation (in \f$ \mu\text{m} \f$).
           *                            If omitted, this defaults to 1.
           *  \param[in] contactArea    The actual device contact area (in
           *                            \f$ \mu\text{m}^2 \f$).  If omitted,
           *                            this defaults to 1.
           *  \param[in] initialVoltage The initial value (in V) to be given to
           *                            the voltage parameter corresponding to
           *                            this constraint.  If omitted, this
           *                            defaults to 0 V.
           *  \param[in] elementBlockId The element block ID corresponding to
           *                            the sideset ID.  If omitted, this
           *                            defaults to the empty `string`.
           */
          ConstantCurrent(
            const double&      currentValue,
            const std::string& sidesetId,
            const double&      contactLength  = 1,
            const double&      contactArea    = 1,
            const double&      initialVoltage = 0,
            const std::string& elementBlockID = "");

          /**
           *  \brief Destructor.
           */
          virtual
          ~ConstantCurrent()
          {
          } // end of Destructor

        protected:

          /**
           *  \brief Print out this object.
           *
           *  \param[in,out] os  The output stream to which you'd like to
           *                     print.
           *  \param[in]     tab An optional string with which to prefix each
           *                     line of the ouput.  For instance, you could
           *                     use four spaces within quotes to indent each
           *                     line, or you could use something like "XYZ::"
           *                     for the sake of more easily grepping through
           *                     the output.
           */
          virtual void
          printImpl(
            std::ostream&      os,
            const std::string& tab = "") const;

        private:

          /**
           *  \brief Get the type of constraint.
           *
           *  \returns "ConstantCurrent".
           */
          virtual std::string
          typeImpl() const;

          /**
           *  \brief Is this a Constant Current constraint?
           *
           *  \returns True.
           */
          virtual bool
          isConstantCurrentImpl() const;

          /**
           *  \brief Get the current value (in A).
           *
           *  \returns The amount of current (in A) being applied with this
           *           constant current source.
           */
          virtual double
          currentValueImpl() const;

          /**
           *  \brief How much current (in A) is being applied with the constant
           *         current source.
           */
          double currentValue_;
      }; // end of class ConstantCurrent

      /**
       *  \brief Resistor Contact Constraint.
       *
       *  A Resistor Contact constraint consists of:
       *  - `resistorValue_`, which is the resistance value (in \f$ \Omega \f$)
       *    of the resistor being added; and
       *  - `appliedVoltage_`, which is the value (in V) of the voltage being
       *    applied on the far side of the resistor.
       */
      class ResistorContact
        :
        public ConstraintBase
      {
        public:

          /**
           *  \brief Default Constructor.
           *
           *  Creates a `ResistorContact` object, storing the inputs.
           *
           *  \param[in] resistorValue  The resistance value (in
           *                            \f$ \Omega \f$) of the resistor that's
           *                            hooked up to this contact.
           *  \param[in] appliedVoltage The voltage (in V) being applied on the
           *                            far side of the resistor.
           *  \param[in] sidesetId      The sideset ID denoting the terminal to
           *                            which we're attaching the constraint.
           *  \param[in] contactLength  The length of the contact in the
           *                            simulation (in \f$ \mu\text{m} \f$).
           *                            If omitted, this defaults to 1.
           *  \param[in] contactArea    The actual device contact area (in
           *                            \f$ \mu\text{m}^2 \f$).  If omitted,
           *                            this defaults to 1.
           *  \param[in] initialVoltage The initial value (in V) to be given to
           *                            the voltage parameter corresponding to
           *                            this constraint.  If omitted, this
           *                            defaults to 0.
           *  \param[in] elementBlockId The element block ID corresponding to
           *                            the sideset ID.  If omitted, this
           *                            defaults to the empty `string`.
           */
          ResistorContact(
            const double&      resistorValue,
            const double&      appliedVoltage,
            const std::string& sidesetId,
            const double&      contactLength  = 1,
            const double&      contactArea    = 1,
            const double&      initialVoltage = 0,
            const std::string& elementBlockID = "");

          /**
           *  \brief Destructor.
           */
          virtual
          ~ResistorContact()
          {
          } // end of Destructor

        protected:

          /**
           *  \brief Print out this object.
           *
           *  \param[in,out] os  The output stream to which you'd like to
           *                     print.
           *  \param[in]     tab An optional string with which to prefix each
           *                     line of the ouput.  For instance, you could
           *                     use four spaces within quotes to indent each
           *                     line, or you could use something like "XYZ::"
           *                     for the sake of more easily grepping through
           *                     the output.
           */
          virtual void
          printImpl(
            std::ostream&      os,
            const std::string& tab = "") const;

        private:

          /**
           *  \brief Get the type of constraint.
           *
           *  \returns "ResistorContact".
           */
          virtual std::string
          typeImpl() const;

          /**
           *  \brief Is this a Resistor Contact constraint?
           *
           *  \returns True.
           */
          virtual bool
          isResistorContactImpl() const;

          /**
           *  \brief Get the resistor value (in \f$ \Omega \f$).
           *
           *  \returns The resistance value (in \f$ \Omega \f$)
           *           of the resistor being attached to this contact.
           */
          virtual double
          resistorValueImpl() const;

          /**
           *  \brief Get the applied voltage (in V).
           *
           *  \returns The amount of voltage (in V) being applied on the far
           *           side of the resistor.
           */
          virtual double
          appliedVoltageImpl() const;

          /**
           *  \brief The resistance value (in \f$ \Omega \f$) of the resistor
           *         that you're attaching to the device.
           */
          double resistorValue_;

          /**
           *  \brief The voltage (in V) that you'll be applying on the far side
           *         of the resistor.
           */
          double appliedVoltage_;
      }; // end of class ResistorContact

      /** @} */

    public:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Public Constructors and Destructors.                         */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Default Constructor.
       *
       *  Constructs an empty `CurrentConstraintList`.
       */
      CurrentConstraintList();

      /**
       *  \brief Copy Constructor.
       *
       *  Constructs a copy of the input `CurrentConstraintList`.
       *
       *  \param[in] original The `CurrentConstraintList` that you'd like to
       *                      copy.
       */
      CurrentConstraintList(
        const CurrentConstraintList& original);

      /**
       *  \brief Destructor.
       *
       *  Destroys the `CurrentConstraintList`.
       */
      ~CurrentConstraintList()
      {
      } // end of Destructor

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Public Accessors.                                            */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Do we have a Constant Current constraint?
       *
       *  \returns Whether or not we have a Constant Current constraint.
       */
      bool
      hasConstantCurrent() const;

      /**
       *  \brief Do we have a Resistor Contact constraint?
       *
       *  \returns Whether or not we have a Resistor Contact constraint.
       */
      bool
      hasResistorContact() const;

      /**
       *  \brief Is this list empty?
       *
       *  Returns true if we don't have any constraints.  The intended use case
       *  would be if you have a `CurrentConstraintList` called
       *  `currentConstraints`, you could do something along the lines of:
       *
       *      if (not currentConstraints.empty())
       *      {
       *        // Do something with the constraints...
       *      }
       *
       *  \returns Whether or not this list is empty.
       */
      bool
      empty() const;

      /**
       *  \brief How many Constant Current constraints do we have?
       *
       *  \returns The number of Constant Current constraints, which can be
       *           either 0 or 1.
       */
      int
      numConstantCurrents() const;

      /**
       *  \brief How many Resistor Contact constraints do we have?
       *
       *  \returns The number of Resistor Contact constraints, which can be 0
       *           or greater.
       */
      int
      numResistorContacts() const;

      /**
       *  \brief How many constraints do we have?
       *
       *  \returns The total number of constraints, which can be 0 or greater.
       */
      int
      size() const;

      /**
       *   \brief The beginning of the list.
       *
       *   \returns An `iterator` pointing to the first element in the
       *            `CurrentConstraintList`.
       */
      std::vector<Teuchos::RCP<ConstraintBase>>::const_iterator
      begin() const;

      /**
       *  \brief The end of the list.
       *
       *  \returns An `iterator` pointing to one past the last element in the
       *           `CurrentConstraintList`.
       */
      std::vector<Teuchos::RCP<ConstraintBase>>::const_iterator
      end() const;

      /**
       *  \brief Get the `i`-th constraint in the list.
       *
       *  \throws std::out_of_range If there is no `i`-th constraint in the
       *                            list.
       *
       *  \returns The `i`-th constraint in the list, if it exists; undefined
       *           behavior otherwise.
       */
      Teuchos::RCP<const ConstraintBase>
      operator[](int i) const;

      /**
       *  \brief Get the `i`-th constraint in the list.
       *
       *  \throws std::out_of_range If there is no `i`-th constraint in the
       *                            list.
       *
       *  \returns The `i`-th constraint in the list, if it exists; undefined
       *           behavior otherwise.
       */
      Teuchos::RCP<ConstraintBase>
      operator[](int i);

      /**
       *  \brief Print out this object.
       *
       *  \param[in,out] os  The output stream to which you'd like to print.
       *  \param[in]     tab An optional string with which to prefix each line
       *                     of the ouput.  For instance, you could use four
       *                     spaces within quotes to indent each line, or you
       *                     could use something like "XYZ::" for the sake of
       *                     more easily grepping through the output.
       */
      void
      print(
        std::ostream&      os,
        const std::string& tab = "") const;

      /**
       *  \brief Output streaming operator.
       *
       *  \param[in,out] os The output stream to which you'd like to print.
       *  \param[in]     in The `CurrentConstraintList` object that you'd like
       *                    to print.
       */
      friend std::ostream& operator<<(
        std::ostream& os,
        const CurrentConstraintList& in)
      {
        in.print(os);
        return os;
      } // end of operator<<()

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Public Mutators.                                             */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Clear out the list.
       *
       *  Removes all constraints from this `CurrentConstraintList`, leaving it
       *  as if it had just been constructed.
       */
      void
      clear();

      /**
       *  \brief Add a Constant Current Constraint.
       *
       *  Attempts to add a Constant Current constraint.  There can only be at
       *  most one Constant Current constraint on a device, though, so if a
       *  user attempts to add a second, this operation will fail and return
       *  `false`.  Additionally, we can only attach a single constraint to any
       *  given terminal of the device, so if a user attempts to add this
       *  constraint to a terminal that's already occupied, this operation will
       *  also fail and return `false`.  If there was not yet a Constant
       *  Current constraint on the device, and if the terminal given by
       *  `sidesetId` was empty, then this operation should succeed and return
       *  `true`.
       *
       *  \param[in] currentValue   The current value (in A) for the constant
       *                            current source.
       *  \param[in] sidesetId      The sideset ID denoting the terminal to
       *                            which we're attaching this constraint.
       *  \param[in] contactLength  The length of the contact in the simulation
       *                            (in \f$ \mu\text{m} \f$).  If omitted, the
       *                            default value is 1.
       *  \param[in] contactArea    The actual device contact area (in
       *                            \f$ \mu\text{m}^2 \f$).  If omitted, the
       *                            default value is 1.
       *  \param[in] initialVoltage The initial value (in V) given to the
       *                            voltage parameter that corresponds to this
       *                            constraint.  If omitted, its default value
       *                            is 0.
       *  \param[in] elementBlockId The element block ID for this boundary
       *                            condition.  If omitted, this defaults to an
       *                            empty `string`.
       *
       *  \throws std::logic_error If we attempt to add more than one Constant
       *                           Current constraint to a device, or if we try
       *                           to add more than one constraint to a
       *                           terminal.
       *
       *  \returns Whether or not the operation was successful.
       */
      bool
      addConstantCurrent(
        const double&      currentValue,
        const std::string& sidesetId,
        const double&      contactLength  = 1,
        const double&      contactArea    = 1,
        const double&      initialVoltage = 0,
        const std::string& elementBlockId = "");

      /**
       *  \brief Add a Resistor Contact Constraint.
       *
       *  Attempts to add a Resistor Contact constraint.  There can only be at
       *  most one constraint on each terminal of a device, though, so if a
       *  user attempts to add a second with a `sidesetId` matching a previous
       *  one, this operation will fail and return `false`.  If there was not
       *  yet a constaint on the `sidesetId` terminal, then this operation
       *  should succeed and return `true`.
       *
       *  \param[in] resistorValue  The resistance value (in \f$ \Omega \f$) of
       *                            the resistor being added.
       *  \param[in] appliedVoltage The value (in V) of the voltage being
       *                            applied on the far side of the resistor.
       *  \param[in] sidesetId      The sideset ID denoting the terminal to
       *                            which we're attaching this constraint.
       *  \param[in] contactLength  The length of the contact in the simulation
       *                            (in \f$ \mu\text{m} \f$).  If omitted, this
       *                            defaults to 1.
       *  \param[in] contactArea    The actual device contact area (in
       *                            \f$ \mu\text{m}^2 \f$).  If omitted, this
       *                            defaults to 1.
       *  \param[in] initialVoltage The initial value (in V) given to the
       *                            voltage parameter that corresponds to this
       *                            constraint.  If omitted, its default value
       *                            is 0.
       *  \param[in] elementBlockId The element block ID for this boundary
       *                            condition.  If omitted, this default to an
       *                            empty `string`.
       *
       *  \throws std::logic_error If we attempt to add more than one
       *                           constraint to a terminal.
       *
       *  \returns Whether or not the operation was successful.
       */
      bool
      addResistorContact(
        const double&      resistorValue,
        const double&      appliedVoltage,
        const std::string& sidesetId,
        const double&      contactLength  = 1,
        const double&      contactArea    = 1,
        const double&      initialVoltage = 0,
        const std::string& elementBlockId = "");

      /** @} */

    private:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Private Accessors.                                           */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Is there already a constraint here?
       *
       *  Determine whether or not a constraint is already present on the
       *  terminal denoted by `sidesetId`.  There can only be at most one
       *  constraint on any given terminal of a device.
       *
       *  \returns Whether or not a constraint already exists on `sidesetId`.
       */
      bool
      constraintOnContact(
        const std::string& sidesetId);

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Private Data.                                                */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief A `vector` containing all the current constraints for the
       *         device.
       */
      std::vector<Teuchos::RCP<ConstraintBase>>
      constraints_;

      /**
       *  \brief The number of Constant Current constraints we have.
       */
      int
      numConstantCurrents_;

      /**
       *  \brief The number of Resistor Contact constraints we have.
       */
      int
      numResistorContacts_;

      /** @} */

  }; // end of class CurrentConstraintList

} // end of namespace charon

#endif // __Charon_CurrentConstraintList_hpp__
