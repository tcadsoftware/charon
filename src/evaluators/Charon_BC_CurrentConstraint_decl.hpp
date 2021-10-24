
#ifndef   Charon_BC_CurrentConstraint_decl_hpp
#define   Charon_BC_CurrentConstraint_decl_hpp

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

// Panzer
#include "Panzer_ScalarParameterEntry.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_MDField.hpp"

namespace charon
{
  /**
   *  \brief Evaluate the potential and electron and hole densities at current
   *         constraint contacts.
   *
   *  Used for either constant current or resistor contact constraints, this
   *  class evaluates the potential and both the electron and hole densities
   *  for the boundary condition.
   */
  template<typename EvalT, typename Traits>
  class BC_CurrentConstraint
    :
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
  {
    public:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Constructors and Destructors.                                */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Default Constructor.
       *
       *  Validates the input `ParameterList`, reads in "Voltage Control",
       *  determines whether or not we're using Fermi-Dirac, gets the
       *  Incomplete Ionization parameters, and then creates and adds the
       *  evalauted and dependent fields.
       *
       *  \throws std::runtime_error If the "Voltage Control" parameter is
       *                             null.
       *
       *  \param[in] p The input `ParameterList`, the contents of which are
       *               specified in `getValidParameters()`.
       */
      BC_CurrentConstraint(
        const Teuchos::ParameterList& p);

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Public Mutators.                                             */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Post-Registration Setup.
       *
       *  For this class, this routine does nothing.
       *
       *  \param[in] d  A `vector` of `panzer::Workset`s.
       *  \param[in] fm The Field Manager.
       */
      void
      postRegistrationSetup(
        typename Traits::SetupData d,
        PHX::FieldManager<Traits>& fm);

      /**
       *  \brief Evaluate Fields.
       *
       *  Evaluate the `OhmicContact` boundary condition.
       *
       *  \param[in] d A `panzer::Workset`, which will be passed to the
       *               `OhmicContact` BC.
       */
      void
      evaluateFields(
        typename Traits::EvalData d);

      /** @} */

    private:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Private Enums, Typedefs, Structs, etc.                       */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief The scalar type.
       */
      typedef typename EvalT::ScalarT ScalarT;

      /**
       *  \brief The type for an evaluated field.
       */
      typedef PHX::MDField<typename EvalT::ScalarT, panzer::Cell,
        panzer::BASIS> EvaluatedField;

      /**
       *  \brief The type for a dependent field.
       */
      typedef PHX::MDField<const typename EvalT::ScalarT, panzer::Cell,
        panzer::BASIS> DependentField;

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Private Accessors.                                           */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Get the valid parameters.
       *
       *  Get the list of all the valid parameters so we can validate the
       *  `ParameterList` passed into the Default Constructor.
       */
      Teuchos::RCP<Teuchos::ParameterList>
      getValidParameters() const;

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Private Data.                                                */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Potential.
       *
       *  An evaluated field, representing the scaled potential, which will be
       *  computed by the `OhmicContact` BC.
       */
      EvaluatedField potential_;

      /**
       *  \brief Electron Density.
       *
       *  An evaluated field, representing the electron density, which will be
       *  computed by the `OhmicContact` BC.
       */
      EvaluatedField eDensity_;

      /**
       *  \brief Hole Density.
       *
       *  An evaluated field, representing the hole density, which will be
       *  computed by the `OhmicContact` BC.
       */
      EvaluatedField hDensity_;

      /**
       *  \brief Doping Concentration.
       *
       *  A dependent field, representing the scaled doping concentration.
       */
      DependentField doping_;

      /**
       *  \brief Acceptor Concentration.
       *
       *  A dependent field, representing the concentration of acceptors.
       */
      DependentField acceptor_;

      /**
       *  \brief Donor Concentration.
       *
       *  A dependent field, representing the concentration of donors.
       */
      DependentField donor_;

      /**
       *  \brief Intrinsic Concentration.
       *
       *  A dependent field, representing the intrinsic carrier concentration
       *  -- that is, the number of electrons in the conduction band (and also
       *  the number of holes in the valence band) -- per unit volume.
       */
      DependentField intrinConc_;

      /**
       *  \brief Electron Effective Density of States.
       *
       *  A dependent field, representing the number of effectively available
       *  electron states per unit volume.
       */
      DependentField eEffDos_;

      /**
       *  \brief Hole Effective Density of States.
       *
       *  A dependent field, representing the number of effectively available
       *  hole states per unit volume.
       */
      DependentField hEffDos_;

      /**
       *  \brief Effective Affinity.
       *
       *  A dependent field, representing the effective electron affinity, or
       *  the energy obtained by moving an electron from the vacuum just
       *  outside the semiconductor to the bottom of the conduction band, in
       *  electron-Volts (eV).
       */
      DependentField effAffinity_;

      /**
       *  \brief Effective Bandgap.
       *
       *  A dependent field, representing the effective band gap, or separation
       *  between the valence and conduction bands, in electron-Volts (eV).
       */
      DependentField effBandgap_;

      /**
       *  \brief Lattice Temperature.
       *
       *  A dependent field, representing the temperature of the semiconductor
       *  lattice.
       */
      DependentField lattTemp_;

      /**
       *  \brief Reference Energy.
       *
       *  A dependent field, representing the energy level in electron-Volts
       *  (eV) from which the potential is referenced.
       */
      DependentField refEnergy_;

      /**
       *  \brief Scaling parameters.
       */
      Teuchos::RCP<charon::Scaling_Parameters> scaleParams_;

      /**
       *  \brief Voltage Scaling [V].
       */
      double V0_;

      /**
       *  \brief Density Scaling [1/cm^3].
       */
      double C0_;

      /**
       *  \brief Temperature Scaling.
       */
      double T0_;

      /**
       *  \brief Voltage.
       */
      Teuchos::RCP<panzer::ScalarParameterEntry<EvalT>> voltageParameter_;

      /**
       *  \brief A flag indicating whether or not to use Fermi-Dirac.
       */
      bool bUseFD_;

      /**
       *  \brief A flag indicating whether or not this current constraint will
       *         enforce a BJT 1-D Base Contact boundary condition.
       */
      bool bjt1DBaseContact_;

      /**
       * \brief A flag indicating whether or not to use Reference Energy
       */
      bool bUseRefE_;

      /**
       *  \brief Incomplete Ionization Parameters.
       */
      Teuchos::ParameterList incmplIoniz_;

      /** @} */

  }; // end of class BC_CurrentConstraint

} // end of namespace charon

#endif // Charon_BC_CurrentConstraint_decl_hpp
