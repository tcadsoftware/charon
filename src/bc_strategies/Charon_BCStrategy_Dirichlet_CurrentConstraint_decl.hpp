
#ifndef   Charon_BCStrategy_Dirichlet_CurrentConstraint_decl_hpp
#define   Charon_BCStrategy_Dirichlet_CurrentConstraint_decl_hpp

// C++
#include <vector>
#include <string>

// Teuchos
#include "Teuchos_RCP.hpp"

// Panzer
#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_Traits.hpp"

// Phalanx
#include "Phalanx_FieldManager.hpp"

namespace charon
{
  /**
   *  \brief Dirichlet Current Constraint BC.
   *
   *  This class sets up the Dirichlet conditions for any Current Constraint
   *  boundary conditions.  The two possible types are either the Constant
   *  Current or Resistor Contact cases.
   */
  template <typename EvalT>
  class BCStrategy_Dirichlet_CurrentConstraint
    :
    public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>
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
       *  Ensures the boundary condition strategy is either "Constant Current"
       *  or "Resistor Contact", and then grabs the `voltageParameter` from the
       *  input `globalData`.
       *
       *  \param[in] bc         The boundary condition information as specified
       *                        in the input XML file.
       *  \param[in] globalData The global data object, which provides access
       *                        to the parameter library, from which we can get
       *                        the `voltageParameter`.
       */
      BCStrategy_Dirichlet_CurrentConstraint(
        const panzer::BC&                       bc,
        const Teuchos::RCP<panzer::GlobalData>& globalData);

      //-----------------------------------------------------------------------
      /** @}                                                                 */
      /** \name Public Mutators.                                             */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Set up the Dirichlet conditions.
       *
       *  Loops over the degrees of freedom (potential and carrier density),
       *  stores their values, creates the corresponding unique residual name,
       *  maps that residual to the degree of freedom and to the target field,
       *  and then stores the degree of freedom basis.
       *
       *  \throws std::runtime_error If `basis_` is still null by the end of
       *                             the routine.
       *
       *  \param[in] sidePB   The physics block from which we can get the
       *                      provided degrees of freedom.
       *  \param[in] userData This actually isn't used in this routine, but
       *                      it's required because this class inherits from
       *                      `panzer::BCStrategy`.
       */
      void
      setup(
        const panzer::PhysicsBlock&   sidePB,
        const Teuchos::ParameterList& userData);

      /**
       *  \brief Build and Register Evaluators.
       *
       *  This routine builds and registers the `BC_CurrentConstraint`
       *  evaluator.  It does so as follows:
       *  - Builds and registers all of the closure models.
       *  - Ensures the element block IDs for the physics block and boundary
       *    condition match.
       *  - Ensures there's only one equation set per physics block.
       *  - Gets the equation set `ParameterList`.
       *  - Determines whether or not Fermi-Dirac is turned on.
       *  - Determines if either Donor or Acceptor Incomplete Ionization is
       *    turned on.
       *  - Creates a `ParameterList` with the appropriate information and
       *    passes it to the `BC_CurrentConstraint` Default Constructor.
       *
       *  \throws std::runtime_error    If the element block IDs from the
       *                                physics block and boundary condition
       *                                don't match or if the physics block has
       *                                more than one equation set.
       *  \throws std::invalid_argument If Donor or Acceptor Incomplete
       *                                Ionization is on, but there is no
       *                                "Model ID", or the "Model ID" isn't a
       *                                list.
       *
       *  \param[in,out] fm       The Field Manager to which the Evaluators
       *                          will be registered.
       *  \param[in,out] pb       The Physics Block.
       *  \param[in,out] factory  The Closure Model Factory, which will allow
       *                          us to build and register the closure models.
       *  \param[in,out] models   The "Closure Models" `ParameterList` from the
       *                          input XML file.
       *  \param[in,out] userData A `ParameterList`, which is used in building
       *                          the registering the closure models, and which
       *                          supplies the scaling parameters.
       */
      void
      buildAndRegisterEvaluators(
        PHX::FieldManager<panzer::Traits>& fm,
        const panzer::PhysicsBlock&        pb,
        const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>&
                                           factory,
        const Teuchos::ParameterList&      models,
        const Teuchos::ParameterList&      userData) const;

      /** @} */

    private:

      //-----------------------------------------------------------------------
      /*                                                                     */
      /** \name Private Data.                                                */
      /** @{                                                                 */
      //-----------------------------------------------------------------------

      /**
       *  \brief Basis.
       *
       *  The basis for the degrees of freedom.  For the time being, we are
       *  assuming all degrees of freedom use the same basis.
       */
      Teuchos::RCP<panzer::PureBasis> basis_;

      /**
       *  \brief Voltage Parameter.
       *
       *  The voltage parameter, whose value is ultimately passed to the Ohmic
       *  Contact boundary condition.
       */
      Teuchos::RCP<panzer::ScalarParameterEntry<EvalT>> voltageParameter_;

      /**
       *  \brief A flag indicating whether or not this current constraint will
       *         also enforce a BJT 1-D Base Contact boundary condition.
       */
      bool bjt1DBaseContact_;

      /** @} */

  }; // end of class BCStrategy_Dirichlet_CurrentConstraint

} // end of namespace charon

#endif // Charon_BCStrategy_Dirichlet_CurrentConstraint_decl_hpp
