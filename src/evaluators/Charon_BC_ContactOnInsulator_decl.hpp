
#ifndef CHARON_BC_CONTACTONINSULATOR_DECL_HPP
#define CHARON_BC_CONTACTONINSULATOR_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;

// Evaluate the potential at contacts on insulators
namespace charon {

template<typename EvalT, typename Traits>
class BC_ContactOnInsulator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_ContactOnInsulator(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential;   // scaled, no unit

  // input
  PHX::MDField<const ScalarT> ref_energy;  // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]

  std::string basis_name;
  std::size_t basis_index;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  double work_func;

private:
  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters(bool stochasticParams) const;

}; // end of class BC_ContactOnInsulator


}

#endif
