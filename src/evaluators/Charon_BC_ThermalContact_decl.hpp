
#ifndef CHARON_BC_THERMALCONTACT_DECL_HPP
#define CHARON_BC_THERMALCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;


namespace charon {

template<typename EvalT, typename Traits>
class BC_ThermalContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_ThermalContact(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,BASIS> latt_temp;   // scaled, no unit

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // Temperature Scaling in [K]

  int num_basis;

  ScalarT user_value;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BC_ThermalContact


}

#endif
