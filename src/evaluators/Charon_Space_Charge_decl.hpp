
#ifndef CHARON_SPACE_CHARGE_DECL_HPP
#define CHARON_SPACE_CHARGE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief Compute the scaled space charge rho = p-n+Nd-Na = p-n+doping
 */

template<typename EvalT, typename Traits>
class Space_Charge
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Space_Charge(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> space_charge;

  // input 
  PHX::MDField<ScalarT,Cell,Point> elec_density;
  PHX::MDField<ScalarT,Cell,Point> hole_density;
  PHX::MDField<ScalarT,Cell,Point> doping;
  int num_points;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Space_Charge


}

#endif
