
#ifndef CHARON_SUBTRACT_DECL_HPP
#define CHARON_SUBTRACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;

// Subtract valueB from valueA, i.e., diff = valueA - valueB.
namespace charon {

template<typename EvalT, typename Traits>
class Subtract
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Subtract(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,Point> diff;

  // input
  PHX::MDField<const ScalarT,Cell,Point> valueA;
  PHX::MDField<const ScalarT,Cell,Point> valueB;

  int num_points;

  bool enableA, enableB;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Subtract


}

#endif
