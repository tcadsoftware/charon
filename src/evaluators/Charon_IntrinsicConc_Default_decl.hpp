
#ifndef CHARON_INTRINSICCONC_DEFAULT_DECL_HPP
#define CHARON_INTRINSICCONC_DEFAULT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! calculate intrinsic concentration
template<typename EvalT, typename Traits>
class IntrinsicConc_Default
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IntrinsicConc_Default(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> ni;
  PHX::MDField<ScalarT,Cell,Point> effEg;
  PHX::MDField<ScalarT,Cell,Point> effChi;

  // input
  PHX::MDField<const ScalarT,Cell,Point> Eg;  // w/o BGN
  PHX::MDField<const ScalarT,Cell,Point> Chi; // w/o BGN

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]

  int num_points;

  // input from closure model factory
  double niValue;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class IntrinsicConc_Default


}

#endif
