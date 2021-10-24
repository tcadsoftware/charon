
#ifndef CHARON_LATTICETEMP_CONSTANT_DECL_HPP
#define CHARON_LATTICETEMP_CONSTANT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {


template<typename EvalT, typename Traits>
class LatticeTemp_Constant
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    LatticeTemp_Constant(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> latt_temp;  // [scaled]

  // input
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;

  int num_points;

  // input from closure model factory
  double temp;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class LatticeTemp_Constant


}

#endif
