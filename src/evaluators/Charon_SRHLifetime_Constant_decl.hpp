
#ifndef CHARON_SRHLIFETIME_CONSTANT_DECL_HPP
#define CHARON_SRHLIFETIME_CONSTANT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain SRH lifetime for electrons and holes
template<typename EvalT, typename Traits>
class SRHLifetime_Constant
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SRHLifetime_Constant(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> lifetime;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]

  int num_points;

  // input from closure model factory
  double ltValue;  // in [s]

  std::string carrType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SRHLifetime_Constant


}

#endif
