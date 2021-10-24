
#ifndef CHARON_SRHLIFETIME_FUNCTION_DECL_HPP
#define CHARON_SRHLIFETIME_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain SRH lifetime for electrons and holes
template<typename EvalT, typename Traits>
class SRHLifetime_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SRHLifetime_Function(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> lifetime;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> acceptor;
  PHX::MDField<const ScalarT,Cell,Point> donor;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]
  double C0; // concentration scaling, [cm^-3]
  double T0; // temperature scaling, [K]

  int num_points;

  // model parameters
  double tau0, nsrh, Tpowlaw, Texp;

  bool isConcDep, isTempDep, isExpTempDep;


}; // end of class SRHLifetime_Function


}

#endif
