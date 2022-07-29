
#ifndef CHARON_PREVPOTENTIALGRAD_DECL_HPP
#define CHARON_PREVPOTENTIALGRAD_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"


using panzer::Cell;
using panzer::IP;
using panzer::Dim;

// Evaluate the grad potential at IPs for a previous time

namespace charon {

template<typename EvalT, typename Traits>
class PrevPotentialGrad
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  using ScalarT = typename EvalT::ScalarT;

public:
  PrevPotentialGrad(
    const Teuchos::ParameterList& p);

  void
  evaluateFields(
    typename Traits::EvalData d);

private:
  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> prev_grad_phi;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi; 
 
  // IP
  std::size_t num_ip;
  std::size_t num_dim;

  double time;

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class PrevPotentialGrad


}

#endif
