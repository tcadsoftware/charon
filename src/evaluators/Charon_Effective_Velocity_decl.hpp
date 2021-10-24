
#ifndef CHARON_EFFECTIVE_VELOCITY_DECL_HPP
#define CHARON_EFFECTIVE_VELOCITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;

namespace charon {

//! obtain effective velocity at IPs
template<typename EvalT, typename Traits>
class Effective_Velocity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Effective_Velocity(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> eff_velocity; // @ IPs

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> velocity; // @ IPs
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_temp;

  PHX::MDField<const ScalarT,Cell,IP> mobility;
  PHX::MDField<const ScalarT,Cell,IP> thermodiff_coeff;

  // IP
  std::size_t num_ip;
  std::size_t num_dim;

  // carrier type
  std::string carrType;

  bool includeSoret;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Effective_Velocity


}

#endif
