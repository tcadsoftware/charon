
#ifndef CHARON_PDE_RESIDUAL_DD_DECL_HPP
#define CHARON_PDE_RESIDUAL_DD_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include <vector>

using panzer::Cell;
using panzer::IP;
using panzer::Dim;


namespace charon {

template<typename EvalT, typename Traits>
class PDE_Residual_DD
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PDE_Residual_DD(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,IP> residual;

  // input
  PHX::MDField<const ScalarT,Cell,IP, Dim> grad_density;
  PHX::MDField<const ScalarT,Cell,IP, Dim> velocity;
  PHX::MDField<const ScalarT,Cell,IP> total_recomb;
  PHX::MDField<const ScalarT,Cell,IP> dof_time_deriv;

  std::size_t num_ip;
  std::size_t num_dim;

  std::string carrType;

  bool bDofTimeDeriv;

  int ionCharge;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class PDE_Residual_DD


}

#endif
