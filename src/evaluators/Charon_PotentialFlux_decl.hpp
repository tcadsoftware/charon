
#ifndef CHARON_POTENTIALFLUX_DECL_HPP
#define CHARON_POTENTIALFLUX_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;

// Evaluates the scaled potential flux: Lambda2*rel_perm*grad(phi)
namespace charon {

template<typename EvalT, typename Traits>
class PotentialFlux
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PotentialFlux(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> phi_flux;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi;
  PHX::MDField<const ScalarT,Cell,IP> rel_perm;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Lambda2;

  std::size_t num_ip;
  std::size_t num_dim;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class PotentialFlux


}

#endif
