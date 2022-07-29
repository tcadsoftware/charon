
#ifndef CHARON_DISPLACEMENTCURRENTDENSITY_DECL_HPP
#define CHARON_DISPLACEMENTCURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"


using panzer::Cell;
using panzer::IP;
using panzer::Dim;

// Evaluate the SGCVFEM version of displacement current density at IPs

namespace charon {

template<typename EvalT, typename Traits>
class DisplacementCurrentDensity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
public:

  DisplacementCurrentDensity(
    const Teuchos::ParameterList& p);

  void
  evaluateFields(
    typename Traits::EvalData d);

private:

  using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> current_density;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi_prev;
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi; 
  PHX::MDField<const ScalarT,Cell,IP> rel_perm;
 
  // IP
  std::size_t num_ip;
  std::size_t num_dim;

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0;
  double E0;
  double J0;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DisplacementCurrentDensity


}

#endif
