
#ifndef CHARON_FEM_CURRENTDENSITY_DECL_HPP
#define CHARON_FEM_CURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;

// Evaluate the FEM version of current density at IPs

namespace charon {

template<typename EvalT, typename Traits>
class FEM_CurrentDensity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    FEM_CurrentDensity(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> current_density;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> electric_field;
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_carr_dens; // gradient of carrier density
  PHX::MDField<const ScalarT,Cell,IP> carr_dens;
  PHX::MDField<const ScalarT,Cell,IP> diff_coeff;
  PHX::MDField<const ScalarT,Cell,IP> mobility;

  // IP
  std::size_t num_ip;
  std::size_t num_dim;

  // carrier type
  std::string carrType;

  // diffusion coefficient sign
  double sign;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class FEM_CurrentDensity


}

#endif
