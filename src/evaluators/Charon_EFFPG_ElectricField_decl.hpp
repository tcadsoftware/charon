
#ifndef CHARON_EFFPG_ELECTRICFIELD_DECL_HPP
#define CHARON_EFFPG_ELECTRICFIELD_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::BASIS;
using panzer::Dim;

// Compute effective electric field at IPs inversely from the current density
// at IPs that has been calculated using the EFFPG method.

namespace charon {

template<typename EvalT, typename Traits>
class EFFPG_ElectricField
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    EFFPG_ElectricField(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point,Dim> efield;
  PHX::MDField<ScalarT,Cell,Point,Dim> grad_qfp;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> curr_density;  // current density
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_density;  // gradient of carrier density
  PHX::MDField<const ScalarT,Cell,Point> density;           // carrier density

  PHX::MDField<const ScalarT,Cell,Point> diff_coeff;
  PHX::MDField<const ScalarT,Cell,Point> mobility;

  // IP
  std::size_t num_points;
  std::size_t num_dim;

  // carrier type
  std::string carrType;
  double sign;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class EFFPG_ElectricField


}

#endif
