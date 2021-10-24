
#ifndef CHARON_AVALANCHE_UNIBONEW_DECL_HPP
#define CHARON_AVALANCHE_UNIBONEW_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! obtain UniBoNew avalanche (impact ionization) generation rate
template<typename EvalT, typename Traits>
class Avalanche_UniBoNew
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Avalanche_UniBoNew(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> avalanche_rate;
  PHX::MDField<ScalarT,Cell,Point> ava_deriv_e;
  PHX::MDField<ScalarT,Cell,Point> ava_deriv_h;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp_h;

  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_pot;

  PHX::MDField<const ScalarT,Cell,Point,Dim> curr_dens_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> curr_dens_h;

  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field_h;

  PHX::MDField<const ScalarT,Cell,Point> dens_e;
  PHX::MDField<const ScalarT,Cell,Point> dens_h;

  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_dens_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_dens_h;

  PHX::MDField<const ScalarT,Cell,Point> diff_coeff_e;
  PHX::MDField<const ScalarT,Cell,Point> diff_coeff_h;
  PHX::MDField<const ScalarT,Cell,Point> mobility_e;
  PHX::MDField<const ScalarT,Cell,Point> mobility_h;

  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // lattice temperature

  // Scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double J0;      // current density scaling
  double R0;      // recomb./gen. scaling
  double E0;      // electric field scaling
  double C0;      // concentration scaling
  double T0;      // temperature scaling

  int num_points, num_dims;

  // UniBoNew model parameters
  double a0_e, a1_e, a2_e, a3_e, c0_e, c1_e, c2_e, c3_e, d0_e, d1_e, d2_e, d3_e;
  double b0_e, b1_e, b2_e, b3_e, b4_e, b5_e, b6_e, b7_e, b8_e, b9_e, b10_e;
  double a0_h, a1_h, a2_h, a3_h, c0_h, c1_h, c2_h, c3_h, d0_h, d1_h, d2_h, d3_h;
  double b0_h, b1_h, b2_h, b3_h, b4_h, b5_h, b6_h, b7_h, b8_h, b9_h, b10_h;

  // minimum field below which avalanche generation = 0
  double minField;

  // driving force name
  std::string driveForce;

  // driving force damping parameters [cm^-3]
  double eDrForceRefDens, hDrForceRefDens;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Avalanche_UniBoNew


}

#endif
