
#ifndef CHARON_AVALANCHE_VANOVERSTRAETEN_DECL_HPP
#define CHARON_AVALANCHE_VANOVERSTRAETEN_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! obtain vanOverstraeten avalanche (impact ionization) generation rate
template<typename EvalT, typename Traits>
class Avalanche_vanOverstraeten
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Avalanche_vanOverstraeten(
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

  // vanOverstraeten model parameters
  double al_e, bl_e, ah_e, bh_e, E0_e, hbarOmega_e;
  double al_h, bl_h, ah_h, bh_h, E0_h, hbarOmega_h;

  // minimum field below which avalanche generation = 0
  double minField;

  // driving force name
  std::string driveForce;

  // driving force damping parameters [cm^-3]
  double eDrForceRefDens, hDrForceRefDens;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Avalanche_vanOverstraeten


}

#endif
