
#ifndef CHARON_INTRINSIC_FERMIENERGY_DECL_HPP
#define CHARON_INTRINSIC_FERMIENERGY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief Calculate the intrinsic Fermi energy in [eV]:
 * Ei = Eref-q*phi-Chieff-0.5*Egeff-0.5*kb*T*log(Nc/Nv)-0.5*kb*T*log(gamma_n/gamma_p).
 *
 * The above expressions hold in general, i.e., valid for Boltzmann and Fermi-Dirac
 * statistics, and for homo- and hetero-junction devices. gamma_n = gamma_p = 1 for
 * Boltzmann statistics, and they are less than 1 for Fermi-Dirac statistics.
 *
 * Reference: D. Schroeder et. al, "Comparison of transport models for the simulation
 * of degenerate semiconductors," Semicond. Sci. Technol. 9, 364 (1994).
 * Note: the formulation in the reference is for isothermal DD simulations only.
 *
 */

template<typename EvalT, typename Traits>
class Intrinsic_FermiEnergy
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Intrinsic_FermiEnergy(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> intrin_fermi;  // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // voltage scaling factor
  double T0; // temperature scaling

  // input
  PHX::MDField<const ScalarT,Cell,Point> ref_energy;    // [eV]
  PHX::MDField<const ScalarT,Cell,Point> eff_affinity;  // [eV]
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;   // [eV]

  PHX::MDField<const ScalarT,Cell,Point> latt_temp;     // scaled by T0
  PHX::MDField<const ScalarT,Cell,Point> potential;     // scaled by V0=kb*T/q

  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;   // scaled by C0
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;

  PHX::MDField<const ScalarT,Cell,Point> elec_degfactor;  // unitless
  PHX::MDField<const ScalarT,Cell,Point> hole_degfactor;

  int num_points;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Intrinsic_FermiEnergy


}

#endif
