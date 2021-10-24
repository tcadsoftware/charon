
#ifndef CHARON_IC_EQUILIBRIUM_DENSITY_DECL_HPP
#define CHARON_IC_EQUILIBRIUM_DENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

namespace charon {

/**
 * @brief Evaluate the equilibrium carrier density:
 * n0 = Nc * exp(-Ec/kbT), and p0 = Nv * exp(Ev/kbT),  Ef = 0 at equilibrium,
 * where Ec and Ev are computed by charon::CondVale_Band.
 *
 * For HBT simulations, Ec and Ev can be computed from Ec = Eref - Chieff - q*Phi,
 * and Ev = Ec - Egeff, when Effective Band Gap and Effective Electron Affinity
 * are given in the input xml.
 */

template<typename EvalT, typename Traits>
class IC_Equilibrium_Density
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IC_Equilibrium_Density(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,BASIS> carrier_density; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos;  // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos;  // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> cond_band;    // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> vale_band;    // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // scaled

  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;    // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> potential;     // scaled

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling [K]

  int num_basis;

  std::string dof_name;

  double effAffinity; // in [eV]
  double effBandGap;  // in [eV]

  bool haveSuffix;

  Teuchos::RCP<const charon::Names>  m_names;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters(bool haveSuffix) const;

}; // end of class IC_Equilibrium_Density


}

#endif
