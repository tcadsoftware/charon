
#ifndef CHARON_RECOMBRATE_TOTAL_DECL_HPP
#define CHARON_RECOMBRATE_TOTAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::BASIS;

namespace charon {

//! obtain total recombination rate = sum of all recomb. rates - generation rates
template<typename EvalT, typename Traits>
class RecombRate_Total
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RecombRate_Total(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
//PHX::MDField<ScalarT,Cell,Point> total_rate;
  PHX::MDField<ScalarT,Cell,Point> total_rate;
  PHX::MDField<ScalarT,Cell,Point> total_deriv_e;  // derivative: dRtot/dn
  PHX::MDField<ScalarT,Cell,Point> total_deriv_h;  // derivative: dRtot/dp

  // input
  PHX::MDField<const ScalarT,Cell,Point> srh_rate;
  PHX::MDField<const ScalarT,Cell,Point> trap_srh_rate;
  PHX::MDField<const ScalarT,Cell,Point> rad_rate;
  PHX::MDField<const ScalarT,Cell,Point> auger_rate;
  PHX::MDField<const ScalarT,Cell,Point> ava_rate;
  PHX::MDField<const ScalarT,Cell,Point> defect_cluster_rate;
  PHX::MDField<const ScalarT,Cell,Point> empirical_defect_rate;
  PHX::MDField<const ScalarT,Cell,Point> ionization_particle_strike_rate;
  PHX::MDField<const ScalarT,Cell,Point> opt_gen_rate;

  PHX::MDField<const ScalarT,Cell,Point> srh_deriv_e;
  PHX::MDField<const ScalarT,Cell,Point> srh_deriv_h;
  PHX::MDField<const ScalarT,Cell,Point> trap_srh_deriv_e;
  PHX::MDField<const ScalarT,Cell,Point> trap_srh_deriv_h;
  PHX::MDField<const ScalarT,Cell,Point> rad_deriv_e;
  PHX::MDField<const ScalarT,Cell,Point> rad_deriv_h;
  PHX::MDField<const ScalarT,Cell,Point> auger_deriv_e;
  PHX::MDField<const ScalarT,Cell,Point> auger_deriv_h;
  PHX::MDField<const ScalarT,Cell,Point> ava_deriv_e;
  PHX::MDField<const ScalarT,Cell,Point> ava_deriv_h;

  int num_points;
  int num_nodes;
  std::size_t basis_index;
  std::string basis_name;

  bool bSRH, bTrapSRH, bRadiative, bAuger, bAvalanche, bOptGen;
  bool bDefect, bEmpiricalDefect,bParticleStrike;
  bool isSGCVFEM; 

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class RecombRate_Total


}

#endif
