
#ifndef CHARON_BC_CONTACTONINSULATOR_DECL_HPP
#define CHARON_BC_CONTACTONINSULATOR_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;

// Evaluate the potential at contacts on insulators
namespace charon {

template<typename EvalT, typename Traits>
class BC_ContactOnInsulator
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_ContactOnInsulator(
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

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential;   // scaled, no unit

  // input
  PHX::MDField<const ScalarT> ref_energy;  // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]
  double t0; // [s]

  std::string basis_name;
  std::size_t basis_index;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > contactVoltage;
  std::string contactVoltageName;
  double small_signal_perturbation; // for freq dom, where a small signal is applied about a steady dc bias

  double work_func = 0.0;
  double initial_time = 0.0, final_time = 0.0;
  double initial_voltage = 0.0, final_voltage = 0.0;   
  double slope = 0.0, yintercept = 0.0; 
  bool bLinRamp = false; 
  
  // for trapezoid pulse
  bool bTrapezoid = false;
  int num_pulses = 1;
  double dc_offset = 0.0;
  double amplitude;
  double period;
  double rise_time;
  double fall_time;
  double delay = 0.0;
  double duty_cycle = 1.0;
  double rising_slope;
  double falling_slope;
  double y_rise_intercept;
  double pulse_time;
  double t1;
  double peak_time;
  double t2;
  double y_fall_intercept;
  double t3;

  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters(bool stochasticParams) const;

}; // end of class BC_ContactOnInsulator


}

#endif
