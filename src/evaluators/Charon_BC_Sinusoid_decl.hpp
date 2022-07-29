
#ifndef CHARON_BC_SINUSOID_DECL_HPP
#define CHARON_BC_SINUSOID_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

using panzer::Cell;
using panzer::BASIS;

// Sinusoidal voltage at ohmic contact

namespace charon {

template<typename EvalT, typename Traits>
class BC_Sinusoid
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_Sinusoid(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential;   // scaled, no unit
  PHX::MDField<ScalarT,Cell,BASIS> edensity;
  PHX::MDField<ScalarT,Cell,BASIS> hdensity;

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> doping;      // scaled (net doping)
  PHX::MDField<const ScalarT,Cell,BASIS> acceptor;
  PHX::MDField<const ScalarT,Cell,BASIS> donor;
  PHX::MDField<const ScalarT,Cell,BASIS> intrin_conc;

  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos; // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos;

  PHX::MDField<const ScalarT,Cell,BASIS> eff_affinity; // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> eff_bandgap;  // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // scaled

  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;   // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]
  double C0; // [cm^-3]
  double T0; // [K]
  double t0; // time scaling, [s]

  int num_basis;
  double dc_offset;
  double amplitude1;
  double frequency1;
  double phaseshift1;
  double amplitude2;
  double frequency2;
  double phaseshift2;
  bool bUseFD;
  Teuchos::ParameterList incmpl_ioniz;

  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > contactVoltage;
  std::string contactVoltageName;


}; // end of class BC_Sinusoid


}

#endif
