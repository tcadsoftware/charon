
#ifndef CHARON_DDLATTICEBC_PERIODIC_DECL_HPP
#define CHARON_DDLATTICEBC_PERIODIC_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_Names.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;

// Evaluate the potential and e/h densities at ohmic contacts for a given
// periodic voltage source

namespace charon {

template<typename EvalT, typename Traits>
class DDLatticeBC_Periodic
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLatticeBC_Periodic(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential;   // scaled, no unit
  PHX::MDField<ScalarT,Cell,BASIS> edensity;
  PHX::MDField<ScalarT,Cell,BASIS> hdensity;

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> doping;      // scaled (net doping)
  PHX::MDField<const ScalarT,Cell,BASIS> acceptor;
  PHX::MDField<const ScalarT,Cell,BASIS> donor;
  PHX::MDField<const ScalarT,Cell,BASIS> gamma_e;
  PHX::MDField<const ScalarT,Cell,BASIS> gamma_h;
  PHX::MDField<const ScalarT,Cell,BASIS> intrin_conc;

  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos; // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos;

  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> ion_density;

  PHX::MDField<const ScalarT,Cell,BASIS> eff_affinity; // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> eff_bandgap;  // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]
  double T0; // [K]
  double t0; // [s]
  double C0; // [cm^-3]

  int num_basis;

  double amplitude;
  double frequency;
  double sign;

  bool bUseFD;
  Teuchos::ParameterList incmpl_ioniz;
  bool bSolveIon;
  bool bUseFermiPin;

  int ion_charge;

  std::string funcType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;

}; // end of class DDLatticeBC_Periodic


}

#endif
