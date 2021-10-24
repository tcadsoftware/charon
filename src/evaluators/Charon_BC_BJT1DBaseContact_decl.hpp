
#ifndef CHARON_BC_BJT1DBASECONTACT_DECL_HPP
#define CHARON_BC_BJT1DBASECONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;

// Evaluate the built-in potential at ohmic contacts
namespace charon {

template<typename EvalT, typename Traits>
class BC_BJT1DBaseContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_BJT1DBaseContact(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

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

  bool bUseFD;
  Teuchos::ParameterList incmpl_ioniz;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  std::string base_doptype;

private:
  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BC_BJT1DBaseContact


}

#endif
