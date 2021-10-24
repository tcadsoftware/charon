
#ifndef CHARON_BC_OHMICCONTACT_DECL_HPP
#define CHARON_BC_OHMICCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_Names.hpp"
#include "Charon_FermiDirac_Integral.hpp"

#include "Charon_Scaling_Parameters.hpp"

// Evaluate the potential and e/h densities at ohmic contacts
namespace charon {

  using panzer::Cell;
  using panzer::BASIS;
  using panzer::Point;
  using panzer::Dim;

  // evaluates the e/h densities and potential at Ohmic contacts
  // can be called by all relevant BC evaluators

  template<typename EvalT, typename Traits>
  class OhmicContact {
  public:
    static
    void evaluateOhmicContactBC(const bool& bBJT1DBase, const bool& bUseFD,
         const bool& bUseRefE, const Teuchos::ParameterList& incmpl_ioniz,
         const typename EvalT::ScalarT& voltage,
         const typename EvalT::ScalarT& Eref,
         const typename EvalT::ScalarT& vScaling,
         const typename EvalT::ScalarT& densScaling,
         const typename EvalT::ScalarT& tempScaling,
         const typename Traits::EvalData& workset,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& doping,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& acceptor,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& donor,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& intrin_conc,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& elec_effdos,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& hole_effdos,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& eff_affinity,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& eff_bandgap,
         const PHX::MDField<const typename EvalT::ScalarT,Cell,BASIS>& latt_temp,
         PHX::MDField<typename EvalT::ScalarT,Cell,BASIS>& potential,
         PHX::MDField<typename EvalT::ScalarT,Cell,BASIS>& edensity,
         PHX::MDField<typename EvalT::ScalarT,Cell,BASIS>& hdensity);
  };


template<typename EvalT, typename Traits>
class BC_OhmicContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_OhmicContact(
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
  PHX::MDField<const ScalarT,Cell,BASIS> intrin_conc;

  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos; // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos;

  PHX::MDField<const ScalarT,Cell,BASIS> eff_affinity; // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> eff_bandgap;  // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // scaled

  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;   // [eV]

  //scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]
  double C0; // [cm^-3]
  double T0; // [K]

  int num_basis;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  double small_signal_perturbation; // for freq dom, where a small signal is applied about a steady dc bias 

  bool bUseFD;
  Teuchos::ParameterList incmpl_ioniz;

  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BC_OhmicContact


}

#endif
