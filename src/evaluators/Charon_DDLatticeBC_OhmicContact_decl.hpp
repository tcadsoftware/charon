
#ifndef CHARON_DDLATTICEBC_OHMICCONTACT_DECL_HPP
#define CHARON_DDLATTICEBC_OHMICCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_Names.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;


namespace charon {


// select a physical solution for carrier density at contact
double select_physical_sol(std::vector<double>& sol, double dop);

// compute carrier density with MB Statistics at an Ohmic contact
// when incomplete ionization model is active
double compute_MB_carrier_dens (
     int q, const Teuchos::ParameterList& incmpl_ioniz_param,
     double kbT, double Nc, double Nv, double ni, double Na,
     double Nd, double Nion,
     double dens_sc);

// compute carrier density with FD Statistics at an Ohmic contact
// when incomplete ionization model is active
double compute_FD_carrier_dens (
                int q, const Teuchos::ParameterList& incmpl_ioniz_param,
                double kbT, double Nc, double Nv, double ni, double Na,
                double Nd, double Nion, double gamma_n, double gamma_p,
                double dens_sc);


// Evaluate the potential and e/h densities at ohmic contacts
template<typename EvalT, typename Traits>
class DDLatticeBC_OhmicContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLatticeBC_OhmicContact(
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

  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;   // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]
  double T0; // [K]
  double C0; // [cm^-3]

  int num_basis;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > contactVoltage;
  std::string contactVoltageName;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;

  bool bUseFD;
  Teuchos::ParameterList incmpl_ioniz;
  bool bSolveIon;
  bool bUseFermiPin;

  int ion_charge;
  double ionDens;
  double initial_voltage = 0.0;  

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;

}; // end of class DDLatticeBC_OhmicContact


}

#endif
