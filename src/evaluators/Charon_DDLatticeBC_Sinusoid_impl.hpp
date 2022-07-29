
#ifndef CHARON_DDLATTICEBC_SINUSOID_IMPL_HPP
#define CHARON_DDLATTICEBC_SINUSOID_IMPL_HPP

#include <cmath>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_DDLatticeBC_OhmicContact_decl.hpp"
#include "Charon_Util.hpp"


/**
 * @brief The Ohmic BC calculations here for electrons and holes are applicable to
 * the DDIon, and DDIonLattice formulations. The potential used to be 
 * the intrinsic Fermi potential, i.e., -q*\phi = Ei. However it was changed to 
 * the shifted vacumm potential in January 2021 such that the solution of the 
 * Poisson equation is the shifted vacumm potential for all equation sets (not just 
 * isothermal DD equation sets). The reason for the change is that the potential 
 * in Charon's framework is always continuous even across a heterojunction, while 
 * the intrinsic Fermi potential is physically discontinuous in a heterojunction. 
 * Therefore, it is not correct to set the Poisson equation solution as the intrinsic
 * Fermi potential in Charon. In some commercial codes such as Silvaco, the potential 
 * is set to be the intrinsic Fermi potential, because they use double points at a 
 * heterojunction so that the potential is allowed to be discontinuous. 
 * The voltage source is a sinusoidal function in time, with input params given in physical units.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLatticeBC_Sinusoid<EvalT, Traits>::
DDLatticeBC_Sinusoid(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names")) ;
  std::string prefix = p.get<std::string>("Prefix");

  // basis
  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");

  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(names.dof.phi);
  RCP<PHX::DataLayout> data_layout = basis->functional;
  num_basis = data_layout->dimension(1);

  //Set up to write the contact voltage to the parameter library
  contactVoltage = rcp(new panzer::ScalarParameterEntry<EvalT>);
  contactVoltage->setRealValue(0);
  contactVoltageName = p.get<std::string>("Sideset ID")+"_Voltage";
  contactVoltage = panzer::createAndRegisterScalarParameter<EvalT>(
						    std::string(contactVoltageName),
						    *p.get<RCP<panzer::ParamLib> >("ParamLib"));

  // read in user-specified values
  dc_offset = p.get<double>("DC Offset");
  amplitude1 = p.get<double>("Amplitude 1");
  frequency1 = p.get<double>("Frequency 1");
  phaseshift1 = p.get<double>("Phase Shift 1");
  amplitude2 = p.get<double>("Amplitude 2");
  frequency2 = p.get<double>("Frequency 2");
  phaseshift2 = p.get<double>("Phase Shift 2");

  bUseFD = false;
  if (p.isParameter("Fermi Dirac"))  bUseFD = p.get<bool>("Fermi Dirac");

  incmpl_ioniz = p.sublist("Incomplete Ionization");
  expandIonizEnParams(incmpl_ioniz);

  bSolveIon = p.get<bool>("Solve Ion");
  ion_charge = p.get<int>("Ion Charge");
  bUseFermiPin = p.get<bool>("Fermi Level Pinning");
  ionDens = p.get<double>("Contact Ion Density"); 

  //Set the parameter library contact voltage to an initial value
  contactVoltage->setValue(dc_offset);

  // evaluated fields
  potential = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.phi, data_layout);
  edensity = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.edensity, data_layout);
  hdensity = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.hdensity, data_layout);

  // add evaluated fields
  this->addEvaluatedField(potential);
  this->addEvaluatedField(edensity);
  this->addEvaluatedField(hdensity);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  T0 = scaleParams->scale_params.T0;
  t0 = scaleParams->scale_params.t0;
  C0 = scaleParams->scale_params.C0;

  // dependent fields
  doping = MDField<const ScalarT,Cell,BASIS>(names.field.doping_raw, data_layout);
  acceptor = MDField<const ScalarT,Cell,BASIS>(names.field.acceptor_raw, data_layout);
  donor = MDField<const ScalarT,Cell,BASIS>(names.field.donor_raw, data_layout);
  gamma_e = MDField<const ScalarT,Cell,BASIS>(names.field.elec_deg_factor, data_layout);
  gamma_h = MDField<const ScalarT,Cell,BASIS>(names.field.hole_deg_factor, data_layout);
  intrin_conc = MDField<const ScalarT,Cell,BASIS>(names.field.intrin_conc, data_layout);
  elec_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.elec_eff_dos, data_layout);
  hole_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.hole_eff_dos, data_layout);
  eff_affinity =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_affinity, data_layout);
  eff_bandgap =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_band_gap, data_layout);
  latt_temp =  MDField<const ScalarT,Cell,BASIS>(names.field.latt_temp, data_layout);
  ref_energy =  MDField<const ScalarT,Cell,BASIS>(names.field.ref_energy, data_layout);
  this->addDependentField(ref_energy);

  // add dependent fields
  this->addDependentField(doping);
  this->addDependentField(acceptor);
  this->addDependentField(donor);
  this->addDependentField(gamma_e);
  this->addDependentField(gamma_h);
  this->addDependentField(intrin_conc);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(latt_temp);

  if (bSolveIon)
  {
    ion_density = MDField<const ScalarT,Cell,BASIS>(names.dof.iondensity, data_layout);
    this->addDependentField(ion_density);
  }

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral =
    Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

  std::string n = "Ohmic Contact for DDLattice";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLatticeBC_Sinusoid<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  ScalarT Eref = ref_energy(0,0);

  // get the current (present) time in [s]
  double curr_time = workset.time * t0;

  // compute the voltage in [V]
  double pi = 4.0*atan(1.0);
  ScalarT voltage = dc_offset
                    + amplitude1 * std::sin(2.0*pi*frequency1*curr_time - 2.0*pi*phaseshift1)
                    + amplitude2 * std::sin(2.0*pi*frequency2*curr_time - 2.0*pi*phaseshift2);

  contactVoltage->setValue(voltage);

  // obtain kb
  charon::PhysicalConstants const& phyConst = charon::PhysicalConstants::Instance();
  double kbBoltz = phyConst.kb;   // Boltzmann constant in [eV/K]

  const bool withAccIncmplIoniz =
    (incmpl_ioniz.sublist("Acceptor").numParams() == 0) ? false : true;
  const bool withDonIncmplIoniz =
    (incmpl_ioniz.sublist("Donor").numParams() == 0) ? false : true;

 // use the Fermi-Dirac (FD) statistics. To compare FD with MB, we should choose
 // the right intrinsic conc. model so that nie = sqrt(Nc*Nv)*exp(-Egeff/2kbT) for MB
 if (bUseFD)
 {

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      // obtain concentrations (scaled)
      const ScalarT& dop = doping(cell,basis);
      const ScalarT& Nc = elec_effdos(cell,basis);
      const ScalarT& Nv = hole_effdos(cell,basis);
      const ScalarT& Na = acceptor(cell,basis);
      const ScalarT& Nd = donor(cell,basis);
      const ScalarT& nie = intrin_conc(cell,basis);
      const ScalarT& e_gamma = gamma_e(cell,basis);
      const ScalarT& h_gamma = gamma_h(cell,basis);

      // obtain energies in [eV]
      const ScalarT& effChi = eff_affinity(cell,basis);
      const ScalarT& effEg = eff_bandgap(cell,basis);

      // obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,basis)*T0;
      //if (lattT <= 0)
      //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
      //      "Error: Lattice temperature must be > 0 ! But its value is " << lattT << " ! \n");

      ScalarT kbT = kbBoltz * lattT;  // [eV]
      // ScalarT qtheta = effChi + 0.5*effEg + 0.5*kbT*log(Nc/Nv);  // [eV]

      // set Neff
      ScalarT Neff = dop;
      ScalarT Nion = 0.;
      if (bSolveIon)
      {
        if (bUseFermiPin)
          // Nion = 5.e20 / C0;  // hard-code 5.e20 for now ???
          Nion = ionDens / C0;   // ionDens is specified in the input
        else
          Nion = ion_density(cell,basis);
        Neff = dop + ion_charge * Nion;
      }

      // compute carrier density and potential (scaled) at ohmic contact
      if (Neff >= 0.0)  // n-type contact
      {
        ScalarT n0 = 0.;
        if(!withDonIncmplIoniz) {
          n0 = Neff;  // assume n0 = Neff
        } else {
          double ND_crit = incmpl_ioniz.sublist("Donor").
            get<double>("Critical Doping Value");
          if(Nd > ND_crit/C0)
            // dopant fully ionized
            n0 = Neff;
          else
            // with incomplete ionized donor
            n0 = compute_FD_carrier_dens(-1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(e_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(h_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        ScalarT Ef_minus_Ec = kbT * (*inverseFermiIntegral)(n0/Nc);  // [eV]
        ScalarT Ef_minus_Ev = Ef_minus_Ec + effEg;  // [eV]
        ScalarT p0 = Nv * exp(-Ef_minus_Ev/kbT);
        // ScalarT potBC = qtheta - effChi + Ef_minus_Ec + voltage;  // [eV]=[V]
        ScalarT potBC = Eref - effChi + Ef_minus_Ec + voltage;  // [eV]=[V]

        edensity(cell,basis) = n0;
        hdensity(cell,basis) = p0;
        potential(cell,basis) = potBC / V0;
      }
      else  // p-type contact
      {
        ScalarT p0 = 0.;
        if(!withAccIncmplIoniz) {
          p0 = -Neff;  // assume p0 = -Neff
        } else {
          double NA_crit = incmpl_ioniz.sublist("Acceptor").
            get<double>("Critical Doping Value");
          if(Na > NA_crit/C0)
            // dopant fully ionized
            p0 = -Neff;
          else
            // with incomplete ionized acceptor
            p0 = compute_FD_carrier_dens(1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(e_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(h_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        ScalarT Ev_minus_Ef = kbT * (*inverseFermiIntegral)(p0/Nv);  // [eV]
        ScalarT Ec_minus_Ef = Ev_minus_Ef + effEg;  // [eV]
        ScalarT n0 = Nc * exp(-Ec_minus_Ef/kbT);
        // ScalarT potBC = qtheta - effChi -effEg - Ev_minus_Ef + voltage;  // [eV]=[V]
        ScalarT potBC = Eref - effChi -effEg - Ev_minus_Ef + voltage;  // [eV]=[V]

        edensity(cell,basis) = n0;
        hdensity(cell,basis) = p0;
        potential(cell,basis) = potBC / V0;
      }

    }
  }
 }  // end of the if (bUseFD) block


 // use the Maxwell-Boltzmann (MB) statistics. The MB implementation assumes
 // nie = sqrt(Nc*Nv)*exp(-Egeff/2kbT).
 else
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      // obtain concentrations (scaled)
      const ScalarT& dop = doping(cell,basis);
      const ScalarT& Na = acceptor(cell,basis);
      const ScalarT& Nd = donor(cell,basis);
      const ScalarT& nie = intrin_conc(cell,basis);
      const ScalarT& Nc = elec_effdos(cell,basis);
      const ScalarT& Nv = hole_effdos(cell,basis);

      // obtain energies in [eV]
      const ScalarT& effChi = eff_affinity(cell,basis);
      const ScalarT& effEg = eff_bandgap(cell,basis);

      // obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,basis)*T0;
      //if (lattT <= 0)
      //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
      //      "Error: Lattice temperature must be > 0 ! But its value is " << lattT << " ! \n");

      ScalarT kbT = kbBoltz * lattT;  // [eV]
      // ScalarT qtheta = effChi + 0.5*effEg + 0.5*kbT*log(Nc/Nv);  // [eV]

      // set Neff
      ScalarT Neff = dop;
      ScalarT Nion = 0.;
      if (bSolveIon)
      {
        ScalarT Nion;
        if (bUseFermiPin)
          // Nion = 5.e20 / C0;  // hard-code 5.e20 for now ???
          Nion = ionDens / C0;   // ionDens is specified in the input
        else
          Nion = ion_density(cell,basis);
        Neff = dop + ion_charge * Nion;
      }

      ScalarT y, potBC;

      // compute carrier density and potential (scaled) at ohmic contact
      if (Neff >= 0.0)  // n-type contact
      {
        ScalarT n0 = 0.;
        if(!withDonIncmplIoniz) {
          n0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) + 0.5*Neff;
        } else {
          double ND_crit = incmpl_ioniz.sublist("Donor").
            get<double>("Critical Doping Value");
          if(Nd > ND_crit/C0)
            // dopant fully ionized
            n0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) + 0.5*Neff;
          else
            // donor incomplete ionization active
            // acceptor incomplete ionization can be active or not
            // depending on the user input
            n0 = compute_MB_carrier_dens(-1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        edensity(cell,basis) = n0;
        hdensity(cell,basis) = nie*nie/n0;

        y = 0.5*Neff/Nc + sqrt(pow(0.5*Neff/Nc,2.0) + Nv/Nc*exp(-effEg/kbT));
        // potBC = qtheta - effChi + kbT*log(y) + voltage;  // [V]
        potBC = Eref - effChi + kbT*log(y) + voltage;  // [V]
        potential(cell,basis) = potBC / V0;  // scaled
      }

      else  // p-type contact
      {
        ScalarT p0 = 0.;
        if(!withAccIncmplIoniz) {
          p0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) - 0.5*Neff;
        } else {
          double NA_crit = incmpl_ioniz.sublist("Acceptor").
            get<double>("Critical Doping Value");
          if(Na > NA_crit/C0)
            p0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) - 0.5*Neff;
          else
            // acceptor incomplete ionization active
            // donor incomplete ionization can be active or not
            // depending on the user input
            p0 = compute_MB_carrier_dens(1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        hdensity(cell,basis) = p0;
        edensity(cell,basis) = pow(nie,2.0)/p0;

        y = -0.5*Neff/Nv + sqrt(pow(0.5*Neff/Nv,2.0) + Nc/Nv*exp(-effEg/kbT));
        // potBC = qtheta - effChi - effEg - kbT*log(y) + voltage;  // [V]
        potBC = Eref - effChi - effEg - kbT*log(y) + voltage;
        potential(cell,basis) = potBC / V0;  // scaled
      }

    }
  }
 }  // end of the else block

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DDLatticeBC_Sinusoid<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<double>("DC Offset", 0.0);
  p->set<double>("Amplitude 1", 0.0);
  p->set<double>("Amplitude 2", 0.0);
  p->set<double>("Frequency 1", 0.0);
  p->set<double>("Frequency 2", 0.0);
  p->set<double>("Phase Shift 1", 0.0);
  p->set<double>("Phase Shift 2", 0.0);

  p->set<bool>("Fermi Dirac", false);
  p->sublist("Incomplete Ionization");
  p->sublist("Incomplete Ionization").sublist("Acceptor");
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Critical Doping Value", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Degeneracy Factor", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Ionization Energy", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<std::string>("Approximation", "None");
  p->sublist("Incomplete Ionization").sublist("Donor");
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Critical Doping Value", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Degeneracy Factor", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Ionization Energy", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<std::string>("Approximation", "None");

  p->set<bool>("Solve Ion", false);
  p->set<int>("Ion Charge", 1);
  p->set<bool>("Fermi Level Pinning", false);
  p->set<double>("Contact Ion Density", 0.0);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

 //Sideset ID
  p->set<std::string>("Sideset ID","");
 
  Teuchos::RCP<panzer::ParamLib> pl;
  p->set("ParamLib", pl);

  return p;
}

}

#endif

