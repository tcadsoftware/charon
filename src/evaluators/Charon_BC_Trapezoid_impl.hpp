
#ifndef CHARON_BC_TRAPEZOID_IMPL_HPP
#define CHARON_BC_TRAPEZOID_IMPL_HPP

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
#include "Charon_BC_OhmicContact.hpp"
#include "Charon_Util.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_Trapezoid<EvalT, Traits>::
BC_Trapezoid(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  std::string prefix = p.get<std::string>("Prefix");

  m_names = p.get< RCP<const charon::Names> >("Names");
  const charon::Names& names = *m_names;

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
  contactVoltage = 
    panzer::createAndRegisterScalarParameter<EvalT>(
						    std::string(contactVoltageName),
						    *p.get<RCP<panzer::ParamLib> >("ParamLib"));

 // read in user-specified values
  dc_offset = p.get<double>("DC Offset");
  amplitude = p.get<double>("Amplitude");
  period = p.get<double>("Period");
  rise_time = p.get<double>("Rise Time");
  fall_time = p.get<double>("Fall Time");
  delay = p.get<double>("Delay");
  duty_cycle = 1.0;
  duty_cycle = p.get<double>("Duty Cycle");
  num_pulses = 1;
  num_pulses = p.get<int>("Number Pulses");

  // compute slope [V/s] for rise and fall
  rising_slope = amplitude / rise_time;
  falling_slope = -amplitude / fall_time;

  // compute rise and fall time y axis intercepts
  // and times of rise, plateau and fall
  y_rise_intercept = dc_offset;
  pulse_time = period * duty_cycle;
  t1 = rise_time;
  peak_time = pulse_time - (rise_time + fall_time);
  t2 = t1 + peak_time;
  y_fall_intercept = -falling_slope * t2 + (amplitude+dc_offset);
  t3 = t2 + fall_time;
  if (peak_time < 0)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
        "Error: period * duty cycle - (rise time + fall time) >= 0 ! But its value is " << peak_time << " ! \n");

  bUseFD = false;
  if (p.isParameter("Fermi Dirac"))  bUseFD = p.get<bool>("Fermi Dirac");

  incmpl_ioniz = p.sublist("Incomplete Ionization");
  expandIonizEnParams(incmpl_ioniz);

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
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;
  t0 = scaleParams->scale_params.t0;

  // dependent fields
  doping = MDField<const ScalarT,Cell,BASIS>(names.field.doping_raw, data_layout);
  acceptor = MDField<const ScalarT,Cell,BASIS>(names.field.acceptor_raw, data_layout);
  donor = MDField<const ScalarT,Cell,BASIS>(names.field.donor_raw, data_layout);
  intrin_conc = MDField<const ScalarT,Cell,BASIS>(names.field.intrin_conc, data_layout);
  elec_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.elec_eff_dos, data_layout);
  hole_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.hole_eff_dos, data_layout);
  eff_affinity =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_affinity, data_layout);
  eff_bandgap =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_band_gap, data_layout);
  latt_temp =  MDField<const ScalarT,Cell,BASIS>(names.field.latt_temp, data_layout);
  ref_energy =  MDField<const ScalarT,Cell,BASIS>(names.field.ref_energy, data_layout);

  // add dependent fields
  this->addDependentField(doping);
  this->addDependentField(acceptor);
  this->addDependentField(donor);
  this->addDependentField(intrin_conc);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(latt_temp);
  this->addDependentField(ref_energy);

  std::string n = "BC at Trapezoidal Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_Trapezoid<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  // retrieve time scaling factor in [s]
  double tscale = Sacado::ScalarValue<ScalarT>::eval(t0);

  // get the current (present) time in [s]
  double curr_time = workset.time * tscale;

  // compute the voltage in [V]
  auto volt = 0.0;
  curr_time = curr_time - delay;
  volt = dc_offset;
  int int_time = static_cast<int>(curr_time/period);
  if (int_time < num_pulses)
    curr_time = curr_time - period*int_time;
  if (curr_time > 0)
  {
    if (curr_time <= t1)
      volt = rising_slope*curr_time + y_rise_intercept;
    else if (curr_time <= t2)
      volt = dc_offset + amplitude;
    else if (curr_time <= t3)
      volt = falling_slope*curr_time + y_fall_intercept;
    else
    {
      volt = dc_offset;
    }
  }
  
  ScalarT voltage = volt;

  ScalarT Eref = ref_energy(0,0);
  ScalarT vScaling = V0;
  ScalarT densScaling = C0;
  ScalarT tempScaling = T0;
  bool bBJT1DBase = false;
  bool bUseRefE = true;

  contactVoltage->setValue(voltage);

  OhmicContact<EvalT, Traits>::evaluateOhmicContactBC(
       bBJT1DBase, bUseFD, bUseRefE, incmpl_ioniz, voltage, Eref, vScaling, densScaling,
       tempScaling, workset, doping, acceptor, donor,
       intrin_conc, elec_effdos, hole_effdos, eff_affinity, eff_bandgap,
       latt_temp, potential, edensity, hdensity);
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BC_Trapezoid<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<double>("DC Offset", 0.0);
  p->set<double>("Amplitude", 0.0);
  p->set<double>("Period", 0.0);
  p->set<double>("Rise Time", 0.0);
  p->set<double>("Fall Time", 0.0);
  p->set<double>("Delay", 0.0);
  p->set<double>("Duty Cycle", 1.0);
  p->set<int>("Number Pulses", 1);
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
    set<std::string>("AccIncmplIoniz File", "");
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
    set<std::string>("DonIncmplIoniz File", "");
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<std::string>("Approximation", "None");

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

