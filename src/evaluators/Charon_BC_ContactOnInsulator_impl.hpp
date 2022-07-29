
#ifndef CHARON_BC_CONTACTONINSULATOR_IMPL_HPP
#define CHARON_BC_CONTACTONINSULATOR_IMPL_HPP

#include <cmath>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_ContactOnInsulator<EvalT, Traits>::
BC_ContactOnInsulator(
  const Teuchos::ParameterList& p)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string; 

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters(false);
  p.validateParameters(*valid_params);

  std::string prefix = p.get<std::string>("Prefix");
  m_names = p.get< RCP<const charon::Names> >("Names");
  // note that m_names never has a fd suffix, even if in a frequency domain simulation,
  // since the closure model is evaluated at the time collocation points
  // so, for frequency domain simulations, find basis from the zero-th harmonic
  Teuchos::RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(1,"","","","_CosH"+std::to_string(0.0)+"_"));

  const charon::Names& names = *m_names;

  // basis
  // This one works too
  // Teuchos::RCP<panzer::PureBasis> basis = p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");
  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(p.get<bool>("Frequency Domain") ? (*fd_names).dof.phi : (*m_names).dof.phi);

  RCP<PHX::DataLayout> data_layout = basis->functional;
  basis_name = basis->name();

  //Set up to write the contact voltage to the parameter library
  contactVoltage = rcp(new panzer::ScalarParameterEntry<EvalT>);
  contactVoltage->setRealValue(0);
  contactVoltageName = p.get<std::string>("Sideset ID")+"_Voltage";
  contactVoltage = 
    panzer::createAndRegisterScalarParameter<EvalT>(
						    std::string(contactVoltageName),
						    *p.get<RCP<panzer::ParamLib> >("ParamLib"));

  // read in user-specified voltage
  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(0);
  this->small_signal_perturbation = 0.0;
  if(p.get<bool>("Frequency Domain")){
    this->small_signal_perturbation = p.get<double>("Small Signal Perturbation");
  }

  if (p.isParameter("Enable Linear Ramp"))
    bLinRamp = p.get<bool>("Enable Linear Ramp");   

  if (p.isParameter("Enable Trapezoid Pulse"))
  {
    bTrapezoid = p.get<bool>("Enable Trapezoid Pulse");
  }
 
  if (p.isType<double>("Voltage"))
    user_value->setRealValue(p.get<double>("Voltage"));

  // Varying voltage
  else if (p.isType<std::string>("Varying Voltage"))
  {
    if (p.get<string>("Varying Voltage") == "Parameter")
    {
	user_value =
	  panzer::createAndRegisterScalarParameter<EvalT>(
							  std::string("Varying Voltage"),
							  *p.get<RCP<panzer::ParamLib> >("ParamLib"));
	Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > IAmLoca;
	std::string locaFlag = contactVoltageName+"IS LOCA";
	IAmLoca =
	  panzer::createAndRegisterScalarParameter<EvalT>(
							  std::string(locaFlag),
							  *p.get<RCP<panzer::ParamLib> >("ParamLib"));
	IAmLoca->setValue(1.0);

        if (p.isParameter("Initial Voltage"))
          initial_voltage = p.get<double>("Initial Voltage"); 
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "BC_ContactOnInsulator():  Error:  Expecting Varying Voltage value "  \
        "of \"Parameter\"; received \""
        << p.get<std::string>("Varying Voltage") << "\".")
  }

  // Time-dependent linear ramp voltage source
  // else if (p.isSublist("Linear Ramp ParameterList"))
  else if (bLinRamp)
  {
    bLinRamp = true; 
    const ParameterList& linRampPL = *(p.get< RCP<ParameterList> >("Linear Ramp ParameterList"));
    initial_time = linRampPL.get<double>("Initial Time");
    initial_voltage = linRampPL.get<double>("Initial Voltage");
    final_time = linRampPL.get<double>("Final Time");
    final_voltage = linRampPL.get<double>("Final Voltage");

    // compute slope [V/s] and intercept [V] of the linear ramp
    slope = (initial_voltage - final_voltage) / (initial_time - final_time);
    yintercept = -slope * initial_time + initial_voltage;
  }

  else if (bTrapezoid)
  {
    bTrapezoid = true; 
    const ParameterList& trapzPulsePL = *(p.get< RCP<ParameterList> >("Trapezoid Pulse ParameterList"));
    dc_offset = trapzPulsePL.get<double>("DC Offset");
    amplitude = trapzPulsePL.get<double>("Amplitude");
    period = trapzPulsePL.get<double>("Period");
    rise_time = trapzPulsePL.get<double>("Rise Time");
    fall_time = trapzPulsePL.get<double>("Fall Time");
    delay = trapzPulsePL.get<double>("Delay");
    duty_cycle = trapzPulsePL.get<double>("Duty Cycle");
    num_pulses = trapzPulsePL.get<int>("Number Pulses");
    user_value->setValue(dc_offset);

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
  }


  //Set the parameter library contact voltage to an initial value
  contactVoltage->setValue(user_value->getValue());

  work_func = p.get<double>("Work Function");

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  t0 = scaleParams->scale_params.t0;   

  // fields
  potential = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.phi, data_layout);
  ref_energy = MDField<const ScalarT>(names.field.ref_energy, data_layout);

  this->addEvaluatedField(potential);
  this->addDependentField(ref_energy);

  std::string n = "BC at Contact On Insulator";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_ContactOnInsulator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getPureBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_ContactOnInsulator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  ScalarT voltage; 

  if (bLinRamp)  // linear ramp voltage source
  { 
    // retrieve time scaling factor in [s]
    double tscale = Sacado::ScalarValue<ScalarT>::eval(t0);

    // get the current (present) time in [s]
    double curr_time = workset.time * tscale;

    // compute the voltage in [V]
    double volt = 0.0;
    if (curr_time <= initial_time)
      volt = initial_voltage;
    else if (curr_time > final_time)
      volt = final_voltage;
    else
      volt = slope * curr_time + yintercept;

    // assign voltage
    voltage = volt; 
    
  }
  else if (bTrapezoid)
  {
    // retrieve time scaling factor in [s]
    double tscale = Sacado::ScalarValue<ScalarT>::eval(t0);

    // get the current (present) time in [s]
    double curr_time = workset.time * tscale;

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
    voltage = volt;
  }

  else  // for Voltage or Varying Voltage
    voltage = user_value->getValue() + this->small_signal_perturbation + initial_voltage;

  //For output later
  contactVoltage->setValue(voltage);
 
  ScalarT offsetDueToWF = (work_func - ref_energy(0,0))/1.0;  // 1.0 converts from [eV] to [V]
  ScalarT bcValue = (voltage - offsetDueToWF)/V0;

  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  using panzer::index_t;
  size_type num_basis = potential.dimension(1);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (size_type basis = 0; basis < num_basis; ++basis)
    {
      potential(cell,basis) = bcValue;
    }
  }

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BC_ContactOnInsulator<EvalT, Traits>::getValidParameters(bool stochasticParams) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector; 

  RCP<ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<string>("Prefix", "?");

  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  RCP<const charon::Names> n;
  p->set("Names", n);

  p->set("Frequency Domain", false);

  if(stochasticParams)
    p->set<string>("Voltage", "0.0");
  else
    p->set<double>("Voltage", 0.0);

  p->set<string>("Varying Voltage", "Parameter");
  p->set<double>("Small Signal Perturbation", 0.0);
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));

  p->set<double>("Initial Voltage", 0.0); // to be used with Varying Voltage
  p->set<double>("Work Function", 0.0);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set<bool>("Enable Linear Ramp", false); 

  RCP<ParameterList> linRampPL = rcp(new ParameterList);
  p->set("Linear Ramp ParameterList", linRampPL);

  linRampPL->set<double>("Initial Time", 0.0, "Initial time in (s)");
  linRampPL->set<double>("Final Time", 0.0, "Final time in (s)");
  linRampPL->set<double>("Initial Voltage", 0.0, "Initial voltage in (V)");
  linRampPL->set<double>("Final Voltage", 0.0, "Final voltage in (V)"); 

  p->set<bool>("Enable Trapezoid Pulse", false); 

  RCP<ParameterList> trapzPulsePL = rcp(new ParameterList);
  p->set("Trapezoid Pulse ParameterList", trapzPulsePL);
  
  trapzPulsePL->set<double>("DC Offset",0.0);
  trapzPulsePL->set<double>("Amplitude",0.0);
  trapzPulsePL->set<double>("Period",0.0);
  trapzPulsePL->set<double>("Rise Time",0.0);
  trapzPulsePL->set<double>("Fall Time",0.0);
  trapzPulsePL->set<double>("Delay",0.0);
  trapzPulsePL->set<double>("Duty Cycle",1.0);
  trapzPulsePL->set<int>("Number Pulses",1);

  //Sideset ID
  p->set<std::string>("Sideset ID","");

  return p;
}

}

#endif

