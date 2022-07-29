
#ifndef CHARON_BC_BJT1DBASECONTACT_IMPL_HPP
#define CHARON_BC_BJT1DBASECONTACT_IMPL_HPP

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
#include "Charon_BC_OhmicContact.hpp"
#include "Charon_Util.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_BJT1DBaseContact<EvalT, Traits>::
BC_BJT1DBaseContact(
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
  // note that m_names never has a fd suffix, even if in a frequency domain simulation,
  // since the closure model is evaluated at the time collocation points
  // so, for frequency domain simulations, find basis from the zero-th harmonic
  Teuchos::RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(3,"","","","_CosH"+std::to_string(0.0)+"_"));

  // basis
  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");

  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(p.get<bool>("Frequency Domain") ? (*fd_names).dof.phi : (*m_names).dof.phi);

  RCP<PHX::DataLayout> data_layout = basis->functional;


  //Set up to write the contact voltage to the parameter library
  contactVoltage = rcp(new panzer::ScalarParameterEntry<EvalT>);
  contactVoltage->setRealValue(0);
  contactVoltageName = p.get<std::string>("Sideset ID")+"_Voltage";
  contactVoltage = 
    panzer::createAndRegisterScalarParameter<EvalT>(
						    std::string(contactVoltageName),
						    *p.get<RCP<panzer::ParamLib> >("ParamLib"));

  // read in user-specified voltage
  user_value = Teuchos::rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(0);
  if (p.isType<double>("Voltage"))
    user_value->setRealValue(p.get<double>("Voltage"));
  else if (p.isType<std::string>("Varying Voltage"))
  {
    if (p.get<std::string>("Varying Voltage") == "Parameter")
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
      }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "BC_BJT1DBaseContact():  Error:  Expecting Varying Voltage value of " \
        "\"Parameter\"; received \"" << p.get<std::string>("Varying Voltage")
        << "\".")
  }

  //Set the parameter library contact voltage to an initial value
  contactVoltage->setValue(user_value->getValue());

  base_doptype = p.get<std::string>("Base Doping Type");
  bUseFD = false;
  if (p.isParameter("Fermi Dirac"))  bUseFD = p.get<bool>("Fermi Dirac");

  incmpl_ioniz = p.sublist("Incomplete Ionization");
  expandIonizEnParams(incmpl_ioniz);

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

  std::string n = "BC at BJT1D Base Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_BJT1DBaseContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  ScalarT voltage = user_value->getValue();
  contactVoltage->setValue(voltage);
   ScalarT Eref = ref_energy(0,0);
  ScalarT vScaling = V0;
  ScalarT densScaling = C0;
  ScalarT tempScaling = T0;
  bool bBJT1DBase = true;
  bool bUseRefE = true;

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
BC_BJT1DBaseContact<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);
  p->set<bool>("Frequency Domain", false);

  p->set<double>("Voltage", 0.0);
  p->set<std::string>("Varying Voltage", "Parameter");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));
  p->set<std::string>("Base Doping Type", "??");
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

  Teuchos::RCP<charon::EmpiricalDamage_Data> dmgdata;
  p->set("empirical damage data", dmgdata);

   //Sideset ID
  p->set<std::string>("Sideset ID","");

  return p;
}

}

#endif

