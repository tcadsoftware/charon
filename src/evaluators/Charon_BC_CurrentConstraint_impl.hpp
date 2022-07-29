
#ifndef   Charon_BC_CurrentConstraint_impl_hpp
#define   Charon_BC_CurrentConstraint_impl_hpp

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_BC_OhmicContact.hpp"
#include "Charon_Util.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

// Teuchos
#include "Teuchos_TestForException.hpp"

// Panzer
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

namespace charon
{

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_CurrentConstraint<EvalT, Traits>::
BC_CurrentConstraint(
  const Teuchos::ParameterList& p)
{
  using charon::Names;
  using charon::Scaling_Parameters;
  using panzer::FieldLibraryBase;
  using panzer::PureBasis;
  using panzer::ScalarParameterEntry;
  using PHX::DataLayout;
  using std::runtime_error;
  using std::string;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Validate the input ParameterList.
  RCP<ParameterList> validParams = this->getValidParameters();
  p.validateParameters(*validParams);

  // Read in the Voltage Control.
  voltageParameter_ =
    p.get<RCP<ScalarParameterEntry<EvalT>>>("Voltage Control");
  TEUCHOS_TEST_FOR_EXCEPTION(voltageParameter_ == null, runtime_error,
    "Error:  \"Voltage Control\" is null in BC_CurrentConstraint.");

  // Set the contactvoltage in the parameter library
  //Set up to write the contact voltage to the parameter library
  contactVoltageName = p.get<std::string>("Sideset ID")+"_Voltage";
  contactVoltage = 
  panzer::createAndRegisterScalarParameter<EvalT>(
  						    std::string(contactVoltageName),
  						    *p.get<RCP<panzer::ParamLib> >("ParamLib"));
  contactVoltage->setValue(voltageParameter_->getValue());

  // Determine whether or not we're using Fermi-Dirac.
  bUseFD_ = false;
  if (p.isParameter("Fermi Dirac"))
    bUseFD_ = p.get<bool>("Fermi Dirac");

  // Determine whether or not we'll be enforcing a BJT1D Base Contact
  // boundary condition.
  bjt1DBaseContact_ = false;
  if (p.isParameter("BJT1D Base Contact"))
    bjt1DBaseContact_ = p.get<bool>("BJT1D Base Contact");

  // Determine whether or not to use Reference Energy
  bUseRefE_ = true;
  if (p.isParameter("Use Reference Energy"))
    bUseRefE_ = p.get<bool>("Use Reference Energy");

  // Get the Incomplete Ionization parameters.
  incmplIoniz_ = p.sublist("Incomplete Ionization");
  expandIonizEnParams(incmplIoniz_);

  // Get the DataLayout for the basis.
  RCP<const FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const FieldLibraryBase>>("Field Library");
  const Names& names = *(p.get<RCP<const Names>>("Names"));
  RCP<const PureBasis> basis = fieldLayoutLibrary->lookupBasis(names.dof.phi);
  RCP<DataLayout> dataLayout = basis->functional;

  // Create evaluated fields.
  string prefix(p.get<string>("Prefix"));
  potential_ = EvaluatedField(prefix + names.dof.phi,      dataLayout);
  eDensity_  = EvaluatedField(prefix + names.dof.edensity, dataLayout);
  hDensity_  = EvaluatedField(prefix + names.dof.hdensity, dataLayout);

  // Add evaluated fields.
  this->addEvaluatedField(potential_);
  this->addEvaluatedField(eDensity_);
  this->addEvaluatedField(hDensity_);

  // scaling parameters
  scaleParams_ = p.get<RCP<Scaling_Parameters>>("Scaling Parameters");
  V0_ = scaleParams_->scale_params.V0;
  C0_ = scaleParams_->scale_params.C0;
  T0_ = scaleParams_->scale_params.T0;

  // Create dependent fields.
  doping_      = DependentField(names.field.doping_raw,      dataLayout);
  acceptor_    = DependentField(names.field.acceptor_raw,    dataLayout);
  donor_       = DependentField(names.field.donor_raw,       dataLayout);
  intrinConc_  = DependentField(names.field.intrin_conc,     dataLayout);
  eEffDos_     = DependentField(names.field.elec_eff_dos,    dataLayout);
  hEffDos_     = DependentField(names.field.hole_eff_dos,    dataLayout);
  effAffinity_ = DependentField(names.field.eff_affinity,    dataLayout);
  effBandgap_  = DependentField(names.field.eff_band_gap,    dataLayout);
  lattTemp_    = DependentField(names.field.latt_temp,       dataLayout);

  // Add dependent fields.
  this->addDependentField(doping_);
  this->addDependentField(acceptor_);
  this->addDependentField(donor_);
  this->addDependentField(intrinConc_);
  this->addDependentField(eEffDos_);
  this->addDependentField(hEffDos_);
  this->addDependentField(effAffinity_);
  this->addDependentField(effBandgap_);
  this->addDependentField(lattTemp_);

  if (bUseRefE_)
  {
    refEnergy_   = DependentField(names.field.ref_energy, dataLayout);
    this->addDependentField(refEnergy_);
  }

  // Set the name of this object.
  string n("Current Constraint Contact");
  this->setName(n);
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_CurrentConstraint<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData /* sd */,
  PHX::FieldManager<Traits>& /* fm */)
{
} // end of postRegistrationSetup()

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_CurrentConstraint<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  // Get the voltage value that is inversely computed.
  ScalarT voltage     = voltageParameter_->getValue();
  contactVoltage->setValue(voltage);
  ScalarT vScaling    = V0_;
  ScalarT densScaling = C0_;
  ScalarT tempScaling = T0_;

  ScalarT eRef = 0.0;
  if (bUseRefE_)  eRef = refEnergy_(0, 0);

  // Evaluate the Ohmic Contact boundary condition.
  OhmicContact<EvalT, Traits>::evaluateOhmicContactBC(
    bjt1DBaseContact_, bUseFD_, bUseRefE_, incmplIoniz_, voltage, eRef, vScaling,
    densScaling, tempScaling, workset, doping_, acceptor_, donor_,
    intrinConc_, eEffDos_, hEffDos_, effAffinity_, effBandgap_,
    lattTemp_, potential_, eDensity_, hDensity_);
} // end of evaluateFields()

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BC_CurrentConstraint<EvalT, Traits>::
getValidParameters() const
{
  using charon::Names;
  using charon::Scaling_Parameters;
  using panzer::FieldLibraryBase;
  using panzer::ScalarParameterEntry;
  using std::string;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  RCP<ParameterList> p = rcp(new ParameterList);
  p->set<string>("Prefix", "");
  RCP<const FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);
  RCP<const Names> n;
  p->set("Names", n);
  p->set<RCP<ScalarParameterEntry<EvalT>>>("Voltage Control", null);
  p->set<bool>("Fermi Dirac",        false);
  p->set<bool>("BJT1D Base Contact", false);
  p->set<bool>("Use Reference Energy", true);
  RCP<Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);
  p->set("Sideset ID","");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));
  RCP<charon::EmpiricalDamage_Data> dmgdata;
  p->set("empirical damage data", dmgdata);

  ParameterList& incompleteIonizationPL = p->sublist("Incomplete Ionization");
  ParameterList& acceptorPL = incompleteIonizationPL.sublist("Acceptor");
  acceptorPL.set<double>("Critical Doping Value", 0.0   );
  acceptorPL.set<double>("Degeneracy Factor",     0.0   );
  acceptorPL.set<double>("Ionization Energy",     0.0   );
  acceptorPL.set<string>("AccIncmplIoniz File",   ""    );
  acceptorPL.set<string>("Approximation",         "None");
  ParameterList& donorPL = incompleteIonizationPL.sublist("Donor");
  donorPL.set<double>("Critical Doping Value", 0.0   );
  donorPL.set<double>("Degeneracy Factor",     0.0   );
  donorPL.set<double>("Ionization Energy",     0.0   );
  donorPL.set<string>("DonIncmplIoniz File",   ""    );
  donorPL.set<string>("Approximation",         "None");

  return p;
} // end of getValidParameters()

} // end of namespace charon

#endif // Charon_BC_CurrentConstraint_impl_hpp
