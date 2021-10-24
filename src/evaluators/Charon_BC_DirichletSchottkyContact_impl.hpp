
#ifndef CHARON_BC_DIRICHLETSCHOTTKYCONTACT_IMPL_HPP
#define CHARON_BC_DIRICHLETSCHOTTKYCONTACT_IMPL_HPP

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
#include "Charon_BC_DirichletSchottkyContact.hpp"
#include "Charon_Util.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_DirichletSchottkyContact<EvalT, Traits>::
BC_DirichletSchottkyContact(
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

  // read in user-specified voltage
  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(0);
  if (p.isType<double>("Voltage")) {
    user_value->setRealValue(p.get<double>("Voltage"));
  } else if (p.isType<std::string>("Varying Voltage")) {        
    if (p.get<std::string>("Varying Voltage") == "Parameter")
      user_value =
        panzer::createAndRegisterScalarParameter<EvalT>(
        std::string("Varying Voltage"),
        *p.get<RCP<panzer::ParamLib> >("ParamLib"));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "BC_SchottkyContact():  Error:  Expecting Varying Voltage value of " \
        "\"Parameter\"; received \"" << p.get<std::string>("Varying Voltage")
        << "\".")
  }

  // read in the workfunction
  cnt_wf = p.get<double>("Work Function");
  if(cnt_wf <= 0) 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
	 "BC_SchottkyContact():  Error:  Work Function parameter must be positive! " )

  // evaluated fields
  potential = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.phi, data_layout);

  // add evaluated fields
  this->addEvaluatedField(potential);
  
  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  
  // dependent fields
  ref_energy = MDField<const ScalarT,Cell,BASIS>(names.field.ref_energy, data_layout);

  // add dependent fields
  this->addDependentField(ref_energy);

  std::string n = "BC Dirichlet at Schottky Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_DirichletSchottkyContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  ScalarT voltage = user_value->getValue();
  ScalarT WF = cnt_wf;
  ScalarT Eref = ref_energy(0,0);
  ScalarT vScaling = V0;

  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (size_type basis = 0; basis < num_basis; ++basis) {
      potential(cell,basis) = (Eref - WF + voltage)/vScaling;
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
BC_DirichletSchottkyContact<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<double>("Voltage", 0.0);
  p->set<std::string>("Varying Voltage", "Parameter");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));

  p->set<double>("Work Function", 0.0);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

