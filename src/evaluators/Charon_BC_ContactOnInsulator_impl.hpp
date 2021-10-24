
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
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters(false);
  p.validateParameters(*valid_params);

  std::string prefix = p.get<std::string>("Prefix");
  m_names = p.get< RCP<const charon::Names> >("Names");
  // note that m_names never has a fd suffix, even if in a frequency domain simulation,
  // since the closure model is evaluated at the time collocation points
  // so, for frequency domain simulations, find basis from the zero-th harmonic
  Teuchos::RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(1,"","","","_CosH0.000000_"));

  const charon::Names& names = *m_names;

  // basis
  // This one works too
  // Teuchos::RCP<panzer::PureBasis> basis = p.get< Teuchos::RCP<panzer::PureBasis> >("Basis");

  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");
  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(p.get<bool>("Frequency Domain") ? (*fd_names).dof.phi : (*m_names).dof.phi);

  RCP<PHX::DataLayout> data_layout = basis->functional;
  basis_name = basis->name();

  // read in user-specified voltage
  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(0);
  if (p.isType<double>("Voltage"))
    user_value->setRealValue(p.get<double>("Voltage"));
  else if (p.isType<std::string>("Varying Voltage"))
  {
    if (p.get<std::string>("Varying Voltage") == "Parameter")
      user_value =
        panzer::createAndRegisterScalarParameter<EvalT>(
        std::string("Varying Voltage"),
        *p.get<RCP<panzer::ParamLib> >("ParamLib"));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "BC_ContactOnInsulator():  Error:  Expecting Varying Voltage value "  \
        "of \"Parameter\"; received \""
        << p.get<std::string>("Varying Voltage") << "\".")
  }

  work_func = p.get<double>("Work Function");

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;

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
  ScalarT voltage = user_value->getValue();

  ScalarT offsetDueToWF = (work_func - ref_energy(0,0))/1.0;  // 1.0 converts from [eV] to [V]
  ScalarT bcValue = (voltage - offsetDueToWF)/V0;
  //std::cout << "Eref=" << ref_energy(0,0) << ",WorkFunc=" << work_func
  //  << ",bcValue=" << bcValue << std::endl;

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
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "?");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set("Frequency Domain", false);

  if(stochasticParams)
    p->set<std::string>("Voltage", "0.0");
  else
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

