
#ifndef CHARON_BC_THERMALCONTACT_IMPL_HPP
#define CHARON_BC_THERMALCONTACT_IMPL_HPP

#include <cmath>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_String_Utilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"


// This Dirichlet Thermal Contact is to set Lattice Temperature to a user-specified
// value. Physically, it represents a perfect heat sink.

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_ThermalContact<EvalT, Traits>::
BC_ThermalContact(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  std::string prefix = p.get<std::string>("Prefix");

  // RCP<const charon::Names> m_names = p.get< RCP<const charon::Names> >("Names");
  // const charon::Names& names = *m_names;
  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // get basis and its data layout
  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");
  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(names.dof.latt_temp);

  RCP<DataLayout> data_layout = basis->functional;
  num_basis = data_layout->dimension(1);

  // read in user-specified temperature in [K]
  user_value = p.get<double>("Temperature");

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // fields
  latt_temp = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.latt_temp, data_layout);
  this->addEvaluatedField(latt_temp);

  std::string n = "BC at Thermal Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_ThermalContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  ScalarT bcValue = user_value / T0;  // scaled

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    for (int basis = 0; basis < num_basis; ++basis)
      latt_temp(cell,basis) = bcValue;

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BC_ThermalContact<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "?");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<double>("Temperature", 0.0);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

