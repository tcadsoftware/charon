
#ifndef CHARON_SPACE_CHARGE_IMPL_HPP
#define CHARON_SPACE_CHARGE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Space_Charge<EvalT, Traits>::
Space_Charge(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // set evaluated field
  space_charge = MDField<ScalarT,Cell,Point>(n.field.space_charge,scalar);
  this->addEvaluatedField(space_charge);

  // add dependent fields
  elec_density = MDField<ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hole_density = MDField<ScalarT,Cell,Point>(n.dof.hdensity,scalar);
  doping = MDField<ScalarT,Cell,Point>(n.field.doping,scalar);

  this->addDependentField(elec_density);
  this->addDependentField(hole_density);
  this->addDependentField(doping);
  
  std::string name = "Space_Charge";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Space_Charge<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    for (int pt = 0; pt < num_points; ++pt)
      space_charge(cell,pt) = hole_density(cell,pt) - elec_density(cell,pt) + doping(cell,pt);
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Space_Charge<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  return p;
}

}

#endif
