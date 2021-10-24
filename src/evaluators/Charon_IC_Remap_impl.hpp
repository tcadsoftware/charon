
#ifndef CHARON_IC_REMAP_IMPL_HPP
#define CHARON_IC_REMAP_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IC_Remap<EvalT, Traits>::
IC_Remap(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;

  dof_name = p.get<string>("DOF Name");
  dof_name_input = p.get<string>("Input DOF Name");

  output_field = MDField<ScalarT,Cell,BASIS>(dof_name, data_layout);
  input_field = MDField<const ScalarT,Cell,BASIS>(dof_name_input, data_layout);

  this->addEvaluatedField(output_field);
  this->addDependentField(input_field);

  std::string name = "IC_Remap";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Remap<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  size_type num_basis = input_field.dimension(1);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (size_type basis = 0; basis < num_basis; ++basis)
    {
      const ScalarT& val = input_field(cell,basis);
      output_field(cell,basis) = val;
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
IC_Remap<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("DOF Name", "?");
  p->set<std::string>("Input DOF Name", "?");

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  return p;
}

}

#endif

