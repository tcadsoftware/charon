
#ifndef CHARON_SUBTRACT_IMPL_HPP
#define CHARON_SUBTRACT_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Subtract<EvalT, Traits>::
Subtract(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // Retrieve the data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  diff = PHX::MDField<ScalarT,Cell,Point>(p.get<std::string>("Difference Name"), scalar);
  this->addEvaluatedField(diff);

  // Check which value(s) is given
  enableA = true;
  enableB = true;

  if (p.get<string>("Value A") == "") enableA = false;
  if (p.get<string>("Value B") == "") enableB = false;

  if (enableA)
  {
    valueA = MDField<const ScalarT,Cell,Point>(p.get<std::string>("Value A"), scalar);
    this->addDependentField(valueA);
  }

  if (enableB)
  {
    valueB = MDField<const ScalarT,Cell,Point>(p.get<std::string>("Value B"), scalar);
    this->addDependentField(valueB);
  }

  std::string name = "Subtract: Value A - Value B";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Subtract<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

 // when Value A and Value B are both given
 if (enableA && enableB)
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      diff(cell,point) = valueA(cell,point) - valueB(cell,point);
      // std::cout << "c=" << cell << ", p=" << point << std::setprecision(20) << ", A=" << valueA(cell,point) << ", B=" << valueB(cell,point) <<", diff=" << diff(cell,point) << std::endl;
    }
  }
 }

 // when only Value A is given
 else if (enableA && !enableB)
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
      diff(cell,point) = valueA(cell,point);
  }
 }

 // when only Value B is given
 else if (!enableA && enableB)
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
      diff(cell,point) = - valueB(cell,point);
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
Subtract<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Difference Name", "?");
  p->set<std::string>("Value A", "?");
  p->set<std::string>("Value B", "?");

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  return p;
}

}

#endif

