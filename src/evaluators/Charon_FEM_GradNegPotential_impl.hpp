
#ifndef CHARON_FEM_GRADNEGPOTENTIAL_IMPL_HPP
#define CHARON_FEM_GRADNEGPOTENTIAL_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"


// Assign the x and y components of the negative potential gradient to scalar fields
// that can be output to Exodus

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
FEM_GradNegPotential<EvalT, Traits>::
FEM_GradNegPotential(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<panzer::IntegrationRule> ir =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  // Evaluated fields
  grad_negpot_x = MDField<ScalarT,Cell,Point>(n.field.grad_negpot_x, scalar);
  grad_negpot_y = MDField<ScalarT,Cell,Point>(n.field.grad_negpot_y, scalar);
  this->addEvaluatedField(grad_negpot_x);
  this->addEvaluatedField(grad_negpot_y);

  // Dependent fields
  grad_potential = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi,vector);
  this->addDependentField(grad_potential);

  std::string name = "FEM_GradNegPotential";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_GradNegPotential<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      grad_negpot_x(cell,point) = -grad_potential(cell,point,0);
      grad_negpot_y(cell,point) = -grad_potential(cell,point,1);
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
FEM_GradNegPotential<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;
}

}

#endif

