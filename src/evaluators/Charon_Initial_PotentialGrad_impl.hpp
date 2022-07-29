
#ifndef CHARON_SAVE_GRADPOTENTIAL_IMPL_HPP
#define CHARON_SAVE_GRADPOTENTIAL_IMPL_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Charon_Names.hpp"


namespace charon {

template<typename EvalT, typename Traits>
Initial_PotentialGrad<EvalT, Traits>::
Initial_PotentialGrad(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using panzer::Dim; 

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  basis_name = basis->name();
  num_basis = data_layout->dimension(1);

  // phalanx fields
  potential = MDField<const ScalarT,Cell,BASIS>(n.dof.phi, data_layout);
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.phi, vector);

  this->addDependentField(grad_phi);
  this->addDependentField(potential); 

  initial_phi = MDField<ScalarT,Cell,BASIS>(n.field.initial_phi, data_layout); 
  initial_grad_phi = MDField<ScalarT,Cell,IP,Dim>(n.field.initial_grad_phi,vector);

  this->addEvaluatedField(initial_phi);
  this->addEvaluatedField(initial_grad_phi);

  std::string name = "Initial_PotentialGrad";
  this->setName(name);

  std::size_t num_wksts = p.get<std::size_t>("Max Worksets");
  initial_phi_wkst.resize(num_wksts);
  initial_grad_phi_wkst.resize(num_wksts);
  bSaveField_wkst.resize(num_wksts); 

  for (std::size_t i = 0;i < num_wksts; ++i)
  {
    std::string phi_name = n.field.initial_phi + std::to_string(i); 
    std::string grad_phi_name = n.field.initial_grad_phi + std::to_string(i);
    initial_phi_wkst[i] = MDField<ScalarT,Cell,BASIS>(phi_name,data_layout); 
    initial_grad_phi_wkst[i] = MDField<ScalarT,Cell,IP,Dim>(grad_phi_name,vector);
    this->addEvaluatedField(initial_phi_wkst[i]); 
    this->addEvaluatedField(initial_grad_phi_wkst[i]);
    bSaveField_wkst[i] = false; 
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Initial_PotentialGrad<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  using panzer::index_t;

  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT,typename Traits>
void Initial_PotentialGrad<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
  // initialize worksetId
  worksetId = 0;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Initial_PotentialGrad<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  using std::cout;
  using std::endl;
  using std::ofstream;
  using std::ios;
  using panzer::index_t;
  using Teuchos::RCP;
  using Teuchos::rcp;

  if (bSaveField_wkst[worksetId])  
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int basis = 0; basis < num_basis; ++basis)
        initial_phi(cell,basis) = initial_phi_wkst[worksetId](cell,basis);
 
      for (int ip = 0; ip < num_ip; ++ip)
        for (int dim = 0; dim < num_dim; ++dim)
          initial_grad_phi(cell,ip,dim) = initial_grad_phi_wkst[worksetId](cell,ip,dim); 
    }
  }

  else
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int basis = 0; basis < num_basis; ++basis)
        initial_phi_wkst[worksetId](cell,basis) = potential(cell,basis);

      for (int ip = 0; ip < num_ip; ++ip)
        for (int dim = 0; dim < num_dim; ++dim)
          initial_grad_phi_wkst[worksetId](cell,ip,dim) = grad_phi(cell,ip,dim); 
    }

    bSaveField_wkst[worksetId] = true;  
  }

  ++worksetId;
}


/////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
//////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Initial_PotentialGrad<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<std::size_t>("Max Worksets", 1);

  return p;
}



} // namespace charon

#endif

