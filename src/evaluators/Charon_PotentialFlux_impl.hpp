
#ifndef CHARON_POTENTIALFLUX_IMPL_HPP
#define CHARON_POTENTIALFLUX_IMPL_HPP

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
PotentialFlux<EvalT, Traits>::
PotentialFlux(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  RCP<panzer::IntegrationRule> ir =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");

  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  phi_flux = MDField<ScalarT,Cell,IP,Dim>(p.get<string>("Flux Name"),vector);
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(p.get<string>("Gradient Name"),vector);

  const RCP<const charon::Names> n = p.get< RCP<const charon::Names> >("Names");
  rel_perm = MDField<const ScalarT,Cell,IP>(n->field.rel_perm,scalar);

  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Lambda2 = scaleParams->scale_params.Lambda2;

  this->addEvaluatedField(phi_flux);

  this->addDependentField(grad_phi);
  this->addDependentField(rel_perm);

  std::string name = "PotentialFlux";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
PotentialFlux<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT& gradphi = grad_phi(cell,ip,dim);
        const ScalarT& relperm = rel_perm(cell,ip);
        phi_flux(cell,ip,dim) = Lambda2*relperm*gradphi;
      }
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
PotentialFlux<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Flux Name", "?");
  p->set<std::string>("Gradient Name", "?");
  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);
  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

