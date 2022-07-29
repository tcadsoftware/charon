
#ifndef CHARON_QUANTUMPOTENTIALFLUX_IMPL_HPP
#define CHARON_QUANTUMPOTENTIALFLUX_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Physical_Constants.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
QuantumPotentialFlux<EvalT, Traits>::
QuantumPotentialFlux(
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

  qp_flux  = MDField<ScalarT,Cell,IP,Dim>(p.get<string>("Flux Name"),vector);
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(p.get<string>("Gradient Phi"),vector);
  grad_qp  = MDField<const ScalarT,Cell,IP,Dim>(p.get<string>("Gradient QP"),vector);
  latt_temp = MDField<const ScalarT,Cell,IP>(p.get<string>("Lattice Temperature"), scalar);


  fitParam = p.get<double>("Fit Parameter");

  this->addEvaluatedField(qp_flux);

  this->addDependentField(grad_phi);
  this->addDependentField(grad_qp);
  this->addDependentField(latt_temp);

  std::string name = "QuantumPotentialFlux";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
QuantumPotentialFlux<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{

  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      const ScalarT& lattTemp  = latt_temp (cell,ip);
      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT& gradphi = grad_phi(cell,ip,dim);
        const ScalarT& gradqp  = grad_qp (cell,ip,dim);

        qp_flux(cell,ip,dim) = (gradphi+gradqp*fitParam)/lattTemp;
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
QuantumPotentialFlux<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Flux Name", "?");
  p->set<std::string>("Gradient Phi", "?");
  p->set<std::string>("Gradient QP", "?");
  p->set<std::string>("Lattice Temperature", "?");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<double>("Fit Parameter", 1.0);

  return p;
}

}

#endif

