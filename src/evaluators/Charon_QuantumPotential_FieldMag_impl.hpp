
#ifndef CHARON_QUANTUMPOTENTIAL_FIELDMAG_IMPL_HPP
#define CHARON_QUANTUMPOTENTIAL_FIELDMAG_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
QuantumPotentialFieldMag<EvalT, Traits>::
QuantumPotentialFieldMag(
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

  qp_fieldmag  = MDField<ScalarT,Cell,IP>(p.get<string>("Field Name"), scalar);
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(p.get<string>("Gradient Phi"),vector);
  grad_qp  = MDField<const ScalarT,Cell,IP,Dim>(p.get<string>("Gradient QP"),vector);
  latt_temp = MDField<const ScalarT,Cell,IP>(p.get<string>("Lattice Temperature"), scalar);

  fitParam = p.get<double>("Fit Parameter");

  this->addEvaluatedField(qp_fieldmag);

  this->addDependentField(grad_phi);
  this->addDependentField(grad_qp);
  this->addDependentField(latt_temp);

  std::string name = "QuantumPotentialFieldMag";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
QuantumPotentialFieldMag<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{

  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      qp_fieldmag(cell,ip) = 0;
      const ScalarT& lattTempScal = pow(latt_temp (cell,ip), 2);
      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT& gradphi = grad_phi(cell,ip,dim);
        const ScalarT& gradqp  = grad_qp (cell,ip,dim);

        qp_fieldmag(cell,ip) += pow((gradphi+gradqp*fitParam),2)/lattTempScal;
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
QuantumPotentialFieldMag<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Field Name", "?");
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

