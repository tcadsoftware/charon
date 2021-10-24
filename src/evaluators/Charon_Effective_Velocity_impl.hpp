
#ifndef CHARON_EFFECTIVE_VELOCITY_IMPL_HPP
#define CHARON_EFFECTIVE_VELOCITY_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"

/*
To include the temperature gradient in the calculation of effective velocity:
(the effective velocity is used in the SUPG stabilization)

veff,n = vn - mun*grad_temp = -mun*Fn - mun*grad_temp for electrons,
veff,p = vp - mup*grad_temp = mup*Fp - mup*grad_temp for holes,
veff,i = vi - mui*T*Si*grad_temp = vi - DTi*grad_temp for ions or vacancies.
mui*T = Di when the Einstein relation holds.

All quantities are in scaled units. See Charon Notebook IV pg.197-198 for details.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Effective_Velocity<EvalT, Traits>::
Effective_Velocity(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  carrType = p.get<string>("Carrier Type");
  includeSoret = p.get<bool>("Include Soret Effect");

  // Carrier dependent fields
  if (carrType == "Electron")
  {
    eff_velocity = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_eff_velocity,vector);
    velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.elec_velocity,vector);
    grad_temp = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.latt_temp,vector);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.elec_mobility,scalar);

    this->addDependentField(velocity);
    this->addDependentField(grad_temp);
    this->addDependentField(mobility);
  }
  else if (carrType == "Hole")
  {
    eff_velocity = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_eff_velocity,vector);
    velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.hole_velocity,vector);
    grad_temp = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.latt_temp,vector);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.hole_mobility,scalar);

    this->addDependentField(velocity);
    this->addDependentField(grad_temp);
    this->addDependentField(mobility);
  }
  else if (carrType == "Ion")
  {
    eff_velocity = MDField<ScalarT,Cell,IP,Dim>(n.field.ion_eff_velocity,vector);
    velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.ion_velocity,vector);
    grad_temp = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.latt_temp,vector);
    thermodiff_coeff = MDField<const ScalarT,Cell,IP>(n.field.ion_thermodiff_coeff,scalar);

    this->addDependentField(velocity);
    this->addDependentField(grad_temp);
    this->addDependentField(thermodiff_coeff);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole or Ion !");

  // Evaluated fields
  this->addEvaluatedField(eff_velocity);

  std::string name = "Effective_Carrier_Velocity";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Effective_Velocity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

 if (includeSoret)  // include the soret effect
 {
  if ( (carrType == "Electron") || (carrType == "Hole") )
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t ip = 0; ip < num_ip; ++ip)
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          eff_velocity(cell,ip,dim) = velocity(cell,ip,dim) - mobility(cell,ip) * grad_temp(cell,ip,dim);
  }

  else if (carrType == "Ion")
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t ip = 0; ip < num_ip; ++ip)
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          eff_velocity(cell,ip,dim) = velocity(cell,ip,dim) - thermodiff_coeff(cell,ip) * grad_temp(cell,ip,dim);
  }
 }

 else  // set eff_velocity to the original velocity
 {
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
     for (std::size_t ip = 0; ip < num_ip; ++ip)
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          eff_velocity(cell,ip,dim) = velocity(cell,ip,dim);
 }

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Effective_Velocity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  p->set<int>("Ion Charge", 1, "The integer number of ion charge");

  p->set<bool>("Include Soret Effect", true);

  return p;
}

}

#endif
