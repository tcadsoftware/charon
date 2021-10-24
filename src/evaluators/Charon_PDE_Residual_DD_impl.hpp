
#ifndef CHARON_PDE_RESIDUAL_DD_IMPL_HPP
#define CHARON_PDE_RESIDUAL_DD_IMPL_HPP

#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"


/*
Compute the PDE residual for SUPG stabilization, assume steady state and neglect
2nd-order terms, then
Rn = grad(n) \dot vn + G for electrons,
Rp = grad(p) \dot vp + G for holes,
Ri = grad(Ni) \dot vi for ions.

When DOF Time Derivative = true, add the transient terms \partial_n / \partial_t,
\partial_p / \partial_t, and \partial_Ni / \partial_t.

Quick tests show that the transient term does not add benefit to the stabilization,
hence, it is always disabled for now.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
PDE_Residual_DD<EvalT, Traits>::
PDE_Residual_DD(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  carrType = p.get<string>("Carrier Type");
  bDofTimeDeriv = false;
  if (p.isParameter("DOF Time Derivative"))
    bDofTimeDeriv = p.get<bool>("DOF Time Derivative");
  ionCharge = 1;

  bool includeSoret = false;
  if (p.isParameter("Include Soret Effect"))
    includeSoret = p.get<bool>("Include Soret Effect");

  // Carrier dependent fields
  if (carrType == "Electron")
  {
    residual = MDField<ScalarT,Cell,IP>(n.field.R_e,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.elec_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.elec_velocity,vector);

    grad_density = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.edensity, vector);
    total_recomb = MDField<const ScalarT,Cell,IP>(n.field.total_recomb,scalar);
    this->addDependentField(total_recomb);
    if (bDofTimeDeriv)
    {
      dof_time_deriv = MDField<const ScalarT,Cell,IP>(n.dxdt.edensity,scalar);
      this->addDependentField(dof_time_deriv);
    }
  }
  else if (carrType == "Hole")
  {
    residual = MDField<ScalarT,Cell,IP>(n.field.R_h,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.hole_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.hole_velocity,vector);

    grad_density = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.hdensity, vector);
    total_recomb = MDField<const ScalarT,Cell,IP>(n.field.total_recomb,scalar);
    this->addDependentField(total_recomb);
    if (bDofTimeDeriv)
    {
      dof_time_deriv = MDField<const ScalarT,Cell,IP>(n.dxdt.hdensity,scalar);
      this->addDependentField(dof_time_deriv);
    }
  }
  else if (carrType == "Ion")
  {
    residual = MDField<ScalarT,Cell,IP>(n.field.R_ion,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.ion_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.ion_velocity,vector);

    grad_density = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.iondensity, vector);
    if (bDofTimeDeriv)
    {
      dof_time_deriv = MDField<const ScalarT,Cell,IP>(n.dxdt.iondensity,scalar);
      this->addDependentField(dof_time_deriv);
    }
    ionCharge = p.get<int>("Ion Charge");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole or Ion !");

  // Evaluated field
  this->addEvaluatedField(residual);

  // Dependent fields
  this->addDependentField(velocity);
  this->addDependentField(grad_density);

  std::string name = "PDE_Residual_DD";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
PDE_Residual_DD<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

 if ( (carrType == "Electron") || (carrType == "Hole") )
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      ScalarT v_dot_gcarr = 0.0;

      // compute velocity \dot gradient(carrier)
      for (std::size_t dim = 0; dim < num_dim; ++dim)
        v_dot_gcarr += velocity(cell,ip,dim) * grad_density(cell,ip,dim);

      // pde residual
      residual(cell,ip) = v_dot_gcarr + total_recomb(cell,ip);

      // add the dof time derivative term
      if (bDofTimeDeriv)  residual(cell,ip) += dof_time_deriv(cell,ip);
    }
  }
 }  // end of if (carrType ...)

 else if (carrType == "Ion")
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      ScalarT v_dot_gcarr = 0.0;

      // compute velocity \dot gradient(carrier)
      for (std::size_t dim = 0; dim < num_dim; ++dim)
        v_dot_gcarr += velocity(cell,ip,dim) * grad_density(cell,ip,dim);

      // pde residual
      residual(cell,ip) = v_dot_gcarr;

      // add the dof time derivative term
      if (bDofTimeDeriv)  residual(cell,ip) += dof_time_deriv(cell,ip);

      // multiply the ion charge
      // residual(cell,ip) *= ionCharge;  // ion charge is removed from the ion continuity eqn on 11/12/2014
    }
  }
 }  // end of else if (carrType ...)

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
PDE_Residual_DD<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  p->set<bool>("DOF Time Derivative", false);
  p->set<int>("Ion Charge", 1);

  p->set<bool>("Include Soret Effect", false);

  return p;
}

}

#endif
