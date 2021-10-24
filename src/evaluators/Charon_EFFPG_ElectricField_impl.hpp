
#ifndef CHARON_EFFPG_ELECTRICFIELD_IMPL_HPP
#define CHARON_EFFPG_ELECTRICFIELD_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"


/*
Compute effective electric field at IPs inversely from the current density
at IPs that has been calculated using the EFFPG method, i.e.,
Fneff = (Jn - Dn*grad(n))/(n*mun), Fpeff = (Jp + Dp*grad(p))/(p*mup).

Also compute gradient of quasi-fermi potential at IPs, i.e.,
e_grad_qfp = -grad(n)/n - Fneff, h_grad_qfp = grad(p)/p - Fpeff ; or
e_grad_qfp = -Jn / (n*mun), h_grad_qfp = -Jp / (p*mup).

The two expressions are equivalent provided the Einstein relation Dn/mun (Dp/mup)
= kb*T/q holds, which is true in the general continuity equations by D. Schroeder,
T. Ostermann and O. Kalz, "Comparison of transport
models for the simulation of degenerate semiconductors," Semicond. Sci. Technol.9
(1994) 364-369.

The Fn\peff and e(h)_grad_qfp are needed in computing avalanche generation rate.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
EFFPG_ElectricField<EvalT, Traits>::
EFFPG_ElectricField(
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
  num_dim = vector->dimension(2);

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    efield = MDField<ScalarT,Cell,Point,Dim>(n.field.elec_efield,vector);
    grad_qfp = MDField<ScalarT,Cell,Point,Dim>(n.field.elec_grad_qfp,vector);
    curr_density = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_curr_density,vector);
    grad_density = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.edensity,vector);
    density = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,Point>(n.field.elec_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,Point>(n.field.elec_mobility,scalar);
    sign = -1.0;
  }
  else if (carrType == "Hole")
  {
    efield = MDField<ScalarT,Cell,Point,Dim>(n.field.hole_efield,vector);
    grad_qfp = MDField<ScalarT,Cell,Point,Dim>(n.field.hole_grad_qfp,vector);
    curr_density = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_curr_density,vector);
    grad_density = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.hdensity,vector);
    density = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,Point>(n.field.hole_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,Point>(n.field.hole_mobility,scalar);
    sign = 1.0;
  }

  // Evaluated field
  this->addEvaluatedField(efield);
  this->addEvaluatedField(grad_qfp);

  // Dependent fields
  this->addDependentField(curr_density);
  this->addDependentField(grad_density);
  this->addDependentField(density);
  this->addDependentField(diff_coeff);
  this->addDependentField(mobility);

  std::string name = "EFFPG_ElectricField";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
EFFPG_ElectricField<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // All quantities used here are scaled !

  // Compute the effective electric field at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t point = 0; point < num_points; ++point)
    {
      const ScalarT& diff = diff_coeff(cell,point);
      const ScalarT& mob = mobility(cell,point);
      const ScalarT& dens = density(cell,point);

      for (std::size_t dim = 0; dim < num_dim; ++ dim)
      {
        const ScalarT& curr_dens = curr_density(cell,point,dim);
        const ScalarT& grad_dens = grad_density(cell,point,dim);

        ScalarT Feff = (curr_dens + sign*diff*grad_dens) / (dens*mob);
        efield(cell,point,dim) = Feff;
        grad_qfp(cell,point,dim) = sign * grad_dens/dens - Feff;
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
EFFPG_ElectricField<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;

}

}

#endif

