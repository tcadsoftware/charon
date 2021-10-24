
#ifndef CHARON_FEM_CURRENTDENSITY_IMPL_HPP
#define CHARON_FEM_CURRENTDENSITY_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"


/*
Compute the FEM version of current density at IPs: Jn = n*mun*Fn + Dn*grad(n),
and Jp = p*mup*Fp - Dp*grad(p), where the field Fn (Fp) includes the BGN effect
when BGN is turned on.

*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
FEM_CurrentDensity<EvalT, Traits>::
FEM_CurrentDensity(
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
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Current density name
  string currentName = p.get<string>("Current Name");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    electric_field = MDField<const ScalarT,Cell,IP,Dim>(n.field.elec_efield,vector);
    grad_carr_dens = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.edensity,vector);
    carr_dens = MDField<const ScalarT,Cell,IP>(n.dof.edensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,IP>(n.field.elec_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.elec_mobility,scalar);
    sign = 1.0;
  }
  else if (carrType == "Hole")
  {
    electric_field = MDField<const ScalarT,Cell,IP,Dim>(n.field.hole_efield,vector);
    grad_carr_dens = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.hdensity,vector);
    carr_dens = MDField<const ScalarT,Cell,IP>(n.dof.hdensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,IP>(n.field.hole_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.hole_mobility,scalar);
    sign = -1.0;
  }

  // Evaluated field
  current_density = MDField<ScalarT,Cell,IP,Dim>(currentName,vector);
  this->addEvaluatedField(current_density);

  // Dependent fields
  this->addDependentField(electric_field);
  this->addDependentField(grad_carr_dens);
  this->addDependentField(carr_dens);
  this->addDependentField(diff_coeff);
  this->addDependentField(mobility);

  std::string name = "FEM_CurrentDensity";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_CurrentDensity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      const ScalarT& carr = carr_dens(cell,ip);
      const ScalarT& diff = diff_coeff(cell,ip);
      const ScalarT& mob = mobility(cell,ip);

      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT& gradcarr = grad_carr_dens(cell,ip,dim);
        const ScalarT& field = electric_field(cell,ip,dim);

        ScalarT current = carr*mob*field + sign*diff*gradcarr;
        current_density(cell,ip,dim) = current;

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
FEM_CurrentDensity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");
  p->set("Current Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;

}

}

#endif

