
#ifndef CHARON_DISPLACEMENTCURRENTDENSITY_IMPL_HPP
#define CHARON_DISPLACEMENTCURRENTDENSITY_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"


/*
Compute the displacement current density at IPs: 

*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DisplacementCurrentDensity<EvalT, Traits>::
DisplacementCurrentDensity(
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

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  E0 = scaleParams->scale_params.E0;
  J0 = scaleParams->scale_params.J0;

  // Current density name
  string currentName = p.get<string>("Current Name");

  // Dependent fields
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.phi,vector);
  grad_phi_prev = MDField<const ScalarT,Cell,IP,Dim>(n.field.grad_phi_prev,vector);
  rel_perm = MDField<const ScalarT,Cell,IP>(n.field.rel_perm,scalar);

  // Dependent fields
  this->addDependentField(grad_phi);
  this->addDependentField(grad_phi_prev);
  this->addDependentField(rel_perm);

  // Evaluated field
  current_density = MDField<ScalarT,Cell,IP,Dim>(currentName,vector);
  this->addEvaluatedField(current_density);

  std::string name = "DisplacementCurrentDensity";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DisplacementCurrentDensity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();

  ScalarT timestep = t0 * 1/workset.alpha; // [s]
  
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT gradphi = E0 * grad_phi(cell,ip,dim); // [V/cm]
	const ScalarT gradphi_prev = E0 * grad_phi_prev(cell,ip,dim); // [V/cm]
	const ScalarT eps = rel_perm(cell,ip) * phyConst.eps0; // [C/(Vcm)]

	current_density(cell,ip,dim) = -(eps * (gradphi - gradphi_prev)/timestep)/J0;
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
DisplacementCurrentDensity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  
  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);
 
  p->set("Current Name", "?");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;

}

}

#endif

