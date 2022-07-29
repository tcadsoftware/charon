
#ifndef CHARON_DISPLACEMENTCURRENTONCONTACT_IMPL_HPP
#define CHARON_DISPLACEMENTCURRENTONCONTACT_IMPL_HPP

#include "Kokkos_ViewFactory.hpp"


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
DisplacementCurrentOnContact<EvalT, Traits>::
DisplacementCurrentOnContact(
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
  rel_perm = MDField<const ScalarT,Cell,IP>(n.field.rel_perm,scalar);

  // Dependent fields
  this->addDependentField(grad_phi);
  this->addDependentField(rel_perm);

  // Evaluated field
  current_density = MDField<ScalarT,Cell,IP,Dim>(currentName,vector);

  this->addEvaluatedField(current_density);

  std::string name = "DisplacementCurrentOnContact";
  this->setName(name);

  prev_time = 0.0;
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void DisplacementCurrentOnContact<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  grad_phi_prev = Kokkos::createDynRankView(grad_phi.get_static_view(),
		    "grad_phi_prev",grad_phi.dimension(0), num_ip, num_dim);

}



///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DisplacementCurrentOnContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();

  double curr_time = t0*workset.time; // [s]
  double timestep = curr_time - prev_time; // [s]

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t ip = 0; ip < num_ip; ++ip)
    {
      for (std::size_t dim = 0; dim < num_dim; ++dim)
      {
        const ScalarT gradphi = E0 * grad_phi(cell,ip,dim); // [V/cm]
	const ScalarT gradphi_prev = E0 * grad_phi_prev(cell,ip,dim); // [V/cm]
	const ScalarT eps = rel_perm(cell,ip) * phyConst.eps0; // [C/(Vcm)]

        // avoid problems with timestep being equal to zero during initial step
        if (timestep > std::numeric_limits<double>::min())
          current_density(cell,ip,dim) = -(eps * (gradphi - gradphi_prev))/timestep/J0;
        else
          current_density(cell,ip,dim) = 0.0;
	
	// update previous time and previous grad potential 
	grad_phi_prev(cell,ip,dim) = grad_phi(cell,ip,dim); // scaled
      }
    }
  }

  prev_time = curr_time; // [s]

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DisplacementCurrentOnContact<EvalT, Traits>::getValidParameters() const
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

