
#ifndef CHARON_PREVPOTENTIALGRAD_IMPL_HPP
#define CHARON_PREVPOTENTIALGRAD_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"
//#include "Sacado_No_Kokkos.hpp"


/*
Compute the previous grad of potential at IPs for a previous time: 

*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
PrevPotentialGrad<EvalT, Traits>::
PrevPotentialGrad(
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
  time = 0.0;

  // Current density name
  string currentName = p.get<string>("Current Name");

  // Dependent fields
  grad_phi = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.phi,vector);
  this->addDependentField(grad_phi);
  
  // Evaluated field
  prev_grad_phi = MDField<ScalarT,Cell,IP,Dim>(currentName,vector);
  this->addEvaluatedField(prev_grad_phi);

  std::string name = "PrevPotentialGrad";
  this->setName(name);

}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
PrevPotentialGrad<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // get the current (present) time in [s]
  double curr_time = workset.time * t0;
  double timestep = t0*1.0/workset.alpha;

  // save gradient of potential for subsequent computations at the 
  // end of a convergent time step
  if ( curr_time >= time + timestep) {
    // save gradient of potential 
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      {
	for (std::size_t ip = 0; ip < num_ip; ++ip)
	  {
	    for (std::size_t dim = 0; dim < num_dim; ++dim)
	      {
		prev_grad_phi(cell,ip,dim) = grad_phi(cell,ip,dim);
	      }
	  }
      }
    // update time
    time = curr_time; 
  }
  
}
 

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
PrevPotentialGrad<EvalT, Traits>::getValidParameters() const
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

