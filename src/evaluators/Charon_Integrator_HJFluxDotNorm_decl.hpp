
#ifndef CVFEM_INTEGRATOR_HJFLUXDOTNORM_DECL_HPP
#define CVFEM_INTEGRATOR_HJFLUXDOTNORM_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Intrepid2_Basis.hpp"

#include "Panzer_Evaluator_Macros.hpp"

using panzer::Cell;
using panzer::Dim;
using panzer::BASIS;
using panzer::IP;
using panzer::Edge;

namespace charon {

// Integrate over a subcontrol volume surface \int_cv_surf Jhj \cdot norm dS
PANZER_EVALUATOR_CLASS(Integrator_HJFluxDotNorm)

  // output
  PHX::MDField<ScalarT,Cell,BASIS> residual;

  // input
  PHX::MDField<const ScalarT,Cell,IP> hj_flux;

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_multipliers;

  int num_nodes;
  int num_ips;
  int num_dims;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  std::string residual_name;
  std::string flux_name;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;


PANZER_EVALUATOR_CLASS_END

}


#endif
