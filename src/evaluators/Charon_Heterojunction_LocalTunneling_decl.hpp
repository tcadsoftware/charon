
#ifndef CHARON_HETEROJUNCTION_LOCALTUNNELING_DECL_HPP
#define CHARON_HETEROJUNCTION_LOCALTUNNELING_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Panzer_Evaluator_Macros.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

PANZER_EVALUATOR_CLASS(Heterojunction_LocalTunneling)

private:

  ScalarT evaluateIntegration(const ScalarT & intUpLimit, const ScalarT & fieldParam);

  // output
  PHX::MDField<ScalarT,Cell,Point> flux_local_tunnel;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;    // scaled
  PHX::MDField<const ScalarT,Cell,Point> normal_dot_grad;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling [K]
  double E0; // electric field scaling [V/cm]

  int num_points;
  int num_dims;
  int detailIndex;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  double bandOffset;  // [eV]
  double tunnelMass;  // in units of [m0]

  std::string dof_name;
  std::string flux_name;
  std::string other_dof_name;

  bool isPrimarySideLeft;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

PANZER_EVALUATOR_CLASS_END

}

#endif
