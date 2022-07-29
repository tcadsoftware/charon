#ifndef CHARON_PERMITTIVITY_NITRIDE_DECL_HPP
#define CHARON_PERMITTIVITY_NITRIDE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {
/**
 *  \brief Computes permittivity of nitrides as a function of mole fraction
 *         O.Ambacher et al. Appl.Phys. Vol87 No1 2000
 */

PHX_EVALUATOR_CLASS(Permittivity_Nitride)

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> rel_perm; 

  // input
  PHX::MDField<const ScalarT,Cell,Point> molefrac;

  int num_points;

  std::string materialName;
  
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
PHX_EVALUATOR_CLASS_END

}

#endif
