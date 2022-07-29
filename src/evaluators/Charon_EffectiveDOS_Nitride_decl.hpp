#ifndef CHARON_EFFECTIVEDOS_NITRIDE_DECL_HPP
#define CHARON_EFFECTIVEDOS_NITRIDE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {
/**
 *  \brief Computes Nitride Density of States as a function of mole fraction
 *         Vanheusden et al. Nature Vol386 (1997):587-589
 */

PHX_EVALUATOR_CLASS(EffectiveDOS_Nitride)

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> elec_effdos; 
  PHX::MDField<ScalarT,Cell,Point> hole_effdos; 

  // input
  PHX::MDField<const ScalarT,Cell,Point> latttemp;
  PHX::MDField<const ScalarT,Cell,Point> molefrac;
  
  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;
  double C0;

  int num_points;

  std::string materialName;
  
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
PHX_EVALUATOR_CLASS_END

}

#endif
