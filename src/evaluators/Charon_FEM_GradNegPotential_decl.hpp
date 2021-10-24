
#ifndef CHARON_FEM_GRADNEGPOTENTIAL_DECL_HPP
#define CHARON_FEM_GRADNEGPOTENTIAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

// Assign the x and y components of the negative potential gradient to scalar fields
// that can be output to Exodus

namespace charon {

template<typename EvalT, typename Traits>
class FEM_GradNegPotential
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    FEM_GradNegPotential(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> grad_negpot_x;
  PHX::MDField<ScalarT,Cell,Point> grad_negpot_y;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_potential;  // gradient of potential

  // IP
  int num_points;
  int num_dims;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class FEM_GradNegPotential


}

#endif
