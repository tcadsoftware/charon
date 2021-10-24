
#ifndef CHARON_SUPG_TAU_TANH_DECL_HPP
#define CHARON_SUPG_TAU_TANH_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

template<typename EvalT, typename Traits>
class SUPG_Tau_Tanh
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SUPG_Tau_Tanh(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;

  // output
  PHX::MDField<ScalarT,Cell,Point> tau;

  // input
  PHX::MDField<const ScalarT,Cell,Point> peclet; // carrier peclet number
  PHX::MDField<const ScalarT,Cell,Point,Dim> velocity; // carrier velocity
  PHX::MDField<const ScalarT,Cell,Point> recomb_deriv;

  int num_points;
  int num_dims;
  int ir_degree;
  int ir_index;  // integration rule index to get gc

  std::string carrType;
  bool add_source_stab;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SUPG_Tau_Tanh


}

#endif
