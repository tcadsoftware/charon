#ifndef CHARON_QUASIFERMILEVELS_DECL_HPP
#define CHARON_QUASIFERMILEVELS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

template<typename EvalT, typename Traits>
class QuasiFermiLevels
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{

public:

  QuasiFermiLevels(const Teuchos::ParameterList& p);

  void evaluateFields(typename Traits::EvalData d);

private:
  using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,Point> eQF; // [eV]
  PHX::MDField<ScalarT,Cell,Point> hQF; // [eV] 

  // input
  PHX::MDField<const ScalarT,Cell,Point> edens;
  PHX::MDField<const ScalarT,Cell,Point> hdens;
  PHX::MDField<const ScalarT,Cell,Point> lT;
  PHX::MDField<const ScalarT,Cell,Point> Ec;
  PHX::MDField<const ScalarT,Cell,Point> Ev;
  PHX::MDField<const ScalarT,Cell,Point> Nc;
  PHX::MDField<const ScalarT,Cell,Point> Nv;
  PHX::MDField<const ScalarT,Cell,Point> e_gamma;
  PHX::MDField<const ScalarT,Cell,Point> h_gamma;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling, [K]
  double C0; // concentration scaling, [cm-3]

  int num_points;
  
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
}; // end of class QuasiFermiLevels

}

#endif
