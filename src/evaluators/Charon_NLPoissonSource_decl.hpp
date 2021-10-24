
#ifndef CHARON_NLPOISSONSOURCE_DECL_HPP
#define CHARON_NLPOISSONSOURCE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

// Evaluate the equilibrium nonlinear Poisson source:
// nie*exp(Ei/kbT)-nie*exp(-Ei/kbT)+dop

namespace charon {

template<typename EvalT, typename Traits>
class NLPoissonSource
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    NLPoissonSource(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,Point> nlpsrc; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,Point> intrin_fermi; // [eV]
  PHX::MDField<const ScalarT,Cell,Point> doping;       // scaled (net doping)
  PHX::MDField<const ScalarT,Cell,Point> intrin_conc;  // scaled (effective ni)
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;    // scaled

 // Scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;

  int num_points;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class NLPoissonSource


}

#endif
