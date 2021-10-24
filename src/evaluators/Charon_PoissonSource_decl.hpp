
#ifndef CHARON_POISSONSOURCE_DECL_HPP
#define CHARON_POISSONSOURCE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

// Evaluates the Poisson source term: p-n+Nd-Na, at any given data layout
namespace charon {

template<typename EvalT, typename Traits>
class PoissonSource
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    PoissonSource(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,Point> poissonSource; // scaled

  // input
  //PHX::MDField<const ScalarT,Cell,Point> intrin_fermi;  // [eV]
  PHX::MDField<const ScalarT,Cell,Point> doping;     // scaled (net doping)
  //PHX::MDField<const ScalarT,Cell,Point> intrin_conc; // scaled (effective intrin. conc.)
  PHX::MDField<const ScalarT,Cell,Point> edensity;   // scaled
  PHX::MDField<const ScalarT,Cell,Point> hdensity;   // scaled
  PHX::MDField<const ScalarT,Cell,Point> iondensity;   // scaled
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;   // scaled

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;

  int num_points;

  std::string solveElectron;
  std::string solveHole;

  bool solveIon;
  int ionCharge;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class PoissonSource


}

#endif
