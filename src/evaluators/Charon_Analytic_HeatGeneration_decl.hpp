
#ifndef CHARON_ANALYTIC_HEATGENERATION_DECL_HPP
#define CHARON_ANALYTIC_HEATGENERATION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! obtain the heat generation
template<typename EvalT, typename Traits>
class Analytic_HeatGeneration
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Analytic_HeatGeneration(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> heat_gen;  // [scaled]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double H0; // heat generation scaling [W/cm^3]
  double T0; // temperature scaling [K]

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;  // [scaled]

  int num_points;

  std::string heatGenType;

  double constValue;
  double linFactor;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Analytic_HeatGeneration


}

#endif
