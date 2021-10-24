
#ifndef CHARON_THERMALCONDUCT_LINEARTEMPDEP_DECL_HPP
#define CHARON_THERMALCONDUCT_LINEARTEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain the lattice thermal conductivity
template<typename EvalT, typename Traits>
class ThermalConduct_LinearTempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ThermalConduct_LinearTempDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> ther_cond;  // [scaled]

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;  // lattice temperature [scaled]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;   // temperature scaling [K]
  double kL0;  // thermal conductivity scaling [W/(K.cm)]

  int num_points;

  // Thermal conductivity model parameters
  double kappa0, lambda, Tref;

  // initialize thermal conductivity parameters
  void initialize(const Teuchos::ParameterList& plist);

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class ThermalConduct_LinearTempDep


}

#endif
