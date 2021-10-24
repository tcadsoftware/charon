
#ifndef CHARON_BANDGAP_TEMPDEP_DECL_HPP
#define CHARON_BANDGAP_TEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! calculate temperature-dependent band gap
template<typename EvalT, typename Traits>
class BandGap_TempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BandGap_TempDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> band_gap;
  PHX::MDField<ScalarT,Cell,Point> affinity;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;

  int num_points;

  // bandgap model parameters
  double Eg300, alpha, beta, Chi300;

  bool isCalcAffinity;


private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BandGap_TempDep


}

#endif
