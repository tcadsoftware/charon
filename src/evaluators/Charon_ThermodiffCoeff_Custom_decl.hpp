
#ifndef CHARON_THERMODIFFCOEFF_CUSTOM_DECL_HPP
#define CHARON_THERMODIFFCOEFF_CUSTOM_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain thermodiffusion coefficient for ions
template<typename EvalT, typename Traits>
class ThermodiffCoeff_Custom
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ThermodiffCoeff_Custom(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // initialize mobility parameters
  void initialize(const Teuchos::ParameterList& plist);

  // output
  PHX::MDField<ScalarT,Cell,Point> thermodiff_coeff; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;   // scaled

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling, [K]
  double D0;

  int num_points;
  int num_edges;

  bool isEdgedl;

  // input parameters
  double multiplier, sign, slope;
  double minTemp, maxTemp, minActE, maxActE;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class ThermodiffCoeff_Custom


}

#endif
