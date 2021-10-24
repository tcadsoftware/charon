
#ifndef CHARON_INTRINSICCONC_OLDSLOTBOOM_DECL_HPP
#define CHARON_INTRINSICCONC_OLDSLOTBOOM_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

template<typename EvalT, typename Traits>
class IntrinsicConc_OldSlotboom
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IntrinsicConc_OldSlotboom(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> nie;
  PHX::MDField<ScalarT,Cell,Point> effEg;
  PHX::MDField<ScalarT,Cell,Point> effChi;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  PHX::MDField<const ScalarT,Cell,Point> Eg;
  PHX::MDField<const ScalarT,Cell,Point> Chi;

  PHX::MDField<const ScalarT,Cell,Point> acceptor;
  PHX::MDField<const ScalarT,Cell,Point> donor;

  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]
  double T0; // temperature scaling, [K]

  int num_points;

  // material parameters
  double V0_BGN, N0_BGN, CON_BGN;

  bool includeBGN;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class IntrinsicConc_Slotboom


}

#endif
