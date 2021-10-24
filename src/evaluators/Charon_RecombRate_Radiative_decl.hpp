
#ifndef CHARON_RECOMBRATE_RADIATIVE_DECL_HPP
#define CHARON_RECOMBRATE_RADIATIVE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain radiative recombination rate
template<typename EvalT, typename Traits>
class RecombRate_Radiative
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RecombRate_Radiative(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> rad_rate; // @ IPs
  PHX::MDField<ScalarT,Cell,Point> rad_deriv_e;
  PHX::MDField<ScalarT,Cell,Point> rad_deriv_h;

  // input
  PHX::MDField<const ScalarT,Cell,Point> intrin_conc;
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;

  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // [cm^-3]
  double R0;
  double T0; // [K]


  int num_points;
  double recomb_coeff;
  bool bUseFD;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class RecombRate_Radiative


}

#endif
