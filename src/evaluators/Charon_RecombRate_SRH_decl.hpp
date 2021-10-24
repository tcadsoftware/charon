
#ifndef CHARON_RECOMBRATE_SRH_DECL_HPP
#define CHARON_RECOMBRATE_SRH_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

  // evaluates the n0*p0 product for Fermi-Dirac statistics
  // can be called by all relevant recomb. evaluators that need n0*p0
  // all ScalarT arguments are in physical units

  template<typename EvalT, typename Traits>
  class FermiDiracIntrinsicDensity {
  public:
    static
    typename EvalT::ScalarT evaluateFDIntrinsicDensity(
         const typename EvalT::ScalarT& n,
         const typename EvalT::ScalarT& p,
         const typename EvalT::ScalarT& niMB,
         const typename EvalT::ScalarT& Nc,
         const typename EvalT::ScalarT& Nv,
         const typename EvalT::ScalarT& Eg,
         const typename EvalT::ScalarT& kbT,
         const Teuchos::RCP<charon::FermiDiracIntegral<EvalT> >& invFDInt);
  };


//! obtain SRH recombination rate
template<typename EvalT, typename Traits>
class RecombRate_SRH
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RecombRate_SRH(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> srh_rate;
  PHX::MDField<ScalarT,Cell,Point> srh_deriv_e;
  PHX::MDField<ScalarT,Cell,Point> srh_deriv_h;

  // input
  PHX::MDField<const ScalarT,Cell,Point> elifetime;
  PHX::MDField<const ScalarT,Cell,Point> hlifetime;
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
  double T0; // [K]

  int num_points;
  bool bUseFD;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class RecombRate_SRH


}

#endif
