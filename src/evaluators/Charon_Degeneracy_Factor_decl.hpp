
#ifndef CHARON_DEGENERACY_FACTOR_DECL_HPP
#define CHARON_DEGENERACY_FACTOR_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_FermiDirac_Integral.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

template<typename EvalT, typename Traits>
class Degeneracy_Factor
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Degeneracy_Factor(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> elec_degfactor;
  PHX::MDField<ScalarT,Cell,Point> hole_degfactor;

  // input
  PHX::MDField<const ScalarT,Cell,Point> elec_density;
  PHX::MDField<const ScalarT,Cell,Point> hole_density;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;

  int num_points;

  bool bUseFD;  // bUseFD=true, use Fermi-Dirac statistics, otherwise, use Boltzmann

  std::string fdFormula;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;

}; // end of class Degeneracy_Factor


}

#endif
