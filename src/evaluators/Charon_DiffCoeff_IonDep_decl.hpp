
#ifndef CHARON_DIFFCOEFF_IONDEP_DECL_HPP
#define CHARON_DIFFCOEFF_IONDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain diffusion coefficient for ions
template<typename EvalT, typename Traits>
class DiffCoeff_IonDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DiffCoeff_IonDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> diffcoeff; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,Point> mobility;  // scaled
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> carr_dens;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;  // temperature scaling, [K]
  double C0;  // conc. scaling, [cm^-3]

  int num_points;

  double maxIonDens;
  double maxFactor;

  std::string funcType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DiffCoeff_IonDep


}

#endif
