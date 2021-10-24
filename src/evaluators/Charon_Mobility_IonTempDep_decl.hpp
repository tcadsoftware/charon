
#ifndef CHARON_MOBILITY_IONTEMPDEP_DECL_HPP
#define CHARON_MOBILITY_IONTEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! electron mobility model
template<typename EvalT, typename Traits>
class Mobility_IonTempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_IonTempDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // initialize mobility parameters
  void initMobilityParams(const Teuchos::ParameterList& mobParamList);

  // output
  PHX::MDField<ScalarT,Cell,Point> mobility;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp; //lattice temperature [scaled]
  PHX::MDField<const ScalarT,Cell,Point> ion_density;
  PHX::MDField<const ScalarT,Cell,Point> elec_density;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
  double T0;  // temperature scaling, [K]

  int num_points;

  // mobility model parameters
  double maxSigma0, minSigma0, maxSigma;
  double minIonDens, maxIonDens, medIonDens;
  double maxActE, minActE;
  double slopeSigma0, slopeActE;
  double maxMobValue;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_IonTempDep


}

#endif
