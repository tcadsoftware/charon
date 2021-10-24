
#ifndef CHARON_THERMALCONDUCT_POWERLAWTEMPDEP_DECL_HPP
#define CHARON_THERMALCONDUCT_POWERLAWTEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief The thermal conductivity is modeled as kappa = kappa300 * (TL/300 K)^alpha,
 * where TL is the lattice temperature, kappa300 is the thermal conductivity at 300 K,
 * and alpha is the power law exponent. kappa300 and alpha parameters can be changed
 * by user in the input file, or come from charon::Material_Properties when not specified/changed.
 * kappa300 is in unit of [W/(K.cm)], TL is in [K], and alpha is unitless.
 *
 *
 * Specification of the model in the input file takes the form of
 * <ParameterList name="Thermal Conductivity">
 *   <Parameter name="Value" type="string" value="PowerLawTempDep" />
 *   <Parameter name="kappa300" type="double" value="1.48" />
 *   <Parameter name="alpha" type="double" value="-1.3" />
 * </ParameterList>
 */

//! obtain the lattice thermal conductivity using the PowerLawTempDep model
template<typename EvalT, typename Traits>
class ThermalConduct_PowerLawTempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ThermalConduct_PowerLawTempDep(
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
  double T0;  // temperature scaling [K]
  double kL0; // thermal conductivity scaling [W/(K.cm)]

  int num_points;

  // Thermal conductivity model parameters
  double kappa300, alpha;

  // initialize thermal conductivity parameters
  void initialize(const std::string& matName, const Teuchos::ParameterList& plist);

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class ThermalConduct_PowerLawTempDep


}

#endif
