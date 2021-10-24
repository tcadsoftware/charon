
#ifndef CHARON_THERMALCONDUCT_TEMPDEP_DECL_HPP
#define CHARON_THERMALCONDUCT_TEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief The thermal conductivity is modeled as kappa = 1 / (a + b*T + c*T^2),
 * where T is the lattice temperature, a, b, and c parameters can be changed
 * by user in the input file or come from charon::Material_Properties when not
 * specified/changed.
 * 
 *
 * Specification of the thermal conductivity model in the input file takes the form of
 * <ParameterList name="Thermal Conductivity">
 *   <Parameter name="Value" type="string" value="TempDep" />
 *   <Parameter name="a" type="double" value="0.03" />
 *   <Parameter name="b" type="double" value="1.56e-3" />
 *   <Parameter name="c" type="double" value="1.65e-6" />
 * </ParameterList>
 */


//! obtain the lattice thermal conductivity using the TempDep model
template<typename EvalT, typename Traits>
class ThermalConduct_TempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ThermalConduct_TempDep(
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
  double a, b, c;

  // initialize thermal conductivity parameters
  void initialize(const std::string& matName, const Teuchos::ParameterList& plist);

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class ThermalConduct_TempDep


}

#endif
