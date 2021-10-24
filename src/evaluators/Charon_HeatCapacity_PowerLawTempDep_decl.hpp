
#ifndef CHARON_HEATCAPACITY_POWERLAWTEMPDEP_DECL_HPP
#define CHARON_HEATCAPACITY_POWERLAWTEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief The heat capacity is modeled as
 * heat_cap = \rho * (c300 + c1 * num/den) in unit of [J/(K.cm^3)],
 * num = (TL / 300 K)^\beta - 1,
 * den = (TL / 300 K)^\beta + c1 / c300.
 * Here \rho is the mass density in [g/cm^3], TL is the lattice temperature in [K],
 * c300 in [J/(K.g)], c1 in [J/(K.g)], and \beta are fitting parameters,
 * which can be changed by user in the input file or come from
 * charon::Material_Properties when not specified/changed.
 *
 *
 *
 * Specification of the heat capacity model in the input file takes the form of
 * <ParameterList name="Heat Capacity">
 *   <Parameter name="Value" type="string" value="PowerLawTempDep" />
 *   <Parameter name="Mass Density"       type="double" value="2.33" />
 *   <Parameter name="Heat Capacity c300" type="double" value="0.711" />
 *   <Parameter name="Heat Capacity c1"   type="double" value="0.255" />
 *   <Parameter name="Heat Capacity beta" type="double" value="1.85" />
 * </ParameterList>
 */


//! obtain the lattice heat capacity using the PowerLawTempDep model
template<typename EvalT, typename Traits>
class HeatCapacity_PowerLawTempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    HeatCapacity_PowerLawTempDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> heat_cap;  // [scaled]

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;  // lattice temperature [scaled]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;  // temperature scaling, [K]
  double cL0; // heat capacity scaling [J/(K.cm^3)]

  int num_points;

  // Heat capacity model parameters
  double rho, c300, c1, beta;

  // initialize heat capacity parameters
  void initialize(const std::string& matName, const Teuchos::ParameterList& plist);

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class HeatCapacity_PowerLawTempDep


}

#endif
