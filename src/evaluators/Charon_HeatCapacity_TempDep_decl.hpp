
#ifndef CHARON_HEATCAPACITY_TEMPDEP_DECL_HPP
#define CHARON_HEATCAPACITY_TEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief The heat capacity is modeled as heat_cap = a + b*T + c*T^2, where T is the
 * lattice temperature, a, b, and c parameters can be changed by user in the input
 * file or come from charon::Material_Properties when not specified/changed.
 *
 *
 * Specification of the heat capacity model in the input file takes the form of
 * <ParameterList name="Heat Capacity">
 *   <Parameter name="Value" type="string" value="TempDep" />
 *   <Parameter name="a" type="double" value="1.63" />
 *   <Parameter name="b" type="double" value="0" />
 *   <Parameter name="c" type="double" value="0" />
 * </ParameterList>
 */


//! obtain the lattice heat capacity using the TempDep model
template<typename EvalT, typename Traits>
class HeatCapacity_TempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    HeatCapacity_TempDep(
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
  double a, b, c;

  // initialize heat capacity parameters
  void initialize(const std::string& matName, const Teuchos::ParameterList& plist);

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class HeatCapacity_TempDep


}

#endif
