
#ifndef CHARON_SORETCOEFF_TEMPDEP_DECL_HPP
#define CHARON_SORETCOEFF_TEMPDEP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain the soret coefficient for ion/vacancy
template<typename EvalT, typename Traits>
class SoretCoeff_TempDep
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SoretCoeff_TempDep(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> soret_coeff;  // [scaled]

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;  // lattice temperature [scaled]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;  // temperature scaling, [K]

  int num_points;
  int num_edges;

  bool isEdgedl;

  // Soret coefficient model parameters
  double Ua;
  double sign;  // sign of the Soret coefficient

  // initialize the Soret coefficient parameters
  void initialize(const std::string& matName, const Teuchos::ParameterList& plist);

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SoretCoeff_TempDep


}

#endif
