
#ifndef CHARON_INTEGRATOR_SUBCVNODESCALAR_HPP
#define CHARON_INTEGRATOR_SUBCVNODESCALAR_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Intrepid2_Basis.hpp"


using panzer::Cell;
using panzer::Dim;
using panzer::Point;


namespace charon {

// Integrate a scalar at BASIS over a subcontrol volume  \int_cv g dV

template<typename EvalT, typename Traits>
class Integrator_SubCVNodeScalar
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_SubCVNodeScalar(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  // output
  PHX::MDField<ScalarT,Cell,Point> residual;

  // input
  PHX::MDField<const ScalarT,Cell,Point> value;

  std::size_t basis_index;
  std::string basis_name;
  std::size_t int_rule_index;
  std::size_t int_rule_degree;
  int num_nodes;
  int num_ips;
  double multiplier;
  bool WithInterpolation;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;


}; // end of class Integrator_SubCVNodeScalar


}


#endif
