#ifndef CHARON_ANALYTIC_COMPARISON_HPP
#define CHARON_ANALYTIC_COMPARISON_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Dimension.hpp"

#include "Kokkos_DynRankView.hpp"

namespace charon {

  using panzer::Cell;
  using panzer::Point;

/**
 * \brief Analiytic comparison of field values.
 *
 * This function compares the analytic and computed solution of a field
 * variable. It uses a straight difference of the two values.
 */
template<typename EvalT, typename Traits>
class AnalyticComparison
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    AnalyticComparison(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<const ScalarT,Cell,Point> simulation_value;
  PHX::MDField<const ScalarT,Cell,Point> analytic_value;
  PHX::MDField<ScalarT,Cell,Point> error;

}; // end of class AnalyticComparison


/**
 * \brief Analiytic comparison of field values.
 *
 * This function compares the analytic and computed solution of a field
 * variable. It uses an L2 calculation for the result.
 */
template<typename EvalT, typename Traits>
class AnalyticComparison_L2Error
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    AnalyticComparison_L2Error(
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

private:
  void postprocess(std::ostream& os);
  void preEvaluate(typename Traits::PreEvalData d);
  void postEvaluate(typename Traits::PostEvalData d);

  PHX::MDField<const ScalarT,Cell,Point> simulation_value;
  PHX::MDField<const ScalarT,Cell,Point> analytic_value;
  PHX::MDField<ScalarT,Cell,Point> error;

  ScalarT l2error;

  Kokkos::DynRankView<ScalarT,PHX::Device> integral;

  int quadOrder, quadIndex;
  int numIPs, numNodes, basisIndex;
  std::string basisName;

  Teuchos::RCP<Teuchos::Comm<int> const> comm;

}; // end of class AnalyticComparison_L2Error


/**
 * \brief Analiytic comparison of field values.
 *
 * This function compares the analytic and computed solution of a field
 * variable. It uses a relative difference of the form
 * \f$\text{err}=\left| \frac{(v-a)}{a} \right| \f$. Where \f$v\f$ is
 * the computed solution and \f$a\f$ is the analytic solution. The user
 * can choose to disable the absolute value fromt he input
 * file. Addtionally the formula breaks down when the values are very
 * close to zero. In that case a simple absolute difference is
 * performed.
 */
template<typename EvalT, typename Traits>
class AnalyticComparison_RelError
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    AnalyticComparison_RelError(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

  PHX::MDField<const ScalarT,Cell,Point> simulation_value;
  PHX::MDField<const ScalarT,Cell,Point> analytic_value;
  PHX::MDField<ScalarT,Cell,Point> error;

  bool useAbs;

}; // end of class AnalyticComparison_RelError


}

#endif
