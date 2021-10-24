
#ifndef CHARON_GRADBASISDOTGRADDOF_DECL_HPP
#define CHARON_GRADBASISDOTGRADDOF_DECL_HPP

#include <string>
#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"

using panzer::Cell;
using panzer::BASIS;
using panzer::IP;
using panzer::Dim;

namespace charon {

template<typename EvalT, typename Traits>
class Integrator_GradBasisDotGradDOF
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_GradBasisDotGradDOF(
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

  PHX::MDField<ScalarT,Cell,BASIS> residual;

  PHX::MDField<const ScalarT,Cell,IP,Dim> flux;

  std::vector<PHX::MDField<const ScalarT,Cell,IP> > field_multipliers;

  std::size_t num_nodes;

  std::size_t num_qp;

  std::size_t num_dim;

  double multiplier;

  std::string basis_name;
  std::size_t basis_index;

  Kokkos::DynRankView<ScalarT,PHX::Device> tmp;

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Integrator_GradBasisDotGradDOF


}

#endif
