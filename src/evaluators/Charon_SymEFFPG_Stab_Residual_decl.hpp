
#ifndef CHARON_SYMEFFPG_STAB_RESIDUAL_DECL_HPP
#define CHARON_SYMEFFPG_STAB_RESIDUAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::BASIS;
using panzer::Edge;

// compute the symmetrized EFFPG stabilization residual

namespace charon {

template<typename EvalT, typename Traits>
class SymEFFPG_Stab_Residual
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SymEFFPG_Stab_Residual(
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

  // output
  PHX::MDField<ScalarT,Cell,BASIS> residual;

  // input
  PHX::MDField<const ScalarT,Cell,Edge> mobility;     // scaled
  PHX::MDField<const ScalarT,Cell,Edge> diff_coeff;
  PHX::MDField<const ScalarT,Cell,Edge> thermodiff_coeff;
  PHX::MDField<const ScalarT,Cell,Edge> ion_velocity;

  PHX::MDField<const ScalarT,Cell,BASIS> carrier_density;
  PHX::MDField<const ScalarT,Cell,BASIS> electric_potential;
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;

  // IP
  std::size_t int_rule_index;
  int int_rule_degree;
  int num_ips;
  int num_dims;

  // BASIS
  std::string basis_name;
  std::size_t basis_index;
  int num_nodes;
  int num_edges;

  // HCURL BASIS
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  double sign;

  // carrier type
  std::string carrType;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SymEFFPG_Stab_Residual


}

#endif
