
#ifndef CHARON_SGCHARON1_SUBCVCURRDENS_DECL_HPP
#define CHARON_SGCHARON1_SUBCVCURRDENS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Edge;
using panzer::Dim;

namespace charon {

// compute the current density vector along primary cell edge, which is approximated
// as subcv edge current density, the same SG implementation as in Charon1.

template<typename EvalT, typename Traits>
class SGCharon1_SubCVCurrDens
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCharon1_SubCVCurrDens(
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
  PHX::MDField<ScalarT,Cell,Edge,Dim> subcv_currdens; // curr. dens. (vector)

  // input
  PHX::MDField<const ScalarT,Cell,Edge> edge_currdens;  // scalar

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  std::string carrType;

  // for Edge
  int num_dims;
  int num_edges;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCharon1_SubCVCurrDens


}

#endif
