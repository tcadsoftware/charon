
#ifndef CHARON_SGCVFEM_SUBCVCURRDENS_DECL_HPP
#define CHARON_SGCVFEM_SUBCVCURRDENS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Edge;
using panzer::Dim;

namespace charon {

// compute the current density evaluated at the center of subcontrol volume subcell
// (edge in 2D and face in 3D) inside a primary mesh cell in the CVFEM-SG formulation

template<typename EvalT, typename Traits>
class SGCVFEM_SubCVCurrDens
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCVFEM_SubCVCurrDens(
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
  PHX::MDField<const ScalarT,Cell,Edge> edge_currdens;

  // for basis points
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  // reference edge length
  double refEdgeLen;

  std::string carrType;

  // for Edge
  int num_dims;
  int num_edges;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCVFEM_SubCVCurrDens


}

#endif
