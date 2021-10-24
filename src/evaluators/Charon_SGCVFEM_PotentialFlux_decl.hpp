
#ifndef CHARON_SGCVFEM_POTENTIALFLUX_DECL_HPP
#define CHARON_SGCVFEM_POTENTIALFLUX_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;
using panzer::Edge;
using panzer::Dim;

namespace charon {

// evaluate the potential flux, Lambda2*rel_perm*grad(phi), at the midpoints of
// subcontrol volume edges(2D)/faces(3D) in the CVFEM-SG formulation

template<typename EvalT, typename Traits>
class SGCVFEM_PotentialFlux
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCVFEM_PotentialFlux(
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
  PHX::MDField<ScalarT,Cell,Edge,Dim> subcv_phi_flux;

  // input
  PHX::MDField<const ScalarT,Cell,Edge,Dim> subcv_edge_midpt;
  PHX::MDField<const ScalarT,Cell,BASIS> potential;
  PHX::MDField<const ScalarT,Cell,BASIS> rel_perm;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Lambda2;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  // for Edge
  int num_nodes;
  int num_dims;
  int num_edges;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCVFEM_PotentialFlux


}

#endif
