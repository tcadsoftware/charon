
#ifndef CHARON_SGCVFEM_POTENTIALFLUX_IMPL_HPP
#define CHARON_SGCVFEM_POTENTIALFLUX_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SGCVFEM_PotentialFlux<EvalT, Traits>::
SGCVFEM_PotentialFlux(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Obtain the BASIS layout
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> node_scalar = basis->functional;
  num_nodes = node_scalar->dimension(1);
  basis_name = basis->name();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
  num_edges = edge_vector->dimension(1);
  num_dims = edge_vector->dimension(2);

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Evaluated field
  subcv_phi_flux = MDField<ScalarT,Cell,Edge,Dim>(p.get<string>("Flux Name"), edge_vector);
  this->addEvaluatedField(subcv_phi_flux);

  // Scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Lambda2 = scaleParams->scale_params.Lambda2;

  // Dependent fields
  potential = MDField<const ScalarT,Cell,BASIS>(p.get<string>("DOF Name"), node_scalar);
  rel_perm = MDField<const ScalarT,Cell,BASIS>(n.field.rel_perm, node_scalar);

  this->addDependentField(potential);
  this->addDependentField(rel_perm);

  std::string name = "CVFEM-SG_SubCV_Edge_Potential_Flux";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_PotentialFlux<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_PotentialFlux<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // zero out flux
    for (int edge = 0; edge < num_edges; ++edge)
    {
      for (int dim = 0; dim < num_dims; ++dim)
         subcv_phi_flux(cell, edge, dim) = 0.0;
    }

    // evaluate phi_flux at the interior subcv edges(2D)/faces(3D) mid-points
    // note: number of interior subcv edges/faces is always equal to the number
    // of primary edges for quad, tri, hex, and tet mesh elements.

    Kokkos::DynRankView<ScalarT,PHX::Device> midptValues = Kokkos::createDynRankView(subcv_phi_flux.get_static_view(),"midptValues",num_edges);
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // loop over primary cell nodes (= # of nodal basis) to compute grad(phi)
      for (int node = 0; node < num_nodes; ++node)
      {
        midptValues(edge) += rel_perm(cell,node)*Lambda2*
                             (workset.bases[basis_index])->basis_scalar(cell,node,edge);
        for (int dim = 0; dim < num_dims; ++dim)
         {
           subcv_phi_flux(cell,edge,dim) += (workset.bases[basis_index])->grad_basis(cell,node,edge,dim)*potential(cell,node);
         }

      }

      // calculate the phi_flux at subcv edges(2D)/faces(3D) mid-points
      for (int dim = 0; dim < num_dims; ++dim)
      {
        subcv_phi_flux(cell,edge,dim) *= midptValues(edge);
        // std::cout << "edge=" << edge << ", dim=" << dim << ", subcv_phi_flux=" << subcv_phi_flux(cell,edge,dim) << std::endl;
      }
    }

  }  // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SGCVFEM_PotentialFlux<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Flux Name", "??");
  p->set<std::string>("DOF Name", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
