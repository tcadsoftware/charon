
#ifndef CHARON_SGCHARON1_SUBCVCURRDENS_IMPL_HPP
#define CHARON_SGCHARON1_SUBCVCURRDENS_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_Basis.hpp"

#include "Charon_EdgeBasisFactory.hpp"
#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SGCharon1_SubCVCurrDens<EvalT, Traits>::
SGCharon1_SubCVCurrDens(
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
  basis_name = basis->name();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
  num_edges = edge_vector->dimension(1);
  num_dims = edge_vector->dimension(2);

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    subcv_currdens = MDField<ScalarT,Cell,Edge,Dim>(n.field.elec_curr_dens_cvedge, edge_vector);
    edge_currdens = MDField<const ScalarT,Cell,Edge>(n.field.elec_edge_currdens, edge_scalar);
  }
  else if (carrType == "Hole")
  {
    subcv_currdens = MDField<ScalarT,Cell,Edge,Dim>(n.field.hole_curr_dens_cvedge, edge_vector);
    edge_currdens = MDField<const ScalarT,Cell,Edge>(n.field.hole_edge_currdens, edge_scalar);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Add evaluated fields
  this->addEvaluatedField(subcv_currdens);

  // Add dependent fields
  this->addDependentField(edge_currdens);

  std::string name = "Charon1-SG_SubCV_Edge_Current_Density";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCharon1_SubCVCurrDens<EvalT, Traits>::
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
SGCharon1_SubCVCurrDens<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // loop over primary edges
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // compute primary cell edge length and tangent
      Kokkos::DynRankView<double,PHX::Device> startpoint = Kokkos::DynRankView<double,PHX::Device>("startpoint",num_dims);
      Kokkos::DynRankView<double,PHX::Device> endpoint = Kokkos::DynRankView<double,PHX::Device>("endpoint",num_dims);
      Kokkos::DynRankView<double,PHX::Device> distance = Kokkos::DynRankView<double,PHX::Device>("distance",num_dims);
      double edgeLen = 0.0;

      for (int dim = 0; dim < num_dims; ++dim)
      {
        // get local node coordinate
        startpoint(dim) = (workset.bases[basis_index])->basis_coordinates(cell,node0,dim);
        endpoint(dim) = (workset.bases[basis_index])->basis_coordinates(cell,node1,dim);

        // get local coordinate's distance
        distance(dim) = endpoint(dim) - startpoint(dim);

        // get primary cell edge length
        edgeLen += distance(dim) * distance(dim);
      }
      edgeLen = std::sqrt(edgeLen);

      // compute curr. dens. vector along primary edge
      for (int dim = 0; dim < num_dims; ++dim)
        subcv_currdens(cell,edge,dim) = edge_currdens(cell,edge) * distance(dim)/edgeLen;

    }  // end of loop over primary edges

  }  // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SGCharon1_SubCVCurrDens<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Carrier Type", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  return p;
}

}

#endif
