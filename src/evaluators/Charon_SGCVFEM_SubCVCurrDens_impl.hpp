
#ifndef CHARON_SGCVFEM_SUBCVCURRDENS_IMPL_HPP
#define CHARON_SGCVFEM_SUBCVCURRDENS_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SGCVFEM_SubCVCurrDens<EvalT, Traits>::
SGCVFEM_SubCVCurrDens(
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
  RCP<BasisIRLayout> hcurl_basis = p.get<RCP<BasisIRLayout> >("Basis");
  hcurl_basis_name = hcurl_basis->name();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = hcurl_basis->getCellTopologyInfo();
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

  std::string name = "CVFEM-SG_SubCV_Edge_Current_Density";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_SubCVCurrDens<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  hcurl_basis_index = panzer::getBasisIndex(hcurl_basis_name,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_SubCVCurrDens<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // zero out flux and put subcv edge(2D)/face(3D) mid-points in array
    for (int edge = 0; edge < num_edges; ++edge)
    {
      for (int dim = 0; dim < num_dims; ++dim)
        subcv_currdens(cell, edge, dim) = 0.0;
    }

    // loop over primary edges (ie, edge basis functions) to sum up edge curr. dens.
    for (int jedge = 0; jedge < num_edges; ++jedge)
    {

      // evaluate curr.dens. at the subcv edges(2D)/faces(3D) mid-points
      // note: number of interior subcv edges/faces is always equal to the number
      // of primary edges for quad, tri, hex, and tet mesh elements.

      for (int iedge = 0; iedge < num_edges; ++iedge)
      {
        for (int dim = 0; dim < num_dims; ++dim)
        {
          // note: previously edge_currdens() contained 1/edgeLen and edgeBasisValues from Intrepid
          // also contains 1/edgeLen, since subcv_currdens should have only one 1/edgeLen,
          // we needed to multiply edgeLen on the RHS here. Now have removed 1/edgeLen in
          // edge_curdens() so no longer need to multiply here.

          subcv_currdens(cell,iedge,dim) += edge_currdens(cell,jedge)
                            * (workset.bases[hcurl_basis_index])->basis_vector(cell,jedge,iedge,dim);
        }
      }

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
SGCVFEM_SubCVCurrDens<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Carrier Type", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> hcurl_basis;
  p->set("Basis", hcurl_basis);

  return p;
}

}

#endif
