
#ifndef CHARON_SGCVFEM_SUBCVCENTROIDVECTOR_IMPL_HPP
#define CHARON_SGCVFEM_SUBCVCENTROIDVECTOR_IMPL_HPP

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
SGCVFEM_SubCVCentroidVector<EvalT, Traits>::
SGCVFEM_SubCVCentroidVector(
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

  // Basis
  RCP<BasisIRLayout> hcurl_basis = p.get<RCP<BasisIRLayout> >("Basis");
  hcurl_basis_name = hcurl_basis->name();

  // Integration rule for subCV centroid data layout
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_vector = ir->dl_vector;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // Edge data layout
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = hcurl_basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  num_edges = edge_scalar->dimension(1);

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    subcv_centvec = MDField<ScalarT,Cell,IP,Dim>(p.get<string>("Vector Name"), ip_vector);
    edge_currdens = MDField<const ScalarT,Cell,Edge>(n.field.elec_edge_currdens, edge_scalar);
  }
  else if (carrType == "Hole")
  {
    subcv_centvec = MDField<ScalarT,Cell,IP,Dim>(p.get<string>("Vector Name"), ip_vector);
    edge_currdens = MDField<const ScalarT,Cell,Edge>(n.field.hole_edge_currdens, edge_scalar);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Add evaluated fields
  this->addEvaluatedField(subcv_centvec);

  // Add dependent fields
  this->addDependentField(edge_currdens);

  std::string name = "CVFEM-SG_SubCV_Centroid_Vector";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_SubCVCentroidVector<EvalT, Traits>::
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
SGCVFEM_SubCVCentroidVector<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // zero out vector at subcontrol volume centroids
    for (int node = 0; node < num_ips; ++node)
    {
      for (int dim = 0; dim < num_dims; ++dim)
        subcv_centvec(cell, node, dim) = 0.0;
    }

    // loop over primary edges (ie, edge basis functions) to sum up edge curr. dens.
    for (int iedge = 0; iedge < num_edges; ++iedge)
    {

      // evaluate curr.dens. at the subcv centroids
      // note: number of subcv centroids is equal to the number of primary
      // nodes for quad, tri, hex, and tet mesh elements.

      for (int ip = 0; ip < num_ips; ++ip)
      {
        for (int dim = 0; dim < num_dims; ++dim)
        {

          subcv_centvec(cell,ip,dim) += edge_currdens(cell,iedge)
                            * (workset.bases[hcurl_basis_index])->basis_vector(cell,iedge,ip,dim);

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
SGCVFEM_SubCVCentroidVector<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Carrier Type", "??");
  p->set<std::string>("Vector Name", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> hcurl_basis;
  p->set("Basis", hcurl_basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;
}

}

#endif
