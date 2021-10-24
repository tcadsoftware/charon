
#ifndef CVFEM_INTEGRATOR_SUBCVFLUXDOTNORM_IMPL_HPP
#define CVFEM_INTEGRATOR_SUBCVFLUXDOTNORM_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Integrator_SubCVFluxDotNorm<EvalT, Traits>::
Integrator_SubCVFluxDotNorm(
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

  // Obtain the BASIS layout
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> node_scalar = basis->functional;
  num_nodes = node_scalar->dimension(1);

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
  num_edges = edge_vector->dimension(1);
  num_dims = edge_vector->dimension(2);

  // Integration rule
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_scalar->dimension(1);

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Evaluated field
  residual = MDField<ScalarT,Cell,BASIS>(p.get<string>("Residual Name"),node_scalar);
  this->addEvaluatedField(residual);

  // Dependent field
  subcv_edge_flux = MDField<const ScalarT,Cell,Edge,Dim>(p.get<string>("Flux Name"),edge_vector);
  this->addDependentField(subcv_edge_flux);

  // get multiplier
  multiplier = p.get<double>("Multiplier");

  std::string name = "Integrator_SubCV_SurfaceFluxDotNorm";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_SubCVFluxDotNorm<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_SubCVFluxDotNorm<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  residual.deep_copy(ScalarT(0.0));

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // loop over edges
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids for the edge
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // compute the subcv surface flux dot the norm
      ScalarT residualVal = 0.0;
      for (int dim = 0; dim < num_dims; ++dim)
      {
        const ScalarT& flux = subcv_edge_flux(cell, edge, dim);
        const ScalarT norm = workset.int_rules[int_rule_index]->weighted_normals(cell,edge,dim);
        residualVal += flux*norm;
      }
      residualVal *= multiplier;

      // scatter residual contributions to nodes
      residual(cell, node0) += residualVal;
      residual(cell, node1) += -residualVal;

    } // end edge loop

  } // end cell loop
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Integrator_SubCVFluxDotNorm<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Flux Name", "?");
  p->set("Multiplier", 1.0);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;
}

}

#endif
