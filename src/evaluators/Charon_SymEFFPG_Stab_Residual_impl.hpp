
#ifndef CHARON_SYMEFFPG_STAB_RESIDUAL_IMPL_HPP
#define CHARON_SYMEFFPG_STAB_RESIDUAL_IMPL_HPP

#include "Kokkos_ViewFactory.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SymEFFPG_Stab_Residual<EvalT, Traits>::
SymEFFPG_Stab_Residual(
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

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // BASIS
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> node_scalar = basis->functional;
  RCP<DataLayout> node_vector = basis->functional_grad;
  basis_name = basis->name();
  num_nodes = node_scalar->dimension(1);

  // HCurl basis for stabilization
  hcurl_basis_name = "HCurl:1:" + ir->getName();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  num_edges = edge_scalar->dimension(1);
  cellType = cellTopoInfo->getCellTopology();

  // obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // evaluated field
  residual = MDField<ScalarT,Cell,BASIS>(p.get<string>("Residual Name"),node_scalar);
  this->addEvaluatedField(residual);

  // carrier-dependent fields
  if (carrType == "Electron")
  {
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.elec_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.elec_diff_coeff,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.edensity,node_scalar);
    electric_potential = MDField<const ScalarT,Cell,BASIS>(n.dof.phi,node_scalar);
    latt_temp = MDField<const ScalarT,Cell,BASIS>(n.dof.latt_temp,node_scalar);
    sign = -1.0;
    this->addDependentField(mobility);
    this->addDependentField(diff_coeff);
    this->addDependentField(carrier_density);
    this->addDependentField(electric_potential);
    this->addDependentField(latt_temp);
  }

  else if (carrType == "Hole")
  {
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.hole_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.hole_diff_coeff,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.hdensity,node_scalar);
    electric_potential = MDField<const ScalarT,Cell,BASIS>(n.dof.phi,node_scalar);
    latt_temp = MDField<const ScalarT,Cell,BASIS>(n.dof.latt_temp,node_scalar);
    sign = 1.0;
    this->addDependentField(mobility);
    this->addDependentField(diff_coeff);
    this->addDependentField(carrier_density);
    this->addDependentField(electric_potential);
    this->addDependentField(latt_temp);
  }

  else if (carrType == "Ion")
  {
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.ion_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.ion_diff_coeff,edge_scalar);
    thermodiff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.ion_thermodiff_coeff,edge_scalar);
    ion_velocity = MDField<const ScalarT,Cell,Edge>(n.field.ion_velocity,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.iondensity,node_scalar);
    latt_temp = MDField<const ScalarT,Cell,BASIS>(n.dof.latt_temp,node_scalar);
    this->addDependentField(mobility);
    this->addDependentField(diff_coeff);
    this->addDependentField(thermodiff_coeff);
    this->addDependentField(ion_velocity);
    this->addDependentField(carrier_density);
    this->addDependentField(latt_temp);
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole or Ion!");

  std::string name = "SymEFFPG_Stab_Residual";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SymEFFPG_Stab_Residual<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
  hcurl_basis_index = panzer::getBasisIndex(hcurl_basis_name,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SymEFFPG_Stab_Residual<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  residual.deep_copy(ScalarT(0.0));  // set to 0

  auto k_view = residual.get_static_view();
  using KDRV = Kokkos::DynRankView<ScalarT,PHX::Device>;
  KDRV stab_flux_carr = Kokkos::createDynRankView(k_view,"my_test", workset.num_cells, num_ips, num_dims);
  KDRV stab_flux_nb = Kokkos::createDynRankView(k_view,"my_test2",workset.num_cells, num_nodes, num_ips, num_dims);
  KDRV stab_flux_dp = Kokkos::createDynRankView(k_view,"my_test3",workset.num_cells, num_nodes, num_ips);

  double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  double x1 = 0.0, y1 = 0.0, z1 = 0.0;

  // loop over cells (hold for 2D and 3D)
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // zero out the arrays
    for (int ip = 0; ip < num_ips; ++ip)
    {
      for (int dim = 0; dim < num_dims; ++dim)
      {
        stab_flux_carr(cell,ip,dim) = 0.0;
        for (int node = 0; node < num_nodes; ++node)
        {
          stab_flux_nb(cell,node,ip,dim) = 0.0;
          stab_flux_dp(cell,node,ip) = 0.0;
        }
      }
    }

    // loop over edges
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // get local node coordinates
      x0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,0);
      x1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,0);
      if (num_dims > 1)  // 2D or 3D
      {
        y0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,1);
        y1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,1);
      }
      if (num_dims > 2)  // 3D
      {
        z0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,2);
        z1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,2);
      }

      // compute the primary cell edge length
      double edgeLen = std::sqrt(std::pow((x0-x1),2.0) + std::pow((y0-y1),2.0) + std::pow((z0-z1),2.0) );

      // get average diff. coeff. along the edge
      const ScalarT& edgeDiffCoeff = diff_coeff(cell,edge);

      // get average mobility at the center of the edge
      const ScalarT& edgeMob = mobility(cell,edge);

      // get the lattice temperature at nodes
      const ScalarT& latt0 = latt_temp(cell,node0);
      const ScalarT& latt1 = latt_temp(cell,node1);

      // edge coefficient (corresponding to 'theta_tilde' in Pavel's paper)
      ScalarT edgeCoeff = 0.0;  // pure diffusion limit
      ScalarT barEdgeVel = 0.0;

      // compute average effective velocity at the center of an edge
      if ((carrType == "Electron") || (carrType == "Hole") )
      {
        const ScalarT& pot0 = electric_potential(cell,node0);
        const ScalarT& pot1 = electric_potential(cell,node1);

        // compute average velocity at the center of an edge (parallel to the edge)
        barEdgeVel = sign * edgeMob * (pot0 - pot1) / edgeLen  - edgeMob * (latt1 - latt0) / edgeLen;

        // barEdgeVel for holes has a different sign from that in EFFPG_DDIonLattice_CurrentDensity,
        // since we need to rewrite the Jp expression such that the diffusion term uses Dp*\grad_p
        // instead of the original -Dp*\grad_p. This is because we need to make sure D*(pcothp-1), the
        // 'theta' parameter in Pavel's sym-EFFPG paper, is positive to allow for sqrt() operation.
        // A similar argument holds for ions.
      }

      else if (carrType == "Ion")
      {
        // compute average velocity at the center of an edge
        barEdgeVel = ion_velocity(cell,edge) - thermodiff_coeff(cell,edge) * (latt1 - latt0)/edgeLen;
      }

      if (std::abs( Sacado::ScalarValue<ScalarT>::eval(barEdgeVel) ) > 0.0)  // diffusion + advection/drift
      {
        // compute edge Peclet number
        ScalarT edgePeclet = barEdgeVel * edgeLen / (2.0 * edgeDiffCoeff);

        ScalarT coth = 1.0 / std::tanh(edgePeclet);
        ScalarT edgeTheta = edgeDiffCoeff * (edgePeclet * coth - 1.0);

        if (Sacado::ScalarValue<ScalarT>::eval(edgeTheta) > 0.0)
        {
          // In theory, sign of edgeTheta is the same as edgeDiffCoeff. We rewrote the current density expressions
          // such that edgeDiffCoeff is always positive, hence edgeTheta should be always >= 0. However, due to
          // numerical errors, edgeTheta can become a small negative number.
          // Most importantly, when edgeTheta = 0, edgeCoeff = std::sqrt(edgeTheta) introduces "inf" and "nan" into
          // the Jacobian matrix through automatic differentiation. Hence, set edgeCoeff=0 when edgeTheta <= 0.0.

          edgeCoeff = std::sqrt(edgeTheta);
        }
      }

      // compute stab_flux_carr and stab_flux_nb (node basis)
      for (int ip = 0; ip < num_ips; ++ip)
      {
        for (int dim = 0; dim < num_dims; ++dim)
        {
          stab_flux_carr(cell,ip,dim) += edgeCoeff * (carrier_density(cell,node1) - carrier_density(cell,node0)) *
                                          (workset.bases[hcurl_basis_index])->basis_vector(cell, edge, ip, dim);

          // no need to loop over nodes, since for nodes other than node0 and node1, endValue-startValue = 0
          stab_flux_nb(cell,node0,ip,dim) += edgeCoeff * (-1.0) * (workset.bases[hcurl_basis_index])->basis_vector(cell, edge, ip, dim);
          stab_flux_nb(cell,node1,ip,dim) += edgeCoeff * (+1.0) * (workset.bases[hcurl_basis_index])->basis_vector(cell, edge, ip, dim);

        }
      }
    } // end of edge loop

    // compute the dot product of stab_flux_carr(cell,ip,dim) and stab_flux_nb(cell,basis,ip,dim)
    for (int basis = 0; basis < num_nodes; ++basis)
      for (int ip = 0; ip < num_ips; ++ip)
        for (int dim = 0; dim < num_dims; ++dim)
          stab_flux_dp(cell,basis,ip) += stab_flux_carr(cell,ip,dim) * stab_flux_nb(cell,basis,ip,dim);

  }  // end of cell loop

  const panzer::IntegrationValues2<double> & iv = *workset.int_rules[int_rule_index];
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_nodes; ++ basis)
    {
      residual(cell,basis) = stab_flux_dp(cell,basis,0) * iv.weighted_measure(cell,0);
      for (int ip = 1; ip < num_ips; ++ip)
        residual(cell,basis) += stab_flux_dp(cell,basis,ip) * iv.weighted_measure(cell,ip);
    }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SymEFFPG_Stab_Residual<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Carrier Type", "?");
  p->set<std::string>("Residual Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  return p;
}

}

#endif

