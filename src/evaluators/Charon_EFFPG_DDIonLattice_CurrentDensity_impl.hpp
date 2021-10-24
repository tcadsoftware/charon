
#ifndef CHARON_EFFPG_DDIONLATTICE_CURRENTDENSITY_IMPL_HPP
#define CHARON_EFFPG_DDIONLATTICE_CURRENTDENSITY_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Shards_CellTopology.hpp"

#include "Charon_Names.hpp"

// Calculate the EFFPG version of current density at IPs for the coupled DD+Ion+Lattice equations
// When comparing suzey's PDE equations to that used in Pavel Bochev's paper,
// the edge velocity is defined as
// vn = -mun*En - mun*\grad_T = mun*\grad_phi - mun*\grad_T for electrons,
// vp = -mup*Ep + mup*\grad_T = mup*\grad_phi + mup*\grad_T for holes,
// vi = -vi + ui*T*Si*\grad_T = -vi + DTi*\grad_T for ions/vacancies, ui*T = Di (Einstein relation).
// Here En = Ep = -\grad_phi is assumed, vi is computed at the center of a primary edge
// by the Mobility_RigidPointIon evaluator.
// See Suzey's Charon Notebook IV pg.202-203 for more details.

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
EFFPG_DDIonLattice_CurrentDensity<EvalT, Traits>::
EFFPG_DDIonLattice_CurrentDensity(
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

  // Hcurl basis for stabilization
  hcurl_basis_name = "HCurl:1:" + ir->getName();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  num_edges = edge_scalar->dimension(1);
  cellType = cellTopoInfo->getCellTopology();

  // Determine if want to include the \grad_T contribution
  bTempGrad = true;  // default
  if (p.isParameter("Temperature Gradient"))
    bTempGrad = p.get<bool>("Temperature Gradient");

  // obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;

  // carrier-dependent fields
  if (carrType == "Electron")
  {
    diffSign = 1.0;
    current_density = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_curr_density,ip_vector);
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.elec_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.elec_diff_coeff,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.edensity,node_scalar);
    electric_potential = MDField<const ScalarT,Cell,BASIS>(n.dof.phi,node_scalar);
    intrin_fermi = MDField<const ScalarT,Cell,BASIS>(n.field.intrin_fermi,node_scalar);
    bandgap = MDField<const ScalarT,Cell,BASIS>(n.field.band_gap,node_scalar);
    effbandgap = MDField<const ScalarT,Cell,BASIS>(n.field.eff_band_gap,node_scalar);
  }
  else if (carrType == "Hole")
  {
    diffSign = -1.0;
    current_density = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_curr_density,ip_vector);
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.hole_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.hole_diff_coeff,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.hdensity,node_scalar);
    electric_potential = MDField<const ScalarT,Cell,BASIS>(n.dof.phi,node_scalar);
    intrin_fermi = MDField<const ScalarT,Cell,BASIS>(n.field.intrin_fermi,node_scalar);
    bandgap = MDField<const ScalarT,Cell,BASIS>(n.field.band_gap,node_scalar);
    effbandgap = MDField<const ScalarT,Cell,BASIS>(n.field.eff_band_gap,node_scalar);
  }
  else if (carrType == "Ion")
  {
    diffSign = -1.0;
    current_density = MDField<ScalarT,Cell,IP,Dim>(n.field.ion_curr_density,ip_vector);
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.ion_mobility,edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.ion_diff_coeff,edge_scalar);
    thermodiff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.ion_thermodiff_coeff,edge_scalar);
    ion_velocity = MDField<const ScalarT,Cell,Edge>(n.field.ion_velocity,edge_scalar);
    carrier_density = MDField<const ScalarT,Cell,BASIS>(n.dof.iondensity,node_scalar);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole or Ion!");

  // evaluated field
  this->addEvaluatedField(current_density);

  // dependent fields
  latt_temp = MDField<const ScalarT,Cell,BASIS>(n.dof.latt_temp,node_scalar);
  this->addDependentField(latt_temp);

  this->addDependentField(mobility);
  this->addDependentField(diff_coeff);
  this->addDependentField(carrier_density);

  if ((carrType == "Electron") || (carrType == "Hole"))
  {
    this->addDependentField(electric_potential);
    this->addDependentField(intrin_fermi);
    this->addDependentField(bandgap);
    this->addDependentField(effbandgap);
  }
  else if (carrType == "Ion")
  {
    this->addDependentField(thermodiff_coeff);
    this->addDependentField(ion_velocity);
  }

  std::string name = "EFFPG_DDIonLattice_Current_Density";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
EFFPG_DDIonLattice_CurrentDensity<EvalT, Traits>::
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
EFFPG_DDIonLattice_CurrentDensity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // const factor in the effPot calculation
  // ScalarT factor = 0.5 / V0 * diffSign;  // + for e and - for h

  // loop over cells (hold for 2D and 3D)
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // zero out the array
    for (int ip = 0; ip < num_ips; ++ip)
      for (int dim = 0; dim < num_dims; ++dim)
        current_density(cell, ip, dim) = 0.0;

    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // get local node coordinates
      double x0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,0);
      double x1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,0);
      double y0 = 0.0, y1 = 0.0;
      double z0 = 0.0, z1 = 0.0;
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
      double edgeLen = std::sqrt(pow((x0-x1),2.0) + pow((y0-y1),2.0) + pow((z0-z1),2.0) );

      // get average diff. coeff. along the edge
      const ScalarT& edgeDiffCoeff = diff_coeff(cell,edge);

      // get average mobility at the center of the edge
      const ScalarT& edgeMob = mobility(cell,edge);

      // get the band gap narrowing in unit of eV
      // ScalarT dEg0 = bandgap(cell,node0) - effbandgap(cell,node0);
      // ScalarT dEg1 = bandgap(cell,node1) - effbandgap(cell,node1);

      // compute the effective potential (scaled)
      // ScalarT effPot0 = electric_potential(cell,node0) ;
      // ScalarT effPot1 = electric_potential(cell,node1) ;
      // The following effPot calculation was added on 12/05/2017 by Suzey Gao
      // ScalarT effPot0 = -intrin_fermi(cell,node0)/V0 + factor * dEg0;
      // ScalarT effPot1 = -intrin_fermi(cell,node1)/V0 + factor * dEg1;
      ScalarT effPot0 = -intrin_fermi(cell,node0)/V0 ;
      ScalarT effPot1 = -intrin_fermi(cell,node1)/V0 ;

      // get the lattice temperature
      const ScalarT& latt0 = latt_temp(cell,node0);
      const ScalarT& latt1 = latt_temp(cell,node1);

      // compute average velocity at the center of the edge (parallel to the edge)
      ScalarT barEdgeVel = 0;
      if ( (carrType == "Electron") || (carrType == "Hole") )
      {
        barEdgeVel = edgeMob*(effPot1-effPot0) /edgeLen - diffSign*edgeMob*(latt1 - latt0) /edgeLen;

        // Note that the barEdgeVel computed above should be the negative of the coefficient of the carrier density in the DDLattice equations.
        // This is because the parameters and formula used below (e.g., edgeAlpha, EFFPG current_density) are taken from Pavel et. al. 2013 which uses
        // J_n = -vel*n + diff*\grad_n in the analysis; Since our electron equation reads J_n = (\mu*E + \mu*\grad_T)*n + D*n, we must calculate
        // barEdgeVel=-\mu*E - \mu*\grad_T for electrons, etc.
      }
      else if (carrType == "Ion")
        barEdgeVel = -ion_velocity(cell,edge) + thermodiff_coeff(cell,edge) * (latt1 - latt0) / edgeLen;

      // compute edge Peclet number
      ScalarT edgeAlpha = diffSign * barEdgeVel * edgeLen / (edgeDiffCoeff * 2.0);

      // compute nodal coefficients for this edge
      ScalarT edgeCoef0 = edgeDiffCoeff * diffSign;  // pure diffusion limit
      ScalarT edgeCoef1 = edgeDiffCoeff * diffSign;

      if (std::abs( Sacado::ScalarValue<ScalarT>::eval(barEdgeVel) ) > 0.0)  // diffusion + advection/drift
      {
        ScalarT coth = 1.0 / tanh(edgeAlpha);
        edgeCoef0 = barEdgeVel * edgeLen / 2.0 * (coth + 1.0);
        edgeCoef1 = barEdgeVel * edgeLen / 2.0 * (coth - 1.0);
      }

      // loop over integration points and set current density
      for (int ip = 0; ip < num_ips; ++ip)
      {
        for (int dim = 0; dim < num_dims; ++dim)
        {
          current_density(cell, ip, dim) += ( carrier_density(cell,node1)*edgeCoef1 - carrier_density(cell,node0)*edgeCoef0 ) *
                                            ( workset.bases[hcurl_basis_index])->basis_vector(cell, edge, ip, dim);
        }
      }

    } // end edge loop

  }  // end of cell loop

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
EFFPG_DDIonLattice_CurrentDensity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<bool>("Temperature Gradient", true, "Turn on the temperature gradient contribution by default !");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
