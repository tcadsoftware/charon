
#ifndef CHARON_DRIVEFORCE_IMPL_HPP
#define CHARON_DRIVEFORCE_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"

#include "Shards_CellTopology.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"


/*
The effective electric field Fn\p,eff = grad(negEffPot_n\p) (scaled), where
negEffPot_n\p = [Ei -\+ 0.5*dEg -\+ kbT*log(nie/nie0)] /V0, (- for n, + for p)
              = [Ei -\+ 0.5*dEg -\+ 0.5*kbT*log(gamma_n*gamma_p)] /V0.
nie = nie0*sqrt(gamma_n*gamma_p), nie0 = equilibrium effective intrinsic conc.,
and gamma_n\p = electron\hole degeneracy factor.
The calculation is valid for Boltzmann and Fermi-Dirac statistics, and
for BGN = On and Off cases. dEg is equal to 0 if BGN = Off.

Reference: D. Schroeder, T. Ostermann and O. Kalz, "Comparison of
transport models for the simulation of degenerate semiconductors,"
Semicond. Sci. Technol.9 (1994) 364-369.

This evaluator basically implements Eqn.(30)-(31) in the paper except that we
evaluate the scaled version. Note that Eqn.(30)-(31) are applicable to general
devices (homo- and hetero-geneous, the latter has constant mole fraction)
for Boltzmann and Fermi-Dirac statistics.

As demonstrated in the paper, the Boltzmann statistics works well in most
situations and even in certain degenerate regions, provided they are charge
neutral and the BGN effect is included through certain BGN models.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SGCVFEM_CentroidDriveForce<EvalT, Traits>::
SGCVFEM_CentroidDriveForce(
  const Teuchos::ParameterList& p) {
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
  RCP<BasisIRLayout> hcurl_basis = p.get<RCP<BasisIRLayout> >("HCurlBasis");
  RCP<BasisIRLayout> hgrad_basis = p.get<RCP<BasisIRLayout> >("HGradBasis");
  RCP<DataLayout> hgrad_basis_scalar = hgrad_basis->functional;
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

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // selectes driving force to compute
  drForceType = p.get<string>("Driving Force");
  if(drForceType != "EffectiveField" && drForceType != "GradQuasiFermi" &&
     drForceType != "GradPotential")
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Driving Force!");

  if(carrType == "Electron") 
    {
      if(drForceType == "EffectiveField") {
	field = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_efield, ip_vector);
      } 
      else if(drForceType == "GradQuasiFermi") 
      {
	field = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_grad_qfp, ip_vector);
	density = MDField<const ScalarT,Cell,Point>(n.dof.edensity,hgrad_basis_scalar);
      } 
      else if(drForceType == "GradPotential") 
      {
	field = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_grad_negpot, ip_vector);
      }
    sign = -1;
  } 
  else if(carrType == "Hole")
  {
    if(drForceType == "EffectiveField") {
      field = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_efield, ip_vector);
    } 
    else if(drForceType == "GradQuasiFermi") 
    {
      field = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_grad_qfp, ip_vector);
      density = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,hgrad_basis_scalar);
    } 
    else if(drForceType == "GradPotential") 
    {
      field = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_grad_negpot, ip_vector);
    }
    sign = 1;
  } 
 
  // Carrier-independent fields
  if(drForceType == "GradPotential")
    //potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,hcurl_basis_scalar);
    potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,hgrad_basis_scalar);
  if(drForceType != "GradPotential") {
    intrinfermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi,hgrad_basis_scalar);
    bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap,hgrad_basis_scalar);
    effbandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,hgrad_basis_scalar);
    latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,hgrad_basis_scalar);
    // Degeneracy factors
    elec_degfactor = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor,hgrad_basis_scalar);
    hole_degfactor = MDField<const ScalarT,Cell,Point>(n.field.hole_deg_factor,hgrad_basis_scalar);
  }

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  T0 = scaleParams->scale_params.T0;

  // Add evaluated fields
  this->addEvaluatedField(field);
 
  // Add dependent fields
  if(drForceType == "GradPotential")
    this->addDependentField(potential);
  if(drForceType == "GradQuasiFermi")
    this->addDependentField(density);
  if(drForceType != "GradPotential") {
    this->addDependentField(intrinfermi);
    this->addDependentField(bandgap);
    this->addDependentField(effbandgap);
    this->addDependentField(latt_temp);
    this->addDependentField(elec_degfactor);
    this->addDependentField(hole_degfactor);
  }

  std::string name = "SGCVFEM_DrivingForce";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_CentroidDriveForce<EvalT, Traits>::
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
SGCVFEM_CentroidDriveForce<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset) {
  using panzer::index_t;

  // Obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb; // Boltzmann constant in [eV/K]

  // loop over cells
  for(index_t cell = 0; cell < workset.num_cells; ++cell) {
    // zero out vector at subcontrol volume centroids
    for(int node = 0; node < num_ips; ++node) {
      for(int dim = 0; dim < num_dims; ++dim)
	field(cell, node, dim) = 0.0;
    }

    // loop over primary edges (ie, edge basis functions) to sum up edge eff field
    for(int iedge = 0; iedge < num_edges; ++iedge) {
      // get local node ids: first index 1 for edge
      // (0 for vertex, 2 for face, 3 for volume)
      const int node0 = cellType->getNodeMap(1,iedge,0);
      const int node1 = cellType->getNodeMap(1,iedge,1);

      // get local node coordinates
      double x0 = workset.cell_vertex_coordinates(cell,node0,0);
      double x1 = workset.cell_vertex_coordinates(cell,node1,0);
      double y0 = 0.0, y1 = 0.0;
      double z0 = 0.0, z1 = 0.0;
      if(num_dims > 1)  { // 2D or 3D
        y0 = workset.cell_vertex_coordinates(cell,node0,1);
        y1 = workset.cell_vertex_coordinates(cell,node1,1);
      }
      if (num_dims > 2) { // 3D
        z0 = workset.cell_vertex_coordinates(cell,node0,2);
        z1 = workset.cell_vertex_coordinates(cell,node1,2);
      }

      // compute the primary cell edge length
      double edgeLen =
        std::sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));

      // variable to store selected driving force on edge
      ScalarT val_edge = 0.;

      if(drForceType == "EffectiveField" || drForceType == "GradQuasiFermi") {
        // obtain fields at local nodes
        const ScalarT& lattT0 = latt_temp(cell,node0)*T0; // [K]
        const ScalarT& lattT1 = latt_temp(cell,node1)*T0; // [K]

        const ScalarT& kbT0 = kbBoltz*lattT0; // [eV]
        const ScalarT& kbT1 = kbBoltz*lattT1; // [eV]

        const ScalarT& Eg0 = bandgap(cell,node0); // [eV]
        const ScalarT& Eg1 = bandgap(cell,node1); // [eV]
        const ScalarT& effEg0 = effbandgap(cell,node0); // [eV]
        const ScalarT& effEg1 = effbandgap(cell,node1); // [eV]
        const ScalarT deltaEg0 = Eg0 - effEg0;
        const ScalarT deltaEg1 = Eg1 - effEg1;

        const ScalarT& Ei0 = intrinfermi(cell,node0); // [eV]
        const ScalarT& Ei1 = intrinfermi(cell,node1); // [eV]

        const ScalarT& gamma_n0 = elec_degfactor(cell,node0);
        const ScalarT& gamma_n1 = elec_degfactor(cell,node1);

        const ScalarT& gamma_p0 = hole_degfactor(cell,node0);
        const ScalarT& gamma_p1 = hole_degfactor(cell,node1);

        ScalarT degterm0 = 0.5*kbT0*std::log(gamma_n0 * gamma_p0);  // = 0 for Boltzmann statistics
        ScalarT degterm1 = 0.5*kbT1*std::log(gamma_n1 * gamma_p1);

        // compute the effective potential at local nodes [scaled]
        ScalarT effPot0 = (-sign*0.5*deltaEg0 - Ei0 - sign*degterm0) / V0;
        ScalarT effPot1 = (-sign*0.5*deltaEg1 - Ei1 - sign*degterm1) / V0;
        ScalarT efield_edge = (effPot1 - effPot0)/edgeLen;
        val_edge = -efield_edge;

        if(drForceType == "GradQuasiFermi") {
          const ScalarT& dens0 = density(cell,node0);
          const ScalarT& dens1 = density(cell,node1);

          // d(Efn) = elec_effField + d(n)/n
          // d(Efp) = hole_effField - d(p)/p
          ScalarT grad_qfp_edge =
            efield_edge - sign*(dens1 - dens0)/edgeLen/(0.5*(dens0 + dens1));
          val_edge = -grad_qfp_edge;
        }
      } else if(drForceType == "GradPotential") {
        const ScalarT& pot0 = potential(cell,node0);
        const ScalarT& pot1 = potential(cell,node1);
        ScalarT grad_pot_edge = (pot1 - pot0)/edgeLen;
        val_edge = -grad_pot_edge;
      }

      // evaluate driving force at the subcv centroids
      // note: number of subcv centroids is equal to the number of primary
      // nodes for quad, tri, hex, and tet mesh elements.
      for(int ip = 0; ip < num_ips; ++ip) {
        for(int dim = 0; dim < num_dims; ++dim) 
	  field(cell,ip,dim) += val_edge
              * (workset.bases[hcurl_basis_index])->basis_vector(cell,iedge,ip,dim)*edgeLen;
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
SGCVFEM_CentroidDriveForce<EvalT, Traits>::getValidParameters() const {
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Carrier Type", "??");
  p->set<std::string>("Driving Force", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> hcurl_basis;
  p->set("HCurlBasis", hcurl_basis);

  Teuchos::RCP<panzer::BasisIRLayout> hgrad_basis;
  p->set("HGradBasis", hgrad_basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
