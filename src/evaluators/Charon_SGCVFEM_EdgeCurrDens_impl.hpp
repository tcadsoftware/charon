
#ifndef CHARON_SGCVFEM_EDGECURRDENS_IMPL_HPP
#define CHARON_SGCVFEM_EDGECURRDENS_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

/*
Use the formulation from the reference by D. Schroeder, T. Ostermann and O. Kalz,
"Comparison of transport models for the simulation of degenerate semiconductors,"
Semicond. Sci. Technol.9 (1994) 364-369.

In this formulation, D/u = kbT/q always holds, so to include Fermi-Dirac statistics,
we just need to add the additional, kbT/q*log(nie/nie0) = 0.5*kbT/q*log(gamma_n*gamma_p)
term to the calculation of nodal effective potentials, provided that Ei (intrinsic
Fermi energy level) already includes the FD effect.

The implementation here is valid for Boltzmann and Fermi-Dirac statistics under
isothermal DD simulations.
*/



namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SGCVFEM_EdgeCurrDens<EvalT, Traits>::
SGCVFEM_EdgeCurrDens(
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
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
  num_edges = edge_vector->dimension(1);
  num_dims = edge_vector->dimension(2);
  //num_edges = cellTopoInfo->getNumEdges();  // also work
  //num_dims = cellTopoInfo->getDimension();

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    edge_currdens = MDField<ScalarT,Cell,Edge>(n.field.elec_edge_currdens, edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.elec_diff_coeff, edge_scalar);
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.elec_mobility, edge_scalar);
    density = MDField<const ScalarT,Cell,BASIS>(n.dof.edensity, basis_scalar);
    sign = 1.0;
    // Density Gradient
    useEQC = p.get< bool >("Use Electron Quantum Correction");
    if( useEQC ) {
      eqp = MDField<const ScalarT,Cell,BASIS>(n.dof.elec_qpotential, basis_scalar);
      this->addDependentField(eqp);
    }
  }
  else if (carrType == "Hole")
  {
    edge_currdens = MDField<ScalarT,Cell,Edge>(n.field.hole_edge_currdens, edge_scalar);
    diff_coeff = MDField<const ScalarT,Cell,Edge>(n.field.hole_diff_coeff, edge_scalar);
    mobility = MDField<const ScalarT,Cell,Edge>(n.field.hole_mobility, edge_scalar);
    density = MDField<const ScalarT,Cell,BASIS>(n.dof.hdensity, basis_scalar);
    sign = -1.0;
    useHQC = p.get< bool >("Use Hole Quantum Correction");
    if( useHQC ) {
      hqp = MDField<const ScalarT,Cell,BASIS>(n.dof.hole_qpotential, basis_scalar);
      this->addDependentField(hqp);
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Carrier-Indepdenent fields
  intrin_fermi = MDField<const ScalarT,Cell,BASIS>(n.field.intrin_fermi, basis_scalar);
  bandgap = MDField<const ScalarT,Cell,BASIS>(n.field.band_gap, basis_scalar);
  eff_bandgap = MDField<const ScalarT,Cell,BASIS>(n.field.eff_band_gap, basis_scalar);
  elec_degfac = MDField<const ScalarT,Cell,BASIS>(n.field.elec_deg_factor, basis_scalar);
  hole_degfac = MDField<const ScalarT,Cell,BASIS>(n.field.hole_deg_factor, basis_scalar);
  latt_temp = MDField<const ScalarT,Cell,BASIS>(n.field.latt_temp, basis_scalar);


  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  T0 = scaleParams->scale_params.T0;

  // Add evaluated field
  this->addEvaluatedField(edge_currdens);

  // Add dependent fields
  this->addDependentField(diff_coeff);
  this->addDependentField(mobility);
  this->addDependentField(density);
  this->addDependentField(intrin_fermi);
  this->addDependentField(bandgap);
  this->addDependentField(eff_bandgap);
  this->addDependentField(elec_degfac);
  this->addDependentField(hole_degfac);
  this->addDependentField(latt_temp);
  

  std::string name = "CVFEM-SG_Primary_Edge_Current_Density";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SGCVFEM_EdgeCurrDens<EvalT, Traits>::
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
SGCVFEM_EdgeCurrDens<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // loop over primary edges
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // obtain mobility and diff. coeff. at the center of a primary edge
      ScalarT mob = mobility(cell,edge);  // [scaled]
      ScalarT diff = sign * diff_coeff(cell,edge);
      ScalarT beta = diff / mob;

      // obtain fields at local nodes
      ScalarT dEg0 = bandgap(cell,node0) - eff_bandgap(cell,node0) ;  // [eV]
      ScalarT dEg1 = bandgap(cell,node1) - eff_bandgap(cell,node1) ;
      ScalarT iEf0 = intrin_fermi(cell,node0);  // [eV]
      ScalarT iEf1 = intrin_fermi(cell,node1);
      ScalarT edegfac0 = elec_degfac(cell,node0);  // [unitless]
      ScalarT edegfac1 = elec_degfac(cell,node1);
      ScalarT hdegfac0 = hole_degfac(cell,node0);
      ScalarT hdegfac1 = hole_degfac(cell,node1);

      ScalarT kbT0 = kbBoltz*latt_temp(cell,node0)*T0;  // [eV]
      ScalarT kbT1 = kbBoltz*latt_temp(cell,node1)*T0;
      ScalarT degterm0 = 0.5*kbT0*log(edegfac0*hdegfac0);  // = 0 for Boltzmann statistics
      ScalarT degterm1 = 0.5*kbT1*log(edegfac1*hdegfac1);
      ScalarT effPot0, effPot1;
      if ( carrType == "Electron" && useEQC ) {
          // Electron Quantum Correction
          ScalarT eqp0 = eqp(cell,node0); 
          ScalarT eqp1 = eqp(cell,node1);
          effPot0 = (sign*0.5*dEg0 - iEf0 + sign*degterm0 ) / V0 + eqp0;
          effPot1 = (sign*0.5*dEg1 - iEf1 + sign*degterm1 ) / V0 + eqp1;
      } else if ( carrType == "Hole" && useHQC ) {
          // Hole Quantum Correction
          ScalarT hqp0 = hqp(cell,node0); 
          ScalarT hqp1 = hqp(cell,node1);
          effPot0 = (sign*0.5*dEg0 - iEf0 + sign*degterm0 ) / V0 - hqp0;
          effPot1 = (sign*0.5*dEg1 - iEf1 + sign*degterm1 ) / V0 - hqp1;
      } else {
          // compute the effective potential at local nodes [scaled]
          effPot0 = (sign*0.5*dEg0 - iEf0 + sign*degterm0 ) / V0;
          effPot1 = (sign*0.5*dEg1 - iEf1 + sign*degterm1 ) / V0;
      }


      //ScalarT effPot0 = (sign*0.5*dEg0 - iEf0) / V0;
      //ScalarT effPot1 = (sign*0.5*dEg1 - iEf1) / V0;
      // std::cout << "degterm0 = " << degterm0 << ", degterm1 = " << degterm1 << std::endl;

      // get edge Reynolds number
      ScalarT edgeReynoldsNo = (effPot0 - effPot1) / (2.0*beta);

      // compute edge coefficients
      ScalarT edgeCoef0 = 1.0;      // for pure diffusion
      ScalarT edgeCoef1 = 1.0;
      // if (std::abs(edgeReynoldsNo) > 1.0e-10) // include drift
      double tol = 100.*std::abs(Teuchos::ScalarTraits<double>::eps());
      if (std::abs(Sacado::ScalarValue<ScalarT>::eval(edgeReynoldsNo)) > tol) // include drift
      {
        ScalarT coth = 1.0/tanh(edgeReynoldsNo);
        edgeCoef0 = edgeReynoldsNo * (coth - 1.0); // equivalent to B(2*edgeReynoldsNo)
        edgeCoef1 = edgeReynoldsNo * (coth + 1.0); // equivalent to B(-2*edgeReynoldsNo)
      }

      // compute edge current density (scaled scalar)
      const ScalarT& dens0 = density(cell,node0);
      const ScalarT& dens1 = density(cell,node1);

      // Remove edge length here and in SubCVCurrDens because they end up cancelling - no need to compute
      //edge_currdens(cell,edge) = diff/edgeLen* (dens1*edgeCoef1 -dens0*edgeCoef0);
      edge_currdens(cell,edge) = diff * (dens1*edgeCoef1 -dens0*edgeCoef0);

    }  // end of loop over edges
  }  // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SGCVFEM_EdgeCurrDens<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Carrier Type", "??");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set("Use Electron Quantum Correction", false);
  p->set("Use Hole Quantum Correction", false);

  return p;
}

}

#endif
