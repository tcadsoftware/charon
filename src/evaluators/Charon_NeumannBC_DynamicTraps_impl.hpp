
#ifndef CHARON_NEUMANNBC_DYNAMICTRAPS_IMPL_HPP
#define CHARON_NEUMANNBC_DYNAMICTRAPS_IMPL_HPP

#include <cmath>
#include <map>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include <iostream>

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
NeumannBC_DynamicTraps<EvalT, Traits>::
NeumannBC_DynamicTraps(const Teuchos::ParameterList& p)
{
  using Teuchos::ParameterList;
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
 
  RCP<ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));
  const string& matName = p.get<string>("Material Name");

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  X0 = scaleParams->scale_params.X0;
  T0 = scaleParams->scale_params.T0;
  C0 = scaleParams->scale_params.C0; 
  R0 = scaleParams->scale_params.R0;
  E0 = scaleParams->scale_params.E0;

  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  num_ips = ip_scalar->dimension(1);
  num_dims = ir->dl_vector->dimension(2);
  int_rule_degree = ir->cubature_degree;

  // get Basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis->functional->dimension(1);
  basis_name = basis->name();

  // get flux names
  fluxDynTrapsCharge = p.get<string>("Flux Dynamic Traps Charge");
  eFluxDynTrapsRecomb = p.get<string>("Flux Dynamic Traps eRecombination");
  hFluxDynTrapsRecomb = p.get<string>("Flux Dynamic Traps hRecombination");

  // initialize trap parameters
  const RCP<ParameterList> dynamicTrapsParamList = 
    p.get< RCP<ParameterList> >("Dynamic Traps ParameterList");
  initDynamicTrapsParams(matName, *dynamicTrapsParamList);
  withField = traps->WithFieldDepXsec();
  
  // determine the data layout
  RCP<DataLayout> output_scalar = ip_scalar;
  RCP<DataLayout> input_scalar = basis_scalar;  
  
  // dependent fields
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);
  latt_temp =  MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  e_gamma = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor,input_scalar);
  h_gamma = MDField<const ScalarT,Cell,Point>(n.field.hole_deg_factor,input_scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,input_scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,input_scalar);
  eff_bandgap =  MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,input_scalar);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);
  this->addDependentField(latt_temp);
  this->addDependentField(e_gamma);
  this->addDependentField(h_gamma);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(eff_bandgap);

  if (withField) {
    elec_field = MDField<ScalarT,Cell,Point>(p.get<string>("EdotNorm"),output_scalar);
    this->addDependentField(elec_field);
  }

  // Evaluated fields
  eFluxRecomb = MDField<ScalarT,Cell,Point>(eFluxDynTrapsRecomb,output_scalar);
  hFluxRecomb = MDField<ScalarT,Cell,Point>(hFluxDynTrapsRecomb,output_scalar);
  fluxCharge = MDField<ScalarT,Cell,Point>(fluxDynTrapsCharge,output_scalar);
  this->addEvaluatedField(eFluxRecomb);
  this->addEvaluatedField(hFluxRecomb);
  this->addEvaluatedField(fluxCharge);
  
  // setup the initial time
  initial_time = 0.0;
  prev_time = initial_time;

  string name = "NeumannBC_DynamicTraps";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void NeumannBC_DynamicTraps<EvalT, Traits>::postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

  // initialize views of solution to be used in trap state computation
  edens = Kokkos::createDynRankView(edensity.get_static_view(),
	"edens", edensity.dimension(0), num_ips);
  hdens = Kokkos::createDynRankView(hdensity.get_static_view(),
	"hdens", hdensity.dimension(0), num_ips);
  lT = Kokkos::createDynRankView(latt_temp.get_static_view(),
	"lT", latt_temp.dimension(0), num_ips);
  egamma = Kokkos::createDynRankView(e_gamma.get_static_view(),
	"egamma", e_gamma.dimension(0), num_ips);
  hgamma = Kokkos::createDynRankView(h_gamma.get_static_view(),
	"hgamma", h_gamma.dimension(0), num_ips);
  e_effdos = Kokkos::createDynRankView(elec_effdos.get_static_view(),
	"e_effdos", elec_effdos.dimension(0), num_ips);
  h_effdos = Kokkos::createDynRankView(hole_effdos.get_static_view(),
	"h_effdos", hole_effdos.dimension(0), num_ips);
  eff_bg =  Kokkos::createDynRankView(eff_bandgap.get_static_view(),
	"eff_bg", eff_bandgap.dimension(0), num_ips);
  if (withField)
    field = Kokkos::createDynRankView(elec_field.get_static_view(),
	    "field", elec_field.dimension(0), num_ips);

  // trap state initialization
  if (withField) 
    traps->initTrapsStateWithField(edens,hdens,lT,egamma,hgamma,e_effdos,h_effdos,eff_bg,field);
  else
    traps->initTrapsState(edens,hdens,lT,egamma,hgamma,e_effdos,h_effdos,eff_bg);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void NeumannBC_DynamicTraps<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;
  using std::vector;
  using Teuchos::RCP;
  using Kokkos::DynRankView;
  using TrapType = typename Trap<EvalT>::Type;
  using EnDistr = typename Trap<EvalT>::EnergyDistribution;

  double curr_time = t0*workset.time; // [s]

  // zero out or initialize the workset slice 
  //of input arrays into be interpolated  
  for (index_t cell = 0; cell < workset.num_cells; ++cell) 
    for (int ip = 0; ip < num_ips; ++ip) {
      edens(cell,ip) = hdens(cell,ip) = 0.0;
      egamma(cell,ip) = hgamma(cell,ip) = 0.0;
      e_effdos(cell,ip) = h_effdos(cell,ip) = 0.0;
      lT(cell,ip) = eff_bg(cell,ip) = 0.0;
      if (withField) 
	field(cell,ip)= elec_field(cell,ip);
    }
  
  // update solution and dependent variables 
  // used to compute traps state
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    // interpolate scalar fields at nodes to ips
    for (int inode = 0; inode < num_nodes; ++inode) {
      for (int ip = 0; ip < num_ips; ++ip) {
	edens(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  edensity(cell,inode);
	hdens(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  hdensity(cell,inode);
	lT(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  latt_temp(cell,inode);
	egamma(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  e_gamma(cell,inode);
	hgamma(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  h_gamma(cell,inode);
	e_effdos(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  elec_effdos(cell,inode);
	h_effdos(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  hole_effdos(cell,inode);
	eff_bg(cell,ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * 
	  eff_bandgap(cell,inode);
      }
    }
  }

  if (curr_time == initial_time) {
    // the initial state is an equilibrium state, therefore
    // a trap is in equilibrium with all its sources or reservoirs
    // of carriers (conduction and valence bands for the time being);
    // the net recombination rate with the carrier sources is zero;
    // for the particular case where we have only conduction and 
    // valence bands as carrier reservoirs, the trap rate of exchange 
    // with one band is equal to the rate of exchange with the other
    // band but with a negative sign 

    // compute traps initial state
    traps->computeTrapsInitialState(initial_time);

    // zero out evaluated fields
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < num_ips; ++point) {
	eFluxRecomb(cell,point) = 0.0;
	hFluxRecomb(cell,point) = 0.0;
	fluxCharge(cell,point) = 0.0;
      }
    }

    for (size_t itrap=0; itrap<traps->GetTrapNo(); itrap++) {
      const RCP<Trap<EvalT>> trap = traps->GetTrap(itrap);
      // retrieve trap features
      ScalarT Nt = trap->GetNt(); // cm^-2
      TrapType type = trap->GetType();
      EnDistr en_dist = trap->GetDistribution();
      double xmin = trap->GetXmin();
      double xmax = trap->GetXmax();
      double ymin = trap->GetYmin();
      double ymax = trap->GetYmax();
      double zmin = trap->GetZmin();
      double zmax = trap->GetZmax();
     
      // retrieve trap state
      DynRankView<ScalarT,PHX::Device> cC;
      DynRankView<ScalarT,PHX::Device> cV;
      DynRankView<ScalarT,PHX::Device> eC;
      DynRankView<ScalarT,PHX::Device> eV;
      DynRankView<ScalarT,PHX::Device> prob;
      DynRankView<ScalarT,PHX::Device> cC_distr;
      DynRankView<ScalarT,PHX::Device> cV_distr;
      DynRankView<ScalarT,PHX::Device> eC_distr;
      DynRankView<ScalarT,PHX::Device> eV_distr;
      DynRankView<ScalarT,PHX::Device> prev_prob_distr;
      DynRankView<ScalarT,PHX::Device> prob_distr;
      if (en_dist == EnDistr::LEVEL) {
	cC = trap->read_cC();
	cV = trap->read_cV();
	eC = trap->read_eC();
	eV = trap->read_eV();
	prob = trap->read_prob();
      } else { // continuous distribution
	cC_distr = trap->read_cC_distr();
        cV_distr = trap->read_cV_distr();
        eC_distr = trap->read_eC_distr();
        eV_distr = trap->read_eV_distr();
	prob_distr = trap->read_prob_distr();
      } // continuous distribution

      // compute rates and trapped charge densities based on traps state
      for (index_t cell = 0; cell < workset.num_cells; ++cell) {
	for (int point = 0; point < num_ips; ++point) {
	  double xcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
	  if (num_dims == 2) {
	    double ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
	    if (xcoord < xmin || xcoord > xmax || ycoord < ymin || ycoord > ymax)
	      continue;
	  }
	  if (num_dims == 3) {
	    double ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
	    double zcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,2);
	    if (xcoord < xmin || xcoord > xmax || ycoord < ymin || ycoord > ymax ||
		zcoord < zmin || zcoord > zmax)
	      continue;
	  }
	  if (en_dist == EnDistr::LEVEL) {
	    if (type == TrapType::ACCEPTOR) {
	      eFluxRecomb(cell,point) += Nt*( 
		(1.0-prob(cell,point)) * cC(cell,point) - 
		prob(cell,point) * eC(cell,point) ) / (R0*X0); 
	      hFluxRecomb(cell,point) += Nt*(
		prob(cell,point) * eV(cell,point) - 
		(1.0-prob(cell,point)) * cV(cell,point) ) / (R0*X0);
	      fluxCharge(cell,point) -= prob(cell,point) * Nt / (C0*X0);
	    } else { // TrapType::DONOR
	      eFluxRecomb(cell,point) += Nt*(
	        prob(cell,point) * eC(cell,point) -
                (1.0-prob(cell,point)) * cC(cell,point) ) / (R0*X0);				    
	      eFluxRecomb(cell,point) += Nt*(
                (1.0-prob(cell,point)) * cV(cell,point) -
                prob(cell,point) * eV(cell,point) ) / (R0*X0);
	      fluxCharge(cell,point) += prob(cell,point) * Nt / (C0*X0);
	    }
	  } else { // continuous distribution
	    for (int k=0; k<trap->GetnL()-1; k++) {
	      ScalarT deltaE = (trap->GetEnLevels())[k+1] - (trap->GetEnLevels())[k];
	      ScalarT dens = Nt*(trap->GetNormDensities())[k];
	      if (type == TrapType::ACCEPTOR) {
		eFluxRecomb(cell,point) += dens*( 
		   (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) - 
		   prob_distr(k,cell,point) * eC_distr(k,cell,point) ) * deltaE / (R0*X0); 
		hFluxRecomb(cell,point) += dens*(
		   prob_distr(k,cell,point) * eV_distr(k,cell,point) - 
		   (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) ) * deltaE  / (R0*X0);
		fluxCharge(cell,point) -= prob_distr(k,cell,point) * dens * deltaE / (C0*X0);
	      } else { // TrapType::DONOR
		eFluxRecomb(cell,point) += dens*(
		   prob_distr(k,cell,point) * eC_distr(k,cell,point) -
		   (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) ) * deltaE / (R0*X0);          
		eFluxRecomb(cell,point) += dens*(
		   (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) -
                   prob_distr(k,cell,point) * eV_distr(k,cell,point) ) * deltaE / (R0*X0);
		fluxCharge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / (C0*X0);
	      }
	    } // levels
	  } // continuous distribution
	} // point
      } // cell
    } // rates and trapped charge

  } else { // compute time-dependent traps states and rates
    // save previous occupation probability and previous time 
    if (curr_time > prev_time) {
      traps->saveTrapsState(prev_time);
    }

    // compute traps state
    traps->computeTrapsState(curr_time);

    // zero out evaluated fields
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < num_ips; ++point) {
	eFluxRecomb(cell,point) = 0.0;
	hFluxRecomb(cell,point) = 0.0;
	fluxCharge(cell,point) = 0.0;
      }
    }

    // compute rates and trapped charge
    for (size_t itrap=0; itrap<traps->GetTrapNo(); itrap++) {
      const RCP<Trap<EvalT>> trap = traps->GetTrap(itrap);
      // retrieve trap features
      ScalarT Nt = trap->GetNt(); // cm^-3
      TrapType type = trap->GetType();
      EnDistr en_dist = trap->GetDistribution();
      double xmin = trap->GetXmin();
      double xmax = trap->GetXmax();
      double ymin = trap->GetYmin();
      double ymax = trap->GetYmax();
      double zmin = trap->GetZmin();
      double zmax = trap->GetZmax();

      // retrieve trap state
      DynRankView<ScalarT,PHX::Device> cC;
      DynRankView<ScalarT,PHX::Device> cV;
      DynRankView<ScalarT,PHX::Device> eC;
      DynRankView<ScalarT,PHX::Device> eV;
      DynRankView<ScalarT,PHX::Device> prob;
      DynRankView<ScalarT,PHX::Device> cC_distr;
      DynRankView<ScalarT,PHX::Device> cV_distr;
      DynRankView<ScalarT,PHX::Device> eC_distr;
      DynRankView<ScalarT,PHX::Device> eV_distr;
      DynRankView<ScalarT,PHX::Device> prev_prob_distr;
      DynRankView<ScalarT,PHX::Device> prob_distr;
      if (en_dist == EnDistr::LEVEL) {
	cC = trap->read_cC();
	cV = trap->read_cV();
	eC = trap->read_eC();
	eV = trap->read_eV();
	prob = trap->read_prob();
      } else { // continuous distribution
	cC_distr = trap->read_cC_distr();
        cV_distr = trap->read_cV_distr();
        eC_distr = trap->read_eC_distr();
        eV_distr = trap->read_eV_distr();
	prob_distr = trap->read_prob_distr();
      } // continuous distribution

      // compute rates and trapped charge densities based on traps state
      for (index_t cell = 0; cell < workset.num_cells; ++cell) {
	for (int point = 0; point < num_ips; ++point) {
	  double xcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
	  if (num_dims == 2) {
	    double ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
	    if (xcoord < xmin || xcoord > xmax || ycoord < ymin || ycoord > ymax)
	      continue;
	  }
	  if (num_dims == 3) {
	    double ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
	    double zcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,2);
	    if (xcoord < xmin || xcoord > xmax || ycoord < ymin || ycoord > ymax ||
		zcoord < zmin || zcoord > zmax)
	      continue;
	  }
	  if (en_dist == EnDistr::LEVEL) {
	    if (type == TrapType::ACCEPTOR) {
	      eFluxRecomb(cell,point) += Nt*( 
		(1.0-prob(cell,point)) * cC(cell,point) - 
		prob(cell,point) * eC(cell,point) ) / (R0*X0); 
	      hFluxRecomb(cell,point) += Nt*(
		prob(cell,point) * eV(cell,point) - 
		(1.0-prob(cell,point)) * cV(cell,point) ) / (R0*X0);
	      fluxCharge(cell,point) -= prob(cell,point) * Nt / (C0*X0);
	    } else { // TrapType::DONOR
	      eFluxRecomb(cell,point) += Nt*(
	        prob(cell,point) * eC(cell,point) -
                (1.0-prob(cell,point)) * cC(cell,point) ) / (R0*X0);				    
	      eFluxRecomb(cell,point) += Nt*(
                (1.0-prob(cell,point)) * cV(cell,point) -
                prob(cell,point) * eV(cell,point) ) / (R0*X0);
	      fluxCharge(cell,point) += prob(cell,point) * Nt / (C0*X0);
	    }
	  } else { // continuous distribution
	    for (int k=0; k<trap->GetnL()-1; k++) {
	      ScalarT deltaE = (trap->GetEnLevels())[k+1] - (trap->GetEnLevels())[k];
	      ScalarT dens = Nt*(trap->GetNormDensities())[k];
	      if (type == TrapType::ACCEPTOR) {
		eFluxRecomb(cell,point) += dens*( 
		   (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) - 
		   prob_distr(k,cell,point) * eC_distr(k,cell,point) ) * deltaE / (R0*X0); 
		hFluxRecomb(cell,point) += dens*(
		   prob_distr(k,cell,point) * eV_distr(k,cell,point) - 
		   (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) ) * deltaE / (R0*X0);
		fluxCharge(cell,point) -= prob_distr(k,cell,point) * dens * deltaE / (C0*X0);
	      } else { // TrapType::DONOR
		eFluxRecomb(cell,point) += dens*(
		   prob_distr(k,cell,point) * eC_distr(k,cell,point) -
		   (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) ) * deltaE/ (R0*X0);	    
		eFluxRecomb(cell,point) += dens*(
		   (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) -
		   prob_distr(k,cell,point) * eV_distr(k,cell,point) ) * deltaE/ (R0*X0);
		fluxCharge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / (C0*X0);
	      }
	    } // level
	  } // continuous distribution
	} // point
      } // cell
    } // rates and trapped charge
  } // compute time-dependent traps states and rates

  prev_time = curr_time; // [s]

}


///////////////////////////////////////////////////////////////////////////////
//
//  initDynamicTrapsParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void NeumannBC_DynamicTraps<EvalT, Traits>::
initDynamicTrapsParams(const std::string& matName,
		       const Teuchos::ParameterList& trapsPL)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  charon::Material_Properties& matProperty = 
    charon::Material_Properties::getInstance();
  double eMass = matProperty.getPropertyValue(matName, "Electron Effective Mass");
  double hMass = matProperty.getPropertyValue(matName, "Hole Effective Mass");

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  double q = phyConst.q;    // elemental charge in [C]
  double m0 = phyConst.m0;  // electron mass in [kg]
  double pi = phyConst.pi;
  double velPrefactor1 = std::sqrt(8.0*kb*q/pi/m0)*100.0;  // [cm/s/K^1/2]
  double velPrefactor2 = std::sqrt(3.0*kb*q/m0)*100.0; // [cm/s/K^1/2]
  double eVelPrefactor = velPrefactor1/std::sqrt(eMass); // [cm/s/K^1/2]
  double hVelPrefactor = velPrefactor1/std::sqrt(hMass); // [cm/s/K^1/2]

  // create traps according to user input / check for validity
  RCP<std::vector<Teuchos::RCP<Trap<EvalT>>>> trap_entries =
    rcp(new std::vector<Teuchos::RCP<Trap<EvalT>>>);
  for (auto itr = trapsPL.begin(); itr != trapsPL.end(); ++itr) {
    const string trap_name = itr->first;
    const Teuchos::ParameterEntry& entry = itr->second;
    const ParameterList& trapPList = Teuchos::getValue<ParameterList>(entry);
    string trap_type = trapPList.get<string>("Trap Type");
    if ((trap_type != "Acceptor") && (trap_type != "Donor"))  // either Acceptor or Donor
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	 "Error ! Either Acceptor or Donor must be specified for Trap Type!");
    // trap energy relative to band edge, for acceptor traps energy from
    // trap to conduction band (positive value), for donor traps from
    // trap to valence band (positive value)
    double Et = trapPList.get<double>("Energy Level"); // [eV]
    // trap density
    double Nt = trapPList.get<double>("Trap Density");
    double g = trapPList.get<double>("Degeneracy Factor");
    double eSigma = trapPList.get<double>("Electron Cross Section");
    if (eSigma <= 0.0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	 "Error ! Electron Cross Section must be strictly positive!");
    double hSigma = trapPList.get<double>("Hole Cross Section");
    if (hSigma <= 0.0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	 "Error ! Hole Cross Section must be strictly positive!");
    string en_dist = "Level";
    double en_sigma = 0.0;
    if (trapPList.isParameter("Energy Distribution")) {
      const string dist = trapPList.get<string>("Energy Distribution");
      if (!(dist == "Uniform") and !(dist == "Exponential") and
	  !(dist == "Gaussian") and  !(dist == "Level")) 
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		"Error ! Energy Distribution must be 'Level','Uniform'" \
                "'Exponential' or 'Gaussian' !");
      if (dist == "Uniform" or dist == "Exponential" or dist == "Gaussian") {
	if (trapPList.isParameter("Energy Width")) {
	  if (trapPList.get<double>("Energy Width") > 0.0) {
	    en_dist = dist;
	    en_sigma = trapPList.get<double>("Energy Width"); 
	  } else {
	    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		"Error ! Energy Width must be strictly positive!");
	  }
	}
      } 
    } 
    if (trapPList.isParameter("Thermal Velocity Calculation")) {
      const string calc_type = trapPList.get<string>("Thermal Velocity Calculation");
      if (!(calc_type == "Mean") and !(calc_type == "Root Mean Square")) 
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		"Error ! Thermal Velocity Calculation must be 'Mean'" \
                "or 'Root Mean Square' !");
      if (calc_type == "Root Mean Square") {
	eVelPrefactor = velPrefactor2/std::sqrt(eMass); // [cm/s/K^1/2]
	hVelPrefactor = velPrefactor2/std::sqrt(hMass); // [cm/s/K^1/2]
      }
    }
    int NL = 20;
    if (trapPList.isParameter("Number of Levels")) {
      if (!trapPList.isParameter("Energy Distribution")) {
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		"Error ! Number of Levels can be specified only for continuous distributions!");
      } else {
	if (trapPList.get<string>("Energy Distribution") == "Level") 
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	     "Error ! Number of Levels can be specified only for continuous distributions!");
      }
      NL = trapPList.get<int>("Number of Levels");
    }
    double e_xdep = 0.0;
    int e_field_dep=0;
    if (trapPList.isParameter("Electron Electric Field Power Dependency")) {
      e_xdep = trapPList.get<double>("Electron Electric Field Power Dependency");
      if (e_xdep <= 0.0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	 "Error !  Electron Electric Field Power Dependency must be strictly positive!");
      if (trapPList.isParameter("Electron Field Dependence")) {
	std::string type = trapPList.get<string>("Electron Field Dependence");
	if (!(type == "Saturation") and !(type == "Kimpton")) 
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		      "Error ! Electron Field Dependence must be 'Saturation'" \
		      "or 'Kimton' !");
	e_field_dep = (type == "Saturation") ? -1 : 1;
      }
    }
    double h_xdep = 0.0;
    int h_field_dep=0;
    if (trapPList.isParameter("Hole Electric Field Power Dependency")) {
      h_xdep = trapPList.get<double>("Hole Electric Field Power Dependency");
      if (h_xdep <= 0.0)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	 "Error !  Hole Electric Field Power Dependency must be strictly positive!");
      if (trapPList.isParameter("Hole Field Dependence")) {
	std::string type = trapPList.get<string>("Hole Field Dependence");
	if (!(type == "Saturation") and !(type == "Kimpton")) 
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		      "Error ! Hole Field Dependence must be 'Saturation'" \
		      "or 'Kimton' !");
	h_field_dep = (type == "Saturation") ? -1 : 1;
      }
    }
    double xmin = -1e20, xmax = 1e20;
    double ymin = -1e20, ymax = 1e20;
    double zmin = -1e20, zmax = 1e20;
    if (trapPList.isParameter("X Min")) xmin = trapPList.get<double>("X Min");
    if (trapPList.isParameter("X Max")) xmax = trapPList.get<double>("X Max");
    if (trapPList.isParameter("Y Min")) ymin = trapPList.get<double>("Y Min");
    if (trapPList.isParameter("Y Max")) ymax = trapPList.get<double>("Y Max");
    if (trapPList.isParameter("Z Min")) zmin = trapPList.get<double>("Z Min");
    if (trapPList.isParameter("Z Max")) zmax = trapPList.get<double>("Z Max");
    
    RCP<Trap<EvalT>> trap = 
      rcp(new Trap<EvalT>(trap_type,"Interface",Et,Nt,g,eSigma,hSigma,
			  en_dist,en_sigma,eVelPrefactor,
			  hVelPrefactor,T0,C0,kb,xmin,xmax,
			  ymin,ymax,zmin,zmax,NL,e_xdep,h_xdep,e_field_dep,h_field_dep));
    trap_entries->push_back(trap);
  }
  traps = rcp(new DynamicTraps<EvalT>(trap_entries));
}



///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
NeumannBC_DynamicTraps<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  //using std::vector; 

  RCP<ParameterList> p = rcp(new ParameterList);
  RCP<const charon::Names> n;
  p->set("Names", n);
  p->set<string>("Material Name", "?");

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<string>("Flux Dynamic Traps Charge", "?", 
		 "Name for the computed scaled dynamic traps charge");
  p->set<string>("Flux Dynamic Traps eRecombination", "?", 
		 "Name for the computed scaled dynamic traps electron rate");
  p->set<string>("Flux Dynamic Traps hRecombination", "?", 
		 "Name for the computed scaled dynamic traps hole rate");
  p->set<string>("EdotNorm", "???");
  
  RCP<ParameterList> trapsPL = rcp(new ParameterList);
  p->set("Dynamic Traps ParameterList", trapsPL);

  return p;
}



}

#endif

