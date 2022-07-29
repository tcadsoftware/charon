
#ifndef CHARON_NEUMANNBC_SURFACECHARGE_IMPL_HPP
#define CHARON_NEUMANNBC_SURFACECHARGE_IMPL_HPP

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

const int DEF_NL = 20;

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
NeumannBC_SurfaceCharge<EvalT, Traits>::
NeumannBC_SurfaceCharge(const Teuchos::ParameterList& p)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using std::string;
  using std::vector; 
  using Teuchos::rcp;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // get output data layout
  auto output_dl = p.get< RCP<DataLayout> >("Output Data Layout");
  num_ips = output_dl->dimension(1);

  // obtain basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> input_dl = basis->functional;
  basis_name = basis->name();
  num_nodes = input_dl->dimension(1);
  
  // get charon::Names
  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // scaling parameters
  auto scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  X0 = scaleParams->scale_params.X0;
  R0 = scaleParams->scale_params.R0;
  T0 = scaleParams->scale_params.T0;

  // get flux names
  fluxSurfCharge = p.get<string>("Flux Surface Charge");
  fluxSurfRecomb = p.get<string>("Flux Surface Recombination");

  // get boolean flags and parameterlists
  bFixCharge = p.get<bool>("Include Fixed Charge");
  bSurfTrap = p.get<bool>("Include Surface Trap"); 
  bSurfRecomb = p.get<bool>("Include Surface Recombination"); 
  bPolar = p.get<bool>("Include Polarization");
  bVaryingCharge = p.get<bool>("Include Varying Charge");


  // read in user-specified charge
  fixedCharge = rcp(new panzer::ScalarParameterEntry<EvalT>);
  fixedCharge->setRealValue(0.0);
  if (bFixCharge)  // Fixed Charge is specified
    fixedCharge->setRealValue(p.get<double>("Fixed Charge"));  // in unit of cm^{-2}

  if (bVaryingCharge) //varying charge is specified
    {
      if (p.get<std::string>("Varying Charge") == "Parameter")
	fixedCharge =
	  panzer::createAndRegisterScalarParameter<EvalT>(
							  std::string("Varying Charge"),
							  *p.get<RCP<panzer::ParamLib> >("ParamLib"));
      else
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
				   "BC_NeumannBC_SurfaceCharge():  Error:  Expecting Varying Charge value of " \
				   "\"Parameter\"; received \"" << p.get<std::string>("Varying Charge")
				   << "\".")
  }
  

  // obtain q
  auto const& phyConst = charon::PhysicalConstants::Instance();
  auto q = phyConst.q;   // electron elemental charge in [C]

  if (bPolar)
  {
    const auto& polarPList = 
      *(p.get< RCP<ParameterList> >("Polarization ParameterList"));
    auto type = polarPList.get<string>("Type");
    auto top = polarPList.get<string>("Top");
    auto bottom = polarPList.get<string>("Bottom");
    auto xcomp = polarPList.get<double>("Xcomp");
    auto scale = polarPList.get<double>("Scale");
    auto topMap = getPolarvals(top);
    auto bottomMap = getPolarvals(bottom);
    if (type == "piezo")
      //fixedCharge = scale*piezo(topMap, bottomMap, xcomp)/1.0e4/q;
      fixedCharge->setRealValue(scale*piezo(topMap, bottomMap, xcomp)/1.0e4/q);
    if (type == "total")
      //fixedCharge = scale*polarization(topMap, bottomMap, xcomp)/1.0e4/q;
      fixedCharge->setRealValue(scale*polarization(topMap, bottomMap, xcomp)/1.0e4/q);
  }

  if (bSurfTrap)   // Surface Trap is specified
  { 
    const ParameterList& surfTrapPList = *(p.get< RCP<ParameterList> >("Surface Trap ParameterList"));
    initSurfTrapParams(surfTrapPList);
  }

  if (bSurfRecomb) // Surface Recombination is specified
  { 
    const ParameterList& surfRecombPList = *(p.get< RCP<ParameterList> >("Surface Recombination ParameterList"));
    eSurfRecombVel = surfRecombPList.get<double>("Electron Surface Velocity");
    hSurfRecombVel = surfRecombPList.get<double>("Hole Surface Velocity");
    if (surfRecombPList.isParameter("Energy Level"))
      surfRecombEnergy = surfRecombPList.get<double>("Energy Level"); 
    else
      surfRecombEnergy = 0.0;
  }  

  // evaluated fields
  if (bFixCharge || bVaryingCharge || bSurfTrap || bPolar)
  {
    surface_charge = MDField<ScalarT,Cell,Point>(fluxSurfCharge, output_dl);
    this->addEvaluatedField(surface_charge);
  }
  
  if (bSurfTrap || bSurfRecomb)
  {
    surface_recomb = MDField<ScalarT,Cell,Point>(fluxSurfRecomb, output_dl);
    this->addEvaluatedField(surface_recomb);
  }
  
  // dependent fields
  if (bSurfTrap || bSurfRecomb)
  {
    edensity =  MDField<const ScalarT,Cell,Point>(names.dof.edensity, input_dl);
    hdensity =  MDField<const ScalarT,Cell,Point>(names.dof.hdensity, input_dl);
    intrin_conc = MDField<const ScalarT,Cell,Point>(names.field.intrin_conc, input_dl);
    latt_temp = MDField<const ScalarT,Cell,Point>(names.field.latt_temp, input_dl);
    this->addDependentField(edensity);
    this->addDependentField(hdensity);
    this->addDependentField(intrin_conc);
    this->addDependentField(latt_temp);
  }
  
  string n = "NeumannBC Surface Charge";
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void NeumannBC_SurfaceCharge<EvalT, Traits>::postRegistrationSetup(
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
void NeumannBC_SurfaceCharge<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  double fixedChargeLocal = fixedCharge->getRealValue();
  double scaling4Charge = C0*X0;
  double scaling4Recomb = R0*X0; 
  double scaledFixCharge = fixedChargeLocal / scaling4Charge; 

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  
  // zero out the arrays
  surface_charge.deep_copy(ScalarT(0.0));
  surface_recomb.deep_copy(ScalarT(0.0));

  if (bFixCharge || bPolar)  // Fixed Charge is specified 
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int point = 0; point < num_ips; ++point)
        surface_charge(cell,point) += scaledFixCharge;
  }

  
  if (bSurfTrap)  // Surface Trap is specified, include charge contribution
  {
    int numofTraps = trapDensity.size();
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      // initialization
      Kokkos::DynRankView<ScalarT,PHX::Device> rate_nodes;
      rate_nodes = createDynRankView(edensity.get_static_view(),"rate_nodes",num_nodes);
      Kokkos::deep_copy(rate_nodes,ScalarT(0.0));
      Kokkos::DynRankView<ScalarT,PHX::Device> charge_nodes;
      charge_nodes = createDynRankView(edensity.get_static_view(),"charge_nodes",num_nodes);
      Kokkos::deep_copy(charge_nodes,ScalarT(0.0));
      Kokkos::DynRankView<ScalarT,PHX::Device> rate_ips;
      rate_ips = createDynRankView(edensity.get_static_view(),"rate_ips",num_ips);
      Kokkos::deep_copy(rate_ips,ScalarT(0.0));
      Kokkos::DynRankView<ScalarT,PHX::Device> charge_ips;
      charge_ips = createDynRankView(edensity.get_static_view(),"charge_ips",num_ips);
      Kokkos::deep_copy(charge_ips,ScalarT(0.0));
      
      
      // loop over nodes
      for (int point = 0; point < num_nodes; ++point)
      {
        ScalarT n, p, nie, latt, kbT;
        n = edensity(cell, point) * C0;       // [cm^(-3)]
        p = hdensity(cell, point) * C0; 
        nie = intrin_conc(cell, point) * C0; 
        latt = latt_temp(cell, point) * T0;   // [K]
        kbT = latt * kb;                      // [eV]
        
        if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
            (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
        {
          // loop over the number of traps
          for (int itrap = 0; itrap < numofTraps; ++itrap)
          {
            ScalarT taun = 1.0/eTrapVelocity[itrap] * std::sqrt(300.0/latt); // [s/cm] or [s.eV/cm]
            ScalarT taup = 1.0/hTrapVelocity[itrap] * std::sqrt(300.0/latt); // [s/cm] or [s.eV/cm]
	    
	    if (trapEnDistr[itrap] == "Level") {
	      double Et = trapEnergy[itrap];  // [eV] measured from Ei, + for above Ei, - for below Ei
	      ScalarT n1 = nie * std::exp(Et/kbT);
	      ScalarT p1 = nie * std::exp(-Et/kbT); 
            
	      // compute the surface recombination rate due to one type of surface traps
	      ScalarT numer = n*p - nie*nie;              // [cm^(-6)]
	      ScalarT denom = taup*(n+n1) + taun*(p+p1);  // [cm^(-4).s]
	      ScalarT Rs = numer / denom;                 // [cm^(-2).s^(-1)]
	      Rs *= 1.0/ scaling4Recomb;                  // scaled
	      
	      // sum up the contribution from all specified traps
	      rate_nodes(point) += Rs; 

	      // compute and sum up the trap charge 
	      ScalarT ft = 0.0;  // trap occupation within [0,1]
	      ScalarT Qt = 0.0;  // trap charge (scaled)
	      if (trapType[itrap] == "Acceptor")      // electron capture (close to Ec)
		{ 
		  ft = (taup*n + taun*p1) / denom; 
		  Qt = - trapDensity[itrap] * ft / scaling4Charge;  // carry -e when occupied
		}
	      else if (trapType[itrap] == "Donor")    // hole capture (close to Ev)
		{
		  ft = (taup*n1 + taun*p) / denom; 
		  Qt = trapDensity[itrap] * ft / scaling4Charge;    // carry +e when occupied
		}
	      charge_nodes(point) += Qt; 
	    } else { // for continuous distributions
	      double Nt = trapDensity[itrap]; // [cm-2 eV-1]
	      std::vector<ScalarT> enLevels;
	      std::vector<ScalarT> levNormDensities; // at mid-points
	      discretizeContDistribution(&enLevels, &levNormDensities, trapEnDistr[itrap], 
				       trapEnergy[itrap], trapEnWidth[itrap], trapNL[itrap]);
	      for (size_t k=0; k<enLevels.size()-1; k++) {
		ScalarT norm_dens = levNormDensities[k]; // [1]
		ScalarT en = 0.5*(enLevels[k]+enLevels[k+1]); // [eV]
		ScalarT deltaE = enLevels[k+1] - enLevels[k]; // [eV]
		ScalarT taun_lev = taun/norm_dens; // [s.eV/cm]
		ScalarT taup_lev = taup/norm_dens; // [s.eV/cm]

		ScalarT n1 = nie * std::exp(en/kbT); // [cm^-3]
		ScalarT p1 = nie * std::exp(-en/kbT); // [cm^-3]

		ScalarT numer = n*p - nie*nie; // [cm^(-6)]
		ScalarT denom = taup_lev*(n+n1) + taun_lev*(p+p1);  // [cm^(-4).s.eV]
		ScalarT Rs_lev = (numer / denom) * deltaE;  // [cm^(-2).s^(-1)]
		Rs_lev *= 1.0 / scaling4Recomb;  // scaled

		// accumulate contribution from levels
		rate_nodes(point) += Rs_lev; 

		// compute and sum up the trap charge 
		ScalarT ft_lev = 0.0;  // trap occupation within [0,1]
		ScalarT Qt_lev = 0.0;  // trap charge (scaled)
		if (trapType[itrap] == "Acceptor")      // electron capture (close to Ec)
		{ 
		  ft_lev = (taup_lev*n + taun_lev*p1) / denom; 
		  Qt_lev = -Nt*norm_dens*deltaE*ft_lev / scaling4Charge;  // carry -e when occupied
		}
		else if (trapType[itrap] == "Donor")    // hole capture (close to Ev)
		{
		  ft_lev = (taup_lev*n1 + taun_lev*p) / denom; 
		  Qt_lev =  Nt*norm_dens*deltaE*ft_lev / scaling4Charge;    // carry +e when occupied
		}
		charge_nodes(point) += Qt_lev; 
	      } // cont distr levels 
	    }
          }  // end of loop over traps           
        }  // end of if (n > 0. && p > 0.)

        else
        {
          rate_nodes(point) = 0.0;
          charge_nodes(point) = 0.0; 
        }
      }  // end of loop over nodes
      
      // interpolate values from nodes to ips
      for (int inode = 0; inode < num_nodes; ++inode)
      {
        for (int ip = 0; ip < num_ips; ++ip)
        {
          rate_ips(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * rate_nodes(inode);
          charge_ips(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * charge_nodes(inode);
        }
      }
      
      // add the contribution
      for (int ip = 0; ip < num_ips; ++ip)
      {
        surface_charge(cell,ip) += charge_ips(ip);
        surface_recomb(cell,ip) += rate_ips(ip); 
      }
      
    }  // end of loop over cells
  }  // end of if (bSurfTrap)

  
  if (bSurfRecomb)  // Surface Recombination is specified, neglect charge contribution
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      // initialization
      Kokkos::DynRankView<ScalarT,PHX::Device> rate_nodes;
      rate_nodes = createDynRankView(edensity.get_static_view(),"rate_nodes",num_nodes);
      Kokkos::deep_copy(rate_nodes,ScalarT(0.0));
      Kokkos::DynRankView<ScalarT,PHX::Device> rate_ips;
      rate_ips = createDynRankView(edensity.get_static_view(),"rate_ips",num_ips);
      Kokkos::deep_copy(rate_ips,ScalarT(0.0));

      for (int point = 0; point < num_nodes; ++point)
      {
        ScalarT n, p, nie, latt, kbT;
        n = edensity(cell, point) * C0;       // [cm^(-3)]
        p = hdensity(cell, point) * C0; 
        nie = intrin_conc(cell, point) * C0; 
        latt = latt_temp(cell, point) * T0;   // [K]
        kbT = latt * kb;                      // [eV]

        if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
            (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
        {
          double taun = 1.0/eSurfRecombVel;  // [s/cm]
          double taup = 1.0/hSurfRecombVel;  // [s/cm]
          double Et = surfRecombEnergy;      // [eV] measured from Ei, + for above Ei, - for below Ei
          
          ScalarT n1 = nie * std::exp(Et/kbT);
          ScalarT p1 = nie * std::exp(-Et/kbT); 

          // compute the surface recombination rate 
          ScalarT numer = n*p - nie*nie;              // [cm^(-6)]
          ScalarT denom = taup*(n+n1) + taun*(p+p1);  // [cm^(-4).s]
          ScalarT Rs = numer / denom;                 // [cm^(-2).s^(-1)]
          Rs *= 1.0/ scaling4Recomb;                  // scaled
          rate_nodes(point) = Rs; 

        }  // end of if (n > 0. && p > 0.)
        else
        {
          rate_nodes(point) = 0.0;
        }
      }  // end of loop over nodes
      
      // interpolate values from nodes to ips
      for (int inode = 0; inode < num_nodes; ++inode)
        for (int ip = 0; ip < num_ips; ++ip)
          rate_ips(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * rate_nodes(inode);
      
      // add the contribution
      for (int ip = 0; ip < num_ips; ++ip)
        surface_recomb(cell,ip) += rate_ips(ip); 

    }  // end of loop over cells
  }  // end of if (bSurfRecomb)       
}





template<typename EvalT, typename Traits>
void NeumannBC_SurfaceCharge<EvalT, Traits>::
discretizeContDistribution(std::vector<ScalarT>* enLevels, 
			  std::vector<ScalarT>* norm_densities,
			  const std::string& enDistr, double Et, 
			  double enSigma, int nL) {
  ScalarT en_step = 2.0*enSigma/(nL-1);
  for (int i=0; i<nL; i++) {
    ScalarT en_lev = Et - enSigma + i*en_step; 
    enLevels->push_back(en_lev);
  }
 
  if (enDistr == "Uniform") {  
    for (int i=0; i<nL-1; i++) {
      // normalized densities at energy mid points
      norm_densities->push_back(1.0);
    }
  } else if (enDistr == "Exponential") {
    for (int i=0; i<nL-1; i++) {
      ScalarT en_lev = 0.5*((*enLevels)[i] + (*enLevels)[i+1]);
      // normalized densities at energy mid points
      norm_densities->push_back(
	std::exp(-abs(en_lev-Et)/enSigma) );
    }
  } else if (enDistr == "Gaussian") {
    for (int i=0; i<nL-1; i++) {
      // normalized densities at energy mid points
      ScalarT en_lev = 0.5*((*enLevels)[i] + (*enLevels)[i+1]);
      norm_densities->push_back(
	std::exp(-(en_lev-Et)*(en_lev-Et)/(2.0*enSigma*enSigma)) );
    }
  } 
}





///////////////////////////////////////////////////////////////////////////////
//
//  initSurfTrapParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void NeumannBC_SurfaceCharge<EvalT, Traits>::initSurfTrapParams(const Teuchos::ParameterList& trapPList)
{
  using Teuchos::ParameterList;
  using std::string;

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  double q = phyConst.q;    // elemental charge in [C]
  double m0 = phyConst.m0;  // electron mass in [kg]
  double velPrefactor = std::sqrt(3.0*kb*q*300./m0) * 100.;  // [cm/s], q is to convert kb from eV/K to J/K
  
  double eMass = trapPList.get<double>("Electron Effective Mass");
  double hMass = trapPList.get<double>("Hole Effective Mass"); 
  
  // loop over the number of specified traps
  for (ParameterList::ConstIterator it = trapPList.begin(); it != trapPList.end(); ++it)
  {
    const string key = it->first;
    if (key.find("Trap") != string::npos)  // find the Trap 0, Trap 1, etc. parameterlists
    { 
      const Teuchos::ParameterEntry& entry = it->second;
      const ParameterList& plist = Teuchos::getValue<ParameterList>(entry);

      // trap related parameters (required)
      trapEnergy.push_back(plist.get<double>("Trap Energy"));  // [eV]

      double Nst = plist.get<double>("Trap Density");  // [cm^(-2)] or [cm^(-2)/eV]
      trapDensity.push_back(Nst);

      string tType = plist.get<string>("Trap Type");
      if ((tType != "Acceptor") && (tType != "Donor"))  // either Acceptor or Donor
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! Either Acceptor or Donor must be specified for Trap Type!");
      trapType.push_back(tType);

      string enDistr = "Level";
      if (plist.isParameter("Energy Distribution")) {
	enDistr = plist.get<std::string>("Energy Distribution");
	if ((enDistr != "Level") && (enDistr != "Uniform") && 
	    (enDistr != "Exponential") && (enDistr != "Gaussian")) 
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
            "Error ! Energy Distribution  must be Level, Uniform, Exponential or Gaussian!");
      }
      trapEnDistr.push_back(enDistr);
      if (enDistr != "Level") {
	if (plist.isParameter("Energy Width")) {
	  trapEnWidth.push_back(plist.get<double>("Energy Width"));
	} else {
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	    "Error ! Energy Distribution  must have a width specified!");
	}
      }
      int NL = DEF_NL;
      if (plist.isParameter("Number of Levels")) {
	if (!plist.isParameter("Energy Distribution")) {
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
		"Error ! Number of Levels can be specified only for continuous distributions!");
	} else {
	  if (plist.get<string>("Energy Distribution") == "Level") 
	    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
	     "Error ! Number of Levels can be specified only for continuous distributions!");
	}
	NL = plist.get<int>("Number of Levels");
      }
      trapNL.push_back(NL);

      // compute electron surface recombination velocity due to traps at 300 K
      double eXsec = plist.get<double>("Electron Cross Section");  // [cm^2]
      double eVel = velPrefactor / std::sqrt(eMass);  // [cm/s]
      double eSurfVel = eXsec * eVel * Nst;  // [cm/s] or [cm/s/eV]
      eTrapVelocity.push_back(eSurfVel); 
    
      // compute hole surface recombination velocity due to traps at 300 K
      double hXsec = plist.get<double>("Hole Cross Section");  // [cm^2]
      double hVel = velPrefactor / std::sqrt(hMass);  // [cm/s]
      double hSurfVel = hXsec * hVel * Nst;  // [cm/s] or [cm/s/eV]
      hTrapVelocity.push_back(hSurfVel); 
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
//  getPolarvals()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
auto 
NeumannBC_SurfaceCharge<EvalT, Traits>::
getPolarvals(const std::string& mat) 
      -> std::map<std::string, double>
{
  auto& matProperty = charon::Material_Properties::getInstance();
  std::map<std::string, double> m;
  m["latt"] = matProperty.getPropertyValue(mat, "Lattice Constant");
  m["e33"] = matProperty.getPropertyValue(mat, "Piezoelectric Constant 33");
  m["e31"] = matProperty.getPropertyValue(mat, "Piezoelectric Constant 31");
  m["psp"] = matProperty.getPropertyValue(mat, "Spontaneous Polarization");
  m["c13"] = matProperty.getPropertyValue(mat, "Elastic Constant 13");
  m["c33"] = matProperty.getPropertyValue(mat, "Elastic Constant 33");

  return m;
}

///////////////////////////////////////////////////////////////////////////////
//
//  Piezo()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
auto 
NeumannBC_SurfaceCharge<EvalT, Traits>::
piezo(const std::map<std::string,double>& top, 
      const std::map<std::string,double>& bottom,
      double xcomp) -> double
{
  for (auto to = top.begin(), bo = bottom.begin(); to != top.end();
       to++, bo++)
  {
    interp[to->first] = (to->second - bo->second)*xcomp + bo->second;
  }

  auto latt = 2.0*(bottom.at("latt") - interp.at("latt"))/interp.at("latt");
  auto result = latt*(interp.at("e31") - interp.at("e33")
                *interp.at("c13")/interp.at("c33"));
  return std::fabs(result); // in C/m2
}

///////////////////////////////////////////////////////////////////////////////
//
//  Polarization()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
auto 
NeumannBC_SurfaceCharge<EvalT, Traits>::
polarization(const std::map<std::string,double>& top, 
      const std::map<std::string,double>& bottom,
      double xcomp) -> double
{
  return  (piezo(top, bottom, xcomp)) + (bottom.at("psp") - interp.at("psp"));
  // in C/m2
}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
NeumannBC_SurfaceCharge<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector; 

  RCP<ParameterList> p = rcp(new ParameterList);

  RCP<PHX::DataLayout> dl;
  p->set("Output Data Layout", dl);
  
  RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  RCP<const charon::Names> n;
  p->set("Names", n);

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set<string>("Flux Surface Charge", "?", "Name for the computed scaled surface charge");
  p->set<string>("Flux Surface Recombination", "?", "Name for the computed scaled surface recombination rate");
  p->set<std::string>("Varying Charge", "Parameter");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));

  p->set<bool>("Include Fixed Charge", false, "Check if Fixed Charge is specified");
  p->set<bool>("Include Varying Charge", false, "Check if Varying Charge is specified");
  p->set<bool>("Include Surface Trap", false, "Check if Surface Trap is specified");
  p->set<bool>("Include Surface Recombination", false, "Check if Surface Recombination is specified");
  p->set<double>("Fixed Charge", 0.0, "Fixed Charge in unit of cm^(-2)"); 

  p->set<bool>("Include Polarization", false, "Check if Polarization is specified?");
  RCP<ParameterList> polarPL = rcp(new ParameterList);
  p->set("Polarization ParameterList", polarPL);

  RCP<ParameterList> trapPL = rcp(new ParameterList);
  p->set("Surface Trap ParameterList", trapPL);

  RCP<ParameterList> surfRecPL = rcp(new ParameterList);
  p->set("Surface Recombination ParameterList", surfRecPL);

  return p;
}

}

#endif

