
#ifndef CHARON_RECOMBRATE_DYNAMICTRAPS_IMPL_HPP
#define CHARON_RECOMBRATE_DYNAMICTRAPS_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Material_Properties.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"


const int MAX_NUM_TRAPS = 50;
const int DEF_NL = 20;



namespace charon {


// Trap 
template<typename EvalT>
Trap<EvalT>::Trap(
	  const std::string trap_type, const std::string sp_loc, double en_level, 
	  double dens, double degen, double eSigma, double hSigma, 
	  const std::string en_dist, double en_sigma,
	  double eVelPrefactor, double hVelPrefactor,
          double T0, double C0, double kb, double xmin, 
	  double xmax, double ymin, double ymax, double zmin, 
	  double zmax, int NL, double e_xdep, double h_xdep,
	  int e_field_dep, int h_field_dep) {
  // trap parameters
  if (trap_type == "Donor")  
    type = Type::DONOR;
  else 
    type = Type::ACCEPTOR;
  if (sp_loc == "Volume") 
    loc = Location::VOLUME;
  else 
    loc = Location::INTERFACE;

  Et = en_level; // [eV]
  Nt = dens; // [cm-3], [cm-2], [cm-3 eV-1], [cm-2 eV-1]
  g = degen;
  esigma = eSigma; // [cm^2]
  hsigma = hSigma; // [cm^2]
  if (en_dist == "Level")
    enDistr = EnergyDistribution::LEVEL;
  else if (en_dist == "Uniform")
    enDistr = EnergyDistribution::UNIFORM;
  else if (en_dist == "Exponential")
    enDistr = EnergyDistribution::EXPONENTIAL;
  else 
    enDistr = EnergyDistribution::GAUSSIAN;
  enSigma = en_sigma; // [eV]
  eVel_pre = eVelPrefactor; // [cm/s/K^1/2]
  hVel_pre = hVelPrefactor; // [cm/s/K^1/2]
  T_sc = T0;
  dens_sc = C0;
  kB = kb;

  nL = NL;
  if (enDistr != EnergyDistribution::LEVEL) 
    dicretizeContDistribution();

  esigma_pwr = e_xdep;
  hsigma_pwr = h_xdep;
  if (e_field_dep == -1)
    eFDepXsec = FieldDepXsec::SATURATION;
  else if (e_field_dep == 1)
    eFDepXsec = FieldDepXsec::KIMPTON;
  else 
    eFDepXsec = FieldDepXsec::NONE;
  if (h_field_dep == -1)
    hFDepXsec = FieldDepXsec::SATURATION;
  else if (h_field_dep == 1)
    hFDepXsec = FieldDepXsec::KIMPTON;
  else 
    hFDepXsec = FieldDepXsec::NONE;

  xMin = xmin; xMax = xmax; 
  yMin = ymin; yMax = ymax;
  zMin = zmin; zMax = zmax;
}


// discretize continuous distribution
template<typename EvalT> void 
Trap<EvalT>::dicretizeContDistribution() {
  ScalarT en_step = 2.0*enSigma/(nL-1);
  for (int i=0; i<nL; i++) {
    ScalarT en_lev = Et - enSigma + i*en_step; 
    enLevels.push_back(en_lev);
  }
  if (enDistr == EnergyDistribution::UNIFORM) {  
    for (int i=0; i<nL-1; i++) {
      // normalized densities at energy mid points
      norm_densities.push_back(1.0);
    }
  } else if (enDistr == EnergyDistribution::EXPONENTIAL) {
    for (int i=0; i<nL-1; i++) {
      ScalarT en_lev = 0.5*(enLevels[i] + enLevels[i+1]);
      // normalized densities at energy mid points
      norm_densities.push_back(
	std::exp(-abs(en_lev-Et)/enSigma) );
    }
  } else if (enDistr == EnergyDistribution::GAUSSIAN) {
    for (int i=0; i<nL-1; i++) {
      ScalarT en_lev = 0.5*(enLevels[i] + enLevels[i+1]);
      // normalized densities at energy mid points
      norm_densities.push_back(
	std::exp(-(en_lev-Et)*(en_lev-Et)/(2.0*enSigma*enSigma)) );
    }
  }
}

// initialize trap state
template<typename EvalT> void
Trap<EvalT>::initTrapStateWithField(
     const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
     const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& field) {
  initTrapState(edens,hdens,lT,egamma,hgamma,e_effdos,h_effdos,eff_bg);
  Emag = field;
}

// initialize trap state
template<typename EvalT> void
Trap<EvalT>::initTrapState(
     const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
     const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
     const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg) {
  using std::string;

  n = edens; p = hdens; T=lT; gamma_e = egamma; gamma_h = hgamma;
  Nc = e_effdos; Nv = h_effdos; Eg = eff_bg; 
  // allocate state variables 
  if (enDistr == LEVEL) {
    // for level
    cC = Kokkos::createDynRankView(n, "cC", n.extent(0), n.extent(1));
    cV = Kokkos::createDynRankView(n, "cV", n.extent(0), n.extent(1));
    eC = Kokkos::createDynRankView(n, "eC", n.extent(0), n.extent(1));
    eV = Kokkos::createDynRankView(n, "eV", n.extent(0), n.extent(1));
    prev_prob = Kokkos::createDynRankView(n, "prev_prob", n.extent(0), n.extent(1));
    prob = Kokkos::createDynRankView(n, "prob", n.extent(0), n.extent(1));
  } else {
    // for continuous distribution 
    cC_distr = Kokkos::createDynRankView(n, "cV_distr", nL, n.extent(0), n.extent(1));
    cV_distr = Kokkos::createDynRankView(n, "cV_distr", nL, n.extent(0), n.extent(1));
    eC_distr = Kokkos::createDynRankView(n, "eC_distr", nL, n.extent(0), n.extent(1));
    eV_distr = Kokkos::createDynRankView(n, "eV_distr", nL, n.extent(0), n.extent(1));
    prev_prob_distr = Kokkos::createDynRankView(n, "prev_prob_distr", nL, n.extent(0), n.extent(1));
    prob_distr = Kokkos::createDynRankView(n, "prob_distr", nL, n.extent(0), n.extent(1));
  }
  prev_time = 0.0;
}

// compute trap initial state
template<typename EvalT> void
Trap<EvalT>::computeTrapInitialState(double t0) {
  // compute steady state solution
  for (size_t cell = 0; cell < n.extent(0); ++cell) {
    // compute trap state at IPs
    for (size_t ip = 0; ip < n.extent(1); ++ip) {
      ScalarT vth_n = eVel_pre*std::sqrt(T(cell,ip)*T_sc); // [cm/s]
      ScalarT vth_p = hVel_pre*std::sqrt(T(cell,ip)*T_sc); // [cm/s]
      ScalarT e_sigma = esigma;
      ScalarT h_sigma = hsigma;
      if (esigma_pwr > 0.0) {
	// compute field-dependent electron 
	// capture cross section
	if (eFDepXsec == FieldDepXsec::SATURATION) {
	  if (Emag(cell,ip) > 1.0e6)
	    e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	  else 
	    e_sigma = esigma;
	} else { // Kimpton
	  ScalarT dens = Nt;
	  if (enDistr == EnergyDistribution::UNIFORM) {
	    dens = Nt * 2.0*enSigma; 
	  } else if (enDistr == EnergyDistribution::EXPONENTIAL) {
	    dens = Nt * 2.0*enSigma*(1.0-1.0/std::exp(1.0)); 
	  } else if (enDistr == EnergyDistribution::GAUSSIAN) {
	    // integral of exp(-x*x) = sqrt(M_PI)*erf(x)/2
	    dens = std::sqrt(2.0*M_PI)* Nt * enSigma *
		( std::erf(0.5*std::sqrt(2.0)) - std::erf(0.0) );
	  }
	  if (loc == Location::VOLUME) {
	    ScalarT tmp = std::pow(0.75*std::sqrt(3.141592653589793)/dens, 1.5);
	    if (Emag(cell,ip) <= 1.0e6*std::pow(esigma/tmp,1.0/esigma_pwr) ) {
	      e_sigma = tmp;
	    } else {
	      e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	    }
	  } else { // Location::INTERFACE
	    if (Emag(cell,ip) <= 1.0e6*std::pow(dens*esigma,1.0/esigma_pwr) ) {
	      e_sigma = 1.0/dens;
	    } else {
	      e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	    }
	  }
	}
      }
      if (hsigma_pwr > 0.0) {
	// compute field-dependent hole 
	// capture cross section
	if (hFDepXsec == FieldDepXsec::SATURATION) {
	  if (Emag(cell,ip) > 1.0e6)
	    h_sigma = hsigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	  else
	    h_sigma = hsigma;
	} else { // Kimpton
	  ScalarT dens = Nt;
	  if (enDistr == EnergyDistribution::UNIFORM) {
	    dens = Nt * 2.0*enSigma; 
	  } else if (enDistr == EnergyDistribution::EXPONENTIAL) {
	    dens = Nt * 2.0*enSigma*(1.0-1.0/std::exp(1.0)); 
	  } else if (enDistr == EnergyDistribution::GAUSSIAN) {
	    // integral of exp(-x*x) = sqrt(M_PI)*erf(x)/2
	    dens = std::sqrt(2.0*M_PI)* Nt * enSigma *
		( std::erf(0.5*std::sqrt(2.0)) - std::erf(0.0) );
	  }
	  if (loc == Location::VOLUME) {
	    ScalarT tmp = std::pow(0.75*std::sqrt(3.141592653589793)/dens, 1.5);
	    if (Emag(cell,ip) <= 1.0e6*std::pow(hsigma/tmp,1.0/hsigma_pwr) ) {
	      h_sigma = tmp;
	    } else {
	      h_sigma = hsigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	    }
	  } else { // Location::INTERFACE
	    if (Emag(cell,ip) <= 1.0e6*std::pow(dens*hsigma,1.0/hsigma_pwr) ) {
	      h_sigma = 1.0/dens;
	    } else {
	      h_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	    }
	  }
	}
      }

      if (enDistr == LEVEL) {
	if (type == ACCEPTOR) {
	  ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	    std::exp(-Et/(kB*T(cell,ip)*T_sc)); // scaled
	  ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	    std::exp((Et-Eg(cell,ip))/(kB*T(cell,ip))); // scaled  
	  cC(cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	  cV(cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	  eC(cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	  eV(cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	  ScalarT tmp1 = cC(cell,ip) + cV(cell,ip); // [1/s]
	  ScalarT tmp2 = tmp1 + eC(cell,ip) + eV(cell,ip); // [1/s]
	  prev_prob(cell,ip) = prob(cell,ip) = tmp1/tmp2; //[1]
	} else if (type == DONOR) {
	  ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	    std::exp((Et-Eg(cell,ip))/(kB*T(cell,ip))); // scaled
	  ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	    std::exp(-Et/(kB*T(cell,ip)*T_sc)); // scaled 
	  cC(cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	  cV(cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	  eC(cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	  eV(cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	  ScalarT tmp1 = cC(cell,ip) + cV(cell,ip); // [1/s]
	  ScalarT tmp2 = tmp1 + eC(cell,ip) + eV(cell,ip); // [1/s]
	  prev_prob(cell,ip) = prob(cell,ip) = tmp1/tmp2; //[1]
	}
      } else { // continuous distribution
	// state is computed at energy mid points 
	// between the discrete energy levels
	for (int i=0; i<nL-1; i++) {
	  ScalarT en = 0.5*(enLevels[i]+enLevels[i+1]);
	  if (type == ACCEPTOR) {
	    ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	      std::exp(-en/(kB*T(cell,ip)*T_sc)); // scaled
	    ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	      std::exp((en-Eg(cell,ip))/(kB*T(cell,ip))); // scaled  
	    cC_distr(i,cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	    cV_distr(i,cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	    eC_distr(i,cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	    eV_distr(i,cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	    ScalarT tmp1 = cC_distr(i,cell,ip) + cV_distr(i,cell,ip); // [1/s]
	    ScalarT tmp2 = tmp1 + eC_distr(i,cell,ip) + eV_distr(i,cell,ip); // [1/s]
	    prev_prob_distr(i,cell,ip) = prob_distr(i,cell,ip) = tmp1/tmp2; //[1]
	  } else if (type == DONOR) {
	    ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	      std::exp((en-Eg(cell,ip))/(kB*T(cell,ip))); // scaled
	    ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	      std::exp(-en/(kB*T(cell,ip)*T_sc)); // scaled 
	    cC_distr(i,cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	    cV_distr(i,cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	    eC_distr(i,cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	    eV_distr(i,cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	    ScalarT tmp1 = cC_distr(i,cell,ip) + cV_distr(i,cell,ip); // [1/s]
	    ScalarT tmp2 = tmp1 + eC_distr(i,cell,ip) +  eV_distr(i,cell,ip); // [1/s]
	    prev_prob_distr(i,cell,ip) = prob_distr(i,cell,ip) = tmp1/tmp2; //[1]
	  }
	}
      } // continuous distribution
    }
  }

  // save time for the next state
  prev_time = t0; 
}

// compute trap state
template<typename EvalT> void
Trap<EvalT>::computeTrapState(double t) {
  ScalarT dt = t - prev_time;
  for (size_t cell = 0; cell < n.extent(0); ++cell) {
    // compute trap state at IPs
    for (size_t ip = 0; ip < n.extent(1); ++ip) {
      ScalarT vth_n = eVel_pre*std::sqrt(T(cell,ip)*T_sc); // [cm/s]
      ScalarT vth_p = hVel_pre*std::sqrt(T(cell,ip)*T_sc); // [cm/s]

      ScalarT e_sigma = esigma;
      ScalarT h_sigma = hsigma;
      if (esigma_pwr > 0.0) {
	// compute field-dependent electron 
	// capture cross section
	if (eFDepXsec == FieldDepXsec::SATURATION) {
	  if (Emag(cell,ip) > 1.0e6)
	    e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	  else 
	    e_sigma = esigma;
	} else { // Kimpton
	  ScalarT dens = Nt;
	  if (enDistr == EnergyDistribution::UNIFORM) {
	    dens = Nt * 2.0*enSigma; 
	  } else if (enDistr == EnergyDistribution::EXPONENTIAL) {
	    dens = Nt * 2.0*enSigma*(1.0-1.0/std::exp(1.0)); 
	  } else if (enDistr == EnergyDistribution::GAUSSIAN) {
	    // integral of exp(-x*x) = sqrt(M_PI)*erf(x)/2
	    dens = std::sqrt(2.0*M_PI)* Nt * enSigma *
		( std::erf(0.5*std::sqrt(2.0)) - std::erf(0.0) );
	  }
	  if (loc == Location::VOLUME) {
	    ScalarT tmp = std::pow(0.75*std::sqrt(3.141592653589793)/dens, 1.5);
	    if (Emag(cell,ip) <= 1.0e6*std::pow(esigma/tmp,1.0/esigma_pwr) ) {
	      e_sigma = tmp;
	    } else {
	      e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	    }
	  } else { // Location::INTERFACE
	    if (Emag(cell,ip) <= 1.0e6*std::pow(dens*esigma,1.0/esigma_pwr) ) {
	      e_sigma = 1.0/dens;
	    } else {
	      e_sigma = esigma * std::pow(Emag(cell,ip)/1.0e6,-esigma_pwr);
	    }
	  }
	}
      }
      if (hsigma_pwr > 0.0) {
	// compute field-dependent hole 
	// capture cross section
	if (hFDepXsec == FieldDepXsec::SATURATION) {
	  if (Emag(cell,ip) > 1.0e6)
	    h_sigma = hsigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	  else
	    h_sigma = hsigma;
	} else { // Kimpton
	  ScalarT dens = Nt;
	  if (enDistr == EnergyDistribution::UNIFORM) {
	    dens = Nt * 2.0*enSigma; 
	  } else if (enDistr == EnergyDistribution::EXPONENTIAL) {
	    dens = Nt * 2.0*enSigma*(1.0-1.0/std::exp(1.0)); 
	  } else if (enDistr == EnergyDistribution::GAUSSIAN) {
	    // integral of exp(-x*x) = sqrt(M_PI)*erf(x)/2
	    dens = std::sqrt(2.0*M_PI)* Nt * enSigma *
		( std::erf(0.5*std::sqrt(2.0)) - std::erf(0.0) );
	  }
	  if (loc == Location::VOLUME) {
	    ScalarT tmp = std::pow(0.75*std::sqrt(3.141592653589793)/dens, 1.5);
	    if (Emag(cell,ip) <= 1.0e6*std::pow(hsigma/tmp,1.0/hsigma_pwr) ) {
	      h_sigma = tmp;
	    } else {
	      h_sigma = hsigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	    }
	  } else { // Location::INTERFACE
	    if (Emag(cell,ip) <= 1.0e6*std::pow(dens*esigma,1.0/hsigma_pwr) ) {
	      h_sigma = 1.0/dens;
	    } else {
	      h_sigma = hsigma * std::pow(Emag(cell,ip)/1.0e6,-hsigma_pwr);
	    }
	  }
	}
      }

      if (enDistr == LEVEL) {
	if (type == ACCEPTOR) {
	  ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	    std::exp(-Et/(kB*T(cell,ip)*T_sc)); // scaled
	  ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	    std::exp((Et-Eg(cell,ip))/(kB*T(cell,ip))); // scaled  
	  cC(cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	  cV(cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	  eC(cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	  eV(cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	  ScalarT tmp1 = cC(cell,ip) + cV(cell,ip); // [1/s]
	  ScalarT tmp2 = tmp1 + eC(cell,ip) + eV(cell,ip); // [1/s]
	  prob(cell,ip) = prev_prob(cell,ip)/(1.0 + dt*tmp2) + 
	    dt*tmp1/(1.0 + dt*tmp2); //[1]
	} else if (type == DONOR) {
	  ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	    std::exp((Et-Eg(cell,ip))/(kB*T(cell,ip))); // scaled
	  ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	    std::exp(-Et/(kB*T(cell,ip)*T_sc)); // scaled 
	  cC(cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	  cV(cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	  eC(cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	  eV(cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	  ScalarT tmp1 = cC(cell,ip) + cV(cell,ip); // [1/s]
	  ScalarT tmp2 = tmp1 + eC(cell,ip) + eV(cell,ip); // [1/s]
	  prob(cell,ip) = prev_prob(cell,ip)/(1.0 + dt*tmp2) + 
	    dt*tmp1/(1.0 + dt*tmp2); //[1]
	}
      } else { // continuous distribution
	// state is computed at energy mid points 
	// between the discrete energy levels
	for (int i=0; i<nL-1; i++) {
	  ScalarT en = 0.5*(enLevels[i]+enLevels[i+1]);
	  if (type == ACCEPTOR) {
	    ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	      std::exp(-en/(kB*T(cell,ip)*T_sc)); // scaled
	    ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	      std::exp((en-Eg(cell,ip))/(kB*T(cell,ip))); // scaled  
	    cC_distr(i,cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	    cV_distr(i,cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	    eC_distr(i,cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	    eV_distr(i,cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	    ScalarT tmp1 = cC_distr(i,cell,ip) + cV_distr(i,cell,ip); // [1/s]
	    ScalarT tmp2 = tmp1 + eC_distr(i,cell,ip) + eV_distr(i,cell,ip); // [1/s]
	    prob_distr(i,cell,ip) = prev_prob_distr(i,cell,ip)/(1.0 + dt*tmp2) + 
	      dt*tmp1/(1.0 + dt*tmp2); //[1]
	  } else if (type == DONOR) {
	    ScalarT n1 = gamma_e(cell,ip)*Nc(cell,ip)*
	      std::exp((enLevels[i]-Eg(cell,ip))/(kB*T(cell,ip))); // scaled
	    ScalarT p1 = gamma_h(cell,ip)*Nv(cell,ip)*
	      std::exp(-enLevels[i]/(kB*T(cell,ip)*T_sc)); // scaled 
	    cC_distr(i,cell,ip) = e_sigma*vth_n*gamma_e(cell,ip)*n1*dens_sc/g; // [1/s]
	    cV_distr(i,cell,ip) = h_sigma*vth_p*p(cell,ip)*dens_sc; // [1/s]
	    eC_distr(i,cell,ip) = e_sigma*vth_n*n(cell,ip)*dens_sc; // [1/s]
	    eV_distr(i,cell,ip) = g*h_sigma*vth_p*gamma_h(cell,ip)*p1*dens_sc; // [1/s]
	    ScalarT tmp1 = cC_distr(i,cell,ip) + cV_distr(i,cell,ip); // [1/s]
	    ScalarT tmp2 = tmp1 + eC_distr(i,cell,ip) + eV_distr(i,cell,ip); // [1/s]
	    prob_distr(i,cell,ip) = prev_prob_distr(i,cell,ip)/(1.0 + dt*tmp2) + 
	      dt*tmp1/(1.0 + dt*tmp2); //[1]
	  }
	}
      } // continuous distribution
    }
  }
}

template<typename EvalT> void
Trap<EvalT>::saveTrapState(double t) {
  // save occupation probability as previous occupation probability
  // and the previous time 
  for (size_t cell = 0; cell < n.extent(0); ++cell) 
    for (size_t ip = 0; ip < n.extent(1); ++ip) {
      if (enDistr == LEVEL) 
	prev_prob(cell,ip) = prob(cell,ip);
      else // continuous distribution
	// save state at energy mid points
	for (int i=0; i<nL-1; i++) 
	  prev_prob_distr(i,cell,ip) = prob_distr(i,cell,ip);
    }
  // save time for the next state
  prev_time = t;
}


// Trap container/manager
template<typename EvalT>
DynamicTraps<EvalT>::DynamicTraps(
    Teuchos::RCP<std::vector<Teuchos::RCP<Trap<EvalT>>>> traps) {
  trap_entries = traps;
}

// initialize traps state

template<typename EvalT> void 
DynamicTraps<EvalT>::initTrapsStateWithField(
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& field) const {
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) {
    (*trap_entries)[itrap]->initTrapStateWithField(edens,hdens,lT,egamma,hgamma,
					e_effdos,h_effdos,eff_bg,field);
  }
}


template<typename EvalT> void 
DynamicTraps<EvalT>::initTrapsState(
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
	   const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg) const {
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) {
    (*trap_entries)[itrap]->initTrapState(edens,hdens,lT,egamma,hgamma,
					e_effdos,h_effdos,eff_bg);
  }
}


// check if traps need field to cross sections
template<typename EvalT> bool 
DynamicTraps<EvalT>::WithFieldDepXsec() {
  bool val = false;
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) 
    if ( (*trap_entries)[itrap]->WithFieldDepXsec() ) {
      val = true;
      break;
    }
  return val;
}


// compute traps initial state
template<typename EvalT> void
DynamicTraps<EvalT>::computeTrapsInitialState(double t0) {
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) 
    (*trap_entries)[itrap]->computeTrapInitialState(t0);
}

// compute traps state
template<typename EvalT> void
DynamicTraps<EvalT>::computeTrapsState(double t) {
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) 
    (*trap_entries)[itrap]->computeTrapState(t);
}

// save traps state
template<typename EvalT> void
DynamicTraps<EvalT>::saveTrapsState(double t) {
  for(size_t itrap=0; itrap<trap_entries->size(); itrap++) 
    (*trap_entries)[itrap]->saveTrapState(t);
}

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_DynamicTraps<EvalT, Traits>::
RecombRate_DynamicTraps(
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

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));
  const string& matName = p.get<string>("Material Name");
  const string& eqnSetType = p.get<string>("Equation Set Type");
  driveForce = p.get<string>("Driving Force");
 
  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  T0 = scaleParams->scale_params.T0;
  E0 = scaleParams->scale_params.E0;
  C0 = scaleParams->scale_params.C0;
  R0 = scaleParams->scale_params.R0;

  // get IR - It is FEM IP for SUPG-FEM and EFFPG-FEM, but is the SubCV centroid for SGCVFEM
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  num_points = ip_scalar->dimension(1);
  num_dims = ip_vector->dimension(2);
  int_rule_degree = ir->cubature_degree;

  // get Basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis->functional->dimension(1);
  basis_name = basis->name();

  // initialize trap parameters
  const Teuchos::ParameterList& dynamicTrapsParamList = p.sublist("Dynamic Traps ParameterList");
  initDynamicTrapsParams(matName, dynamicTrapsParamList);
  withField = traps->WithFieldDepXsec();

  // determine the data layout
  RCP<DataLayout> output_scalar = ip_scalar;
  RCP<DataLayout> input_vector = ip_vector;  // Driving force is always located at IP
  RCP<DataLayout> input_scalar = ip_scalar;  // True for FEM-based discretization, not for SGCVFEM
  isSGCVFEM = false;  // default
  if (eqnSetType == "SGCVFEM Drift Diffusion")
  {
    input_scalar = basis_scalar;
    isSGCVFEM = true;
  }

  // dependent fields
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);
  latt_temp =  MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  e_gamma = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor, input_scalar);
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
    if (isSGCVFEM) 
      elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_negpot,input_vector);
    else 
      elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi,input_vector);
    this->addDependentField(elec_field);
  }
  
  // Evaluated fields
  erecomb_rate = MDField<ScalarT,Cell,Point>(n.field.dynamic_traps_erecomb,output_scalar);
  hrecomb_rate = MDField<ScalarT,Cell,Point>(n.field.dynamic_traps_hrecomb,output_scalar);
  etrapped_charge = MDField<ScalarT,Cell,Point>(n.field.etrapped_charge,output_scalar);
  htrapped_charge = MDField<ScalarT,Cell,Point>(n.field.htrapped_charge,output_scalar);
  trapped_charge = MDField<ScalarT,Cell,Point>(n.field.trapped_charge,output_scalar);
  this->addEvaluatedField(erecomb_rate);
  this->addEvaluatedField(hrecomb_rate);
  this->addEvaluatedField(etrapped_charge);
  this->addEvaluatedField(htrapped_charge);
  this->addEvaluatedField(trapped_charge);

  // setup the initial time
  initial_time = 0.0;
  prev_time = initial_time;

  std::string name = "DynamicTraps_Recombination";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_DynamicTraps<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  using Teuchos::rcp;

  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

  // initialize views of solution to be used in trap state computation
  edens = Kokkos::createDynRankView(edensity.get_static_view(),
	"edens", edensity.dimension(0), num_points);
  hdens = Kokkos::createDynRankView(hdensity.get_static_view(),
	"hdens", hdensity.dimension(0), num_points);
  lT = Kokkos::createDynRankView(latt_temp.get_static_view(),
	"lT", latt_temp.dimension(0), num_points);
  egamma = Kokkos::createDynRankView(e_gamma.get_static_view(),
	"egamma", e_gamma.dimension(0), num_points);
  hgamma = Kokkos::createDynRankView(h_gamma.get_static_view(),
	"hgamma", h_gamma.dimension(0), num_points);
  e_effdos = Kokkos::createDynRankView(elec_effdos.get_static_view(),
	"e_effdos", elec_effdos.dimension(0), num_points);
  h_effdos = Kokkos::createDynRankView(hole_effdos.get_static_view(),
	"h_effdos", hole_effdos.dimension(0), num_points);
  eff_bg = Kokkos::createDynRankView(eff_bandgap.get_static_view(),
	"eff_bg", eff_bandgap.dimension(0), num_points);
  if (withField)
    field = Kokkos::createDynRankView(elec_field.get_static_view(),
	    "field", elec_field.dimension(0), num_points);

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
void
RecombRate_DynamicTraps<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Kokkos::DynRankView;
  using std::string;
  using std::vector;
  using panzer::index_t;
  using TrapType = typename Trap<EvalT>::Type;
  using EnDistr = typename Trap<EvalT>::EnergyDistribution;

  double curr_time = t0*workset.time; // [s]
  
  // zero out or initialize the workset slice of input 
  // arrays to be interpolated for SGCVFEM  
  for (index_t cell = 0; cell < workset.num_cells; ++cell) 
    if (isSGCVFEM) {
      for (int ip = 0; ip < num_points; ++ip) {
	edens(cell,ip) = hdens(cell,ip) = 0.0;
	egamma(cell,ip) = hgamma(cell,ip) = 0.0;
	e_effdos(cell,ip) = h_effdos(cell,ip) = 0.0;
	lT(cell,ip) = eff_bg(cell,ip) = 0.0;
	if (withField) {
	  ScalarT E_mag = 0.0;
	  for (int dim = 0; dim < num_dims; ++dim) 
	    E_mag += elec_field(cell,ip,dim) * elec_field(cell,ip,dim);
	  E_mag = std::sqrt(E_mag) * E0; // [V/cm]
	  field(cell,ip) = E_mag; // [V/cm]
	}
      }
    }

  // update solution and dependent variables 
  // used to compute traps state
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    if (isSGCVFEM) {
      // interpolate scalar fields at nodes to the 
      // centroids of subcv for SGCVFEM
      const int num_ips = num_points;
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
    } else {
      // use values at the IPs of finite elements for SUPG-FEM and EFFPG-FEM
      for (int point = 0; point < num_points; ++point) {
	edens(cell,point) = edensity(cell,point);
	hdens(cell,point) = hdensity(cell,point);
	lT(cell,point) = latt_temp(cell,point);
	egamma(cell,point) = e_gamma(cell,point);
	hgamma(cell,point) = h_gamma(cell,point);
	e_effdos(cell,point) = elec_effdos(cell,point);
	h_effdos(cell,point) = hole_effdos(cell,point);
	eff_bg(cell,point) = eff_bandgap(cell,point);
	if (withField) {
	  ScalarT E_mag = 0.0;
	  for (int dim = 0; dim < num_dims; ++dim) 
	    E_mag += elec_field(cell,point,dim) * elec_field(cell,point,dim);
	  E_mag = std::sqrt(E_mag) * E0; // [V/cm]
	  field(cell,point) = E_mag; // [V/cm]
	}
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
      for (int point = 0; point < num_points; ++point) {
	erecomb_rate(cell,point) = 0.0;
	hrecomb_rate(cell,point) = 0.0;
	etrapped_charge(cell,point) = 0.0;
	htrapped_charge(cell,point) = 0.0;
	trapped_charge(cell,point) = 0.0;
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
	for (int point = 0; point < num_points; ++point) {
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
	      erecomb_rate(cell,point) += Nt*( 
		(1.0-prob(cell,point)) * cC(cell,point) - 
		prob(cell,point) * eC(cell,point) ) / R0; 
	      hrecomb_rate(cell,point) += Nt*(
		prob(cell,point) * eV(cell,point) - 
		(1.0-prob(cell,point)) * cV(cell,point) ) / R0;
	      etrapped_charge(cell,point) += prob(cell,point) * Nt / C0;
	      trapped_charge(cell,point) -= etrapped_charge(cell,point);
	    } else { // TrapType::DONOR
	      erecomb_rate(cell,point) += Nt*(
	        prob(cell,point) * eC(cell,point) -
                (1.0-prob(cell,point)) * cC(cell,point) ) / R0;				    
	      hrecomb_rate(cell,point) += Nt*(
                (1.0-prob(cell,point)) * cV(cell,point) -
                prob(cell,point) * eV(cell,point) ) / R0;
	      htrapped_charge(cell,point) += prob(cell,point) * Nt / C0;
	      trapped_charge(cell,point) += htrapped_charge(cell,point);
	    } 
	  } else { // continuous distribution
	    for (int k=0; k<trap->GetnL()-1; k++) {
	      ScalarT deltaE = (trap->GetEnLevels())[k+1] - (trap->GetEnLevels())[k];
	      ScalarT dens = Nt*(trap->GetNormDensities())[k];
	      if (type == TrapType::ACCEPTOR) {
		erecomb_rate(cell,point) += dens*( 
		  (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) - 
		  prob_distr(k,cell,point) * eC_distr(k,cell,point) ) * deltaE / R0; 
		hrecomb_rate(cell,point) += dens*(
		  prob_distr(k,cell,point) * eV_distr(k,cell,point) - 
		  (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) ) * deltaE/ R0;
		etrapped_charge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / C0;
		trapped_charge(cell,point) -= etrapped_charge(cell,point);
	      } else { // TrapType::DONOR
		erecomb_rate(cell,point) += dens*(
		  prob_distr(k,cell,point) * eC_distr(k,cell,point) -
                  (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) ) * deltaE / R0;   
		hrecomb_rate(cell,point) += dens*(
		  (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) -
		  prob_distr(k,cell,point) * eV_distr(k,cell,point) ) * deltaE / R0;
		htrapped_charge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / C0;
		trapped_charge(cell,point) += htrapped_charge(cell,point);
	      }
	    } // levels
	  } // continuous distribution

	} // point
      } // cell
    } // rates and trapped charge
   
  } else { // compute time-dependent traps states and rates
    // save previous occupation probability and previous
    // time when time step converged
    if (curr_time > prev_time) {
      traps->saveTrapsState(prev_time);
      //std::cout << "Save traps state at t=" << prev_time << std::endl;
    }

    // compute traps state
    traps->computeTrapsState(curr_time);

    // zero out evaluated fields
    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (int point = 0; point < num_points; ++point) {
	erecomb_rate(cell,point) = 0.0;
	hrecomb_rate(cell,point) = 0.0;
	etrapped_charge(cell,point) = 0.0;
	htrapped_charge(cell,point) = 0.0;
	trapped_charge(cell,point) = 0.0;
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
	for (int point = 0; point < num_points; ++point) {
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
	      erecomb_rate(cell,point) += Nt*( 
		(1.0-prob(cell,point)) * cC(cell,point) - 
		prob(cell,point) * eC(cell,point) ) / R0;
	      hrecomb_rate(cell,point) += Nt*(
		prob(cell,point) * eV(cell,point) - 
		(1.0-prob(cell,point)) * cV(cell,point) ) / R0;
	      etrapped_charge(cell,point) += prob(cell,point) * Nt / C0;
	      trapped_charge(cell,point) -= etrapped_charge(cell,point);
	    } else { // TrapType::DONOR
	      erecomb_rate(cell,point) += Nt*(
	        prob(cell,point) * eC(cell,point) -
                (1.0-prob(cell,point)) * cC(cell,point) ) / R0;				    
	      hrecomb_rate(cell,point) += Nt*(
                (1.0-prob(cell,point)) * cV(cell,point) -
                prob(cell,point) * eV(cell,point) ) / R0;
	      htrapped_charge(cell,point) += prob(cell,point) * Nt / C0;
	      trapped_charge(cell,point) += htrapped_charge(cell,point);
	    }
	  } else { // continuous distribution
	    for (int k=0; k<trap->GetnL()-1; k++) {
	      ScalarT deltaE = (trap->GetEnLevels())[k+1] - (trap->GetEnLevels())[k];
	      ScalarT dens = Nt*(trap->GetNormDensities())[k];
	      if (type == TrapType::ACCEPTOR) {
		erecomb_rate(cell,point) += dens*( 
		  (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) - 
		  prob_distr(k,cell,point) * eC_distr(k,cell,point) ) * deltaE / R0;
		hrecomb_rate(cell,point) += dens*(
		  prob_distr(k,cell,point) * eV_distr(k,cell,point) - 
		  (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) ) * deltaE / R0;
		etrapped_charge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / C0;
		trapped_charge(cell,point) -= etrapped_charge(cell,point);
	      } else { // TrapType::DONOR
		erecomb_rate(cell,point) += dens*(
		  prob_distr(k,cell,point) * eC_distr(k,cell,point) -
		  (1.0-prob_distr(k,cell,point)) * cC_distr(k,cell,point) ) * deltaE/ R0;	 
		hrecomb_rate(cell,point) += dens*(
		  (1.0-prob_distr(k,cell,point)) * cV_distr(k,cell,point) -
		  prob_distr(k,cell,point) * eV_distr(k,cell,point) ) * deltaE / R0;
		htrapped_charge(cell,point) += prob_distr(k,cell,point) * dens * deltaE / C0;
		trapped_charge(cell,point) += htrapped_charge(cell,point);
	      }
	    } // level
	  } // continuous distribution
	} // point
      } // cell
    } // rates and trapped charge
  } // compute time-dependent traps states and rates

  // save current time
  prev_time = curr_time;

}


///////////////////////////////////////////////////////////////////////////////
//
//  initDynamicTrapsParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void RecombRate_DynamicTraps<EvalT, Traits>::initDynamicTrapsParams(
const std::string& matName, const Teuchos::ParameterList& trapsPL) {
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
    int NL = DEF_NL;
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
      rcp(new Trap<EvalT>(trap_type,"Volume",Et,Nt,g,eSigma,hSigma,
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
RecombRate_DynamicTraps<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  RCP<const charon::Names> n;
  p->set("Names", n);
  p->set<string>("Material Name", "?");
  p->set<string>("Equation Set Type", "?");
  p->set<string>("Driving Force", "?");

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  ParameterList& trapsPL = p->sublist("Dynamic Traps ParameterList", false, 
			    "Sublist defining the Dynamic Traps ParameterList");
 
  for (int i = 0; i < MAX_NUM_TRAPS; i++)
  {
    std::stringstream ss;
    ss << i;
    string subListName("Trap " + ss.str());
    trapsPL.sublist(subListName, false, "Sublist defining the parameters for an individual trap");

    // check trap parameters
    trapsPL.sublist(subListName).set<string>("Trap Type", "", "Either Acceptor or Donor");
    trapsPL.sublist(subListName).set<double>("Energy Level", 0.0, 
		  "Trap energy level measured from the band edge in [eV]");
    trapsPL.sublist(subListName).set<double>("Trap Density", 0.0, 
		  "Trap density in [cm^-3]");
    trapsPL.sublist(subListName).set<double>("Degeneracy Factor", 0.0, 
		  "Degeneracy factor unitless");
    trapsPL.sublist(subListName).set<double>("Electron Cross Section", 0.0, 
		  "Electron capture cross section in [cm^2]");
    trapsPL.sublist(subListName).set<double>("Hole Cross Section", 0.0, 
		  "Hole capture cross section in [cm^2]");
    trapsPL.sublist(subListName).set<double>("Electron Electric Field Power Dependency",0.0,
		  "Electric Field Power Dependency Exponenent for Electron Capture Cross Section");
    trapsPL.sublist(subListName).set<double>("Hole Electric Field Power Dependency",0.0,
		  "Electric Field Power Dependency Exponenent for Hole Capture Cross Section");
    trapsPL.sublist(subListName).set<string>("Electron Field Dependence", "", 
		  "Electron Capture Cross Section field dependency type");
    trapsPL.sublist(subListName).set<string>("Hole Field Dependence", "", 
		  "Hole Capture Cross Section field dependency type");
    trapsPL.sublist(subListName).set<string>("Energy Distribution", "", "Energy distribution type");
    trapsPL.sublist(subListName).set<double>("Energy Width", 0.0, 
                  "Distribution energy width [eV]");
    trapsPL.sublist(subListName).set<double>("X Min", 0.0, 
                  "X min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("X Max", 0.0, 
                  "X max for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Y Min", 0.0, 
                  "Y min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Y Max", 0.0, 
                  "Y max for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Z Min", 0.0, 
                  "Z min for the spatial box containing trap");
    trapsPL.sublist(subListName).set<double>("Z Max", 0.0, 
                  "Z max for the spatial box containing trap");
    trapsPL.sublist(subListName).set<string>("Thermal Velocity Calculation", "Mean", "Either Mean or Root Mean Square");
    trapsPL.sublist(subListName).set<int>("Number of Levels", 20, 
		  "Number of discrete energy levels for a continuous distribution");
  }

  return p;
}


}

#endif
