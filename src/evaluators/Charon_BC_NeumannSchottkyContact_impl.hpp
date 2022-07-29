
#ifndef CHARON_BC_NEUMANNSCHOTTKYCONTACT_IMPL_HPP
#define CHARON_BC_NEUMANNSCHOTTKYCONTACT_IMPL_HPP

#include <cmath>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_BC_NeumannSchottkyContact.hpp"
#include "Charon_Util.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace {
  // tunneling current integrand
  template<typename ScalarT> 
  ScalarT tunn_integrand(const ScalarT& c1, const ScalarT& kbT, 
	        	 const ScalarT& barrier, const ScalarT& x, 
		         const ScalarT& V) {
    ScalarT TC = std::exp( -c1 * std::pow(barrier-x,1.5) );
    ScalarT N = std::log( (1.0+exp(-x/kbT))/(1.0+exp((-x+V)/kbT)) );
    ScalarT val = TC * N;
    return val;
  }
  
  template<typename ScalarT> 
  ScalarT romberg_integr(
  		 ScalarT (*f/* function to integrate */)(
                           const ScalarT&, const ScalarT&, 
			   const ScalarT&, const ScalarT&, const ScalarT&), 
  	 	 const ScalarT& /*lower limit*/ a, 
  		 const ScalarT& /*upper limit*/ b, 
  		 const size_t& max_steps, 
  		 const ScalarT& /*desired accuracy*/ acc,
  		 const ScalarT& c1,
		 const ScalarT& kbT,
		 const ScalarT&barrier) {
    //buffers
    std::vector<ScalarT> R1(max_steps), R2(max_steps);
    // Rp is previous row
    ScalarT* Rp = &R1[0];
    // Rc is current row
    ScalarT* Rc = &R2[0]; 
    //step size
    ScalarT h = (b-a); 
    // first trapezoidal step
    Rp[0] = (f(c1,kbT,barrier,a,a) + f(c1,kbT,barrier,b,a))*h*.5; 
      
    for(size_t i = 1; i < max_steps; ++i){
      h /= 2.;
      ScalarT c = 0.;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(size_t j = 1; j <= ep; ++j){
	c += f(c1,kbT,barrier,a+(2*j-1)*h,a);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)
      
      for(size_t j = 1; j <= i; ++j){
	ScalarT n_k = std::pow(4, j);
	//compute R(i,j)
	Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); 
      }
      
      if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc){
	return Rc[i-1];
      }

      //swap Rn and Rc as only the last row is needed
      ScalarT *rt = Rp;
      Rp = Rc;
      Rc = rt;
    }
    //return our best guess
    return Rp[max_steps-1]; 
  }
}


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BC_NeumannSchottkyContact<EvalT, Traits>::
BC_NeumannSchottkyContact(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using std::string;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // get output data layout
  auto scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_ips = scalar->dimension(1);

  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // scaling parameters
  auto scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;
  J0 = scaleParams->scale_params.J0;
  E0 = scaleParams->scale_params.E0;

  // get flux names
  string eflux_name = p.get<string>("Electron Flux Name");
  string hflux_name = p.get<string>("Hole Flux Name");

  cnt_type = (p.get<string>("Contact Type") == "Electron") ? -1 : 1;
  An = p.get<double>("An"); // [A/K^2-cm^2]
  Ap = p.get<double>("Ap"); // [A/K^2-cm^2]
  Wf = p.get<double>("Work Function");

  withBL = false;
  if (p.isParameter("BL_alpha")) {
    withBL = true;
    BL_alpha = p.get<double>("BL_alpha");		
    BL_beta = p.get<double>("BL_beta");
    BL_gamma = p.get<double>("BL_gamma");
  }

  withTunneling = false;
  if (p.isParameter("tun_m")) {
    withTunneling = true;
    tun_m = p.get<double>("tun_m");

    // read in user-specified voltage
    user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
    user_value->setRealValue(0);
    if (p.isType<double>("Voltage")) {
      user_value->setRealValue(p.get<double>("Voltage"));
    } else if (p.isType<std::string>("Varying Voltage")) {        
      if (p.get<std::string>("Varying Voltage") == "Parameter")
	user_value =
	  panzer::createAndRegisterScalarParameter<EvalT>(
	    std::string("Varying Voltage"),
	    *p.get<RCP<panzer::ParamLib> >("ParamLib"));
      else
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
	     "BC_SchottkyContact():  Error:  Expecting Varying Voltage value of " \
	     "\"Parameter\"; received \"" << p.get<std::string>("Varying Voltage")
	     << "\".");
    }		
  }

  // evaluated fields
  {
    eSurfCurrent = MDField<ScalarT,Cell,Point>(eflux_name,scalar);
    this->addEvaluatedField(eSurfCurrent);
    hSurfCurrent = MDField<ScalarT,Cell,Point>(hflux_name,scalar);
    this->addEvaluatedField(hSurfCurrent);
  }
  
  // dependent fields 
  {
    edensity = MDField<const ScalarT,Cell,Point>(names.dof.edensity,scalar);
    hdensity = MDField<const ScalarT,Cell,Point>(names.dof.hdensity,scalar);
    elec_effdos = MDField<const ScalarT,Cell,Point>(names.field.elec_eff_dos,scalar);
    hole_effdos = MDField<const ScalarT,Cell,Point>(names.field.hole_eff_dos,scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(names.field.eff_band_gap,scalar);
    latt_temp = MDField<const ScalarT,Cell,Point>(names.field.latt_temp,scalar);
    effChi = MDField<const ScalarT,Cell,Point>(names.field.eff_affinity,scalar);

    this->addDependentField(edensity);
    this->addDependentField(hdensity);
    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
    this->addDependentField(eff_bandgap);
    this->addDependentField(latt_temp);
    this->addDependentField(effChi);
 
    if (withBL) {
      rel_perm = MDField<ScalarT,Cell,Point>(names.field.rel_perm,scalar);
      EdotNorm = MDField<ScalarT,Cell,Point>(p.get<string>("EdotNorm"),scalar);
      this->addDependentField(rel_perm);
      this->addDependentField(EdotNorm);
    }

    if (withTunneling) {
      if (not withBL) {
	EdotNorm = MDField<ScalarT,Cell,Point>(p.get<string>("EdotNorm"),scalar);
	this->addDependentField(EdotNorm);
      }  
    }
  }

  string n = "BC Neumann at Schottky Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void BC_NeumannSchottkyContact<EvalT, Traits>::postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  //int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
}




///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BC_NeumannSchottkyContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kb = phyConst.kb;  // Boltzmann constant in [eV/K]

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {

    // loop over cell ips
    for (int ip = 0; ip < num_ips; ++ip) {
      ScalarT n, p, Nc, Nv, effEg, latt, kbT, chi;
      n = edensity(cell,ip) * C0;       // [cm^(-3)]
      p = hdensity(cell,ip) * C0;       // [cm^(-3)]
      Nc = elec_effdos(cell,ip) * C0;   // [cm^(-3)]
      Nv = hole_effdos(cell,ip) * C0;   // [cm^(-3)]
      effEg = eff_bandgap(cell,ip);     // [eV]
      latt = latt_temp(cell,ip) * T0;   // [K]
      chi = effChi(cell,ip);            // [eV]
      kbT = latt * kb;                  // [eV]

      // compute barrier 
      ScalarT barrier = Wf - chi;
      ScalarT BL = 0.0;
      if(withBL) {
	// compute barrier lowering
	ScalarT Enorm = -EdotNorm(cell,ip) * E0; // [V/cm] 
	ScalarT eps = rel_perm(cell,ip) * phyConst.eps0; // [C/(Vcm)]
	ScalarT fact = BL_alpha * std::sqrt(phyConst.q/(4*phyConst.pi*eps)); // [V] or [eV]
	if (Enorm > 0.0 and cnt_type < 0) {
	  // electrons
	  BL = fact * std::sqrt(Enorm);
	  BL += BL_beta * std::pow(Enorm, BL_gamma);
	} else if (Enorm < 0.0 and cnt_type > 0) {
	  // holes
	  BL = fact * std::sqrt(-Enorm);
	  BL += BL_beta * std::pow(-Enorm, BL_gamma);
	}
      }
      // apply barrier lowering correction
      if(cnt_type < 0) 
	barrier -= BL; // n-type
      else 
	barrier += BL; // p-type

      // compute equilibrium conc
      ScalarT n0_B = Nc * std::exp(-barrier/kbT);
      ScalarT p0_B = Nv * std::exp((-effEg + barrier)/kbT);

      // compute q*vn and q*vp
      ScalarT qvn = An*latt*latt/Nc;
      ScalarT qvp = Ap*latt*latt/Nv;

      ScalarT ecur = qvn * (n - n0_B) / J0;
      ScalarT hcur = -qvp * (p - p0_B) / J0;

      eSurfCurrent(cell,ip) = ecur;
      hSurfCurrent(cell,ip) = hcur;

      // add tunneling current
      ScalarT tun_curr = 0.0;
      if(withTunneling) {
	ScalarT voltage = user_value->getValue();
	ScalarT c1t = 8.0*phyConst.pi*std::sqrt(2.0*tun_m*phyConst.m0);
	ScalarT Enorm = fabs(EdotNorm(cell,ip) * E0); // [V/cm] 
	ScalarT c1b = 3.0*phyConst.h*phyConst.q*Enorm;
	ScalarT c1 = c1t/c1b; // [J^-1.5];
	c1 = c1 * std::pow(phyConst.q, 1.5); // [ev^-1.5]
	// integration parameters
	const int max_steps = 65; 
	const ScalarT acc = 1.48e-08;
	
	if(cnt_type < 0) { // electron tunneling
	  ScalarT barrier_n = barrier;
	  ScalarT Emin = voltage;
	  ScalarT Emax = barrier_n;
	  if( (Emax - Emin) > 0.0 ) { 
	    tun_curr = romberg_integr<ScalarT>(tunn_integrand,Emin,Emax,
	       max_steps,acc,c1,kbT,barrier_n); // [eV]
	    ScalarT pre_fac = -An*latt/kb; // [A/cm^2-eV]
	    tun_curr *= pre_fac/J0;
	    eSurfCurrent(cell,ip) += tun_curr;
	  }
	} else { // hole tunneling
	  ScalarT barrier_p = effEg-barrier;
	  ScalarT Emin = voltage;
	  ScalarT Emax = barrier_p;
	  if( (Emax - Emin) > 0.0 ) { 
	    tun_curr = romberg_integr<ScalarT>(tunn_integrand,Emin,Emax,
	       max_steps,acc,c1,kbT,barrier_p); // [eV]
	    ScalarT pre_fac = -Ap*latt/kb; // [A/cm^2-eV]
	    tun_curr *= pre_fac/J0;
	    hSurfCurrent(cell,ip) += tun_curr;
	  }
	}
      }
    } // ips

  } // cells 

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BC_NeumannSchottkyContact<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  RCP<ParameterList> p = rcp(new ParameterList);
  
  RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  RCP<const charon::Names> n;
  p->set("Names", n);

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set<string>("Electron Flux Name", "???");
  p->set<string>("Hole Flux Name", "???");
  
  p->set<string>("Contact Type", "???");

  p->set<double>("An", 0.0);
  p->set<double>("Ap", 0.0);
  p->set<double>("Work Function", 0.0); 

  p->set<double>("BL_alpha", 1.0);
  p->set<double>("BL_beta", 0.0);
  p->set<double>("BL_gamma", 1.0);
  p->set<string>("EdotNorm", "???");

  p->set<double>("tun_m", 1.0);

  p->set<double>("Voltage", 0.0);
  p->set<string>("Varying Voltage", "Parameter");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));
  
  return p;
}

}

#endif

