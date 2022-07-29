
#ifndef CHARON_RECOMBRATE_DYNAMICTRAPS_DECL_HPP
#define CHARON_RECOMBRATE_DYNAMICTRAPS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;


namespace charon {

// Trap 
template<typename EvalT>
class Trap
{
public:
  using ScalarT = typename EvalT::ScalarT;

  enum Type {
    ACCEPTOR,
    DONOR
  };
  
  enum Location {
    VOLUME,
    INTERFACE
  };

  enum EnergyDistribution {
    LEVEL,
    UNIFORM,
    EXPONENTIAL,
    GAUSSIAN
  };

  enum FieldDepXsec {
    NONE, 
    SATURATION,
    KIMPTON
  };


  Trap(std::string trap_type, std::string sp_loc, double en_level, double dens, 
       double degen, double eSigma, double hSigma, 
       const std::string en_dist, double en_sigma,
       double eVelPrefactor,double hVelPrefactor, double T0, 
       double C0, double kb, double xmin, double xmax, double ymin,
       double ymax, double zmin, double zmax, int NL, double e_xdep, double h_xdep,
       int e_field_dep, int h_field_dep); 

  Trap<EvalT>::Type GetType() { return type; }

  Trap<EvalT>::EnergyDistribution GetDistribution() { return enDistr; } 

  ScalarT GetNt() { return Nt; }

  int GetnL() { return (enDistr == EnergyDistribution::LEVEL) ? 1 : nL; }
  std::vector<ScalarT>& GetEnLevels() { return enLevels; }

  std::vector<ScalarT>& GetNormDensities() { return norm_densities; }

  // check if trap has field-dependent cross sections
  bool WithFieldDepXsec() { 
    if (esigma_pwr > 0.0 or hsigma_pwr > 0.0) return true;
    else return false;
  }

  double GetXmin () { return xMin; }
  double GetXmax () { return xMax; }
  double GetYmin () { return yMin; }
  double GetYmax () { return yMax; }
  double GetZmin () { return zMin; }
  double GetZmax () { return zMax; }
   
  void initTrapStateWithField(const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& field);

  void initTrapState(const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
		     const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg);
  
  void computeTrapInitialState(double t0);

  void computeTrapState(double t);

  void saveTrapState(double t);

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_cC() {
    return cC;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_cV() {
    return cV;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_eC() {
    return eC;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_eV() {
    return eV;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_prob() {
    return prob;
  }

  /*
  const std::vector<Kokkos::DynRankView<ScalarT,PHX::Device>>& read_cC_distr() { 
    return cC_distr;
  }

  const std::vector<Kokkos::DynRankView<ScalarT,PHX::Device>>& read_cV_distr() { 
    return cV_distr;
  }
  
  const std::vector<Kokkos::DynRankView<ScalarT,PHX::Device>>& read_eC_distr() { 
    return eC_distr;
  }

  const std::vector<Kokkos::DynRankView<ScalarT,PHX::Device>>& read_eV_distr() { 
    return eV_distr;
  }

  const std::vector<Kokkos::DynRankView<ScalarT,PHX::Device>>& read_prob_distr() { 
    return prob_distr;
  }
  */
  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_cC_distr() { 
    return cC_distr;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_cV_distr() { 
    return cV_distr;
  }
  
  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_eC_distr() { 
    return eC_distr;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_eV_distr() { 
    return eV_distr;
  }

  const Kokkos::DynRankView<ScalarT,PHX::Device>& read_prob_distr() { 
    return prob_distr;
  }

private:
  double T_sc; // temperature scaling
  double dens_sc; // carrier density scaling
  double kB; // Boltzmann constant [eV/K]
  double eVel_pre;
  double hVel_pre;

  // characterization
  Trap<EvalT>::Type type;
   Trap<EvalT>::Location loc;
  Trap<EvalT>::EnergyDistribution enDistr;
  ScalarT enSigma;// energy distribution spread [eV]
  ScalarT esigma; // electron capture cross section [cm^2]
  ScalarT hsigma; // hole capture cross section [cm^2]
  double esigma_pwr;
  double hsigma_pwr;
  Trap<EvalT>::FieldDepXsec eFDepXsec; 
  Trap<EvalT>::FieldDepXsec hFDepXsec; 
  
  ScalarT Et; 
  ScalarT Nt;
  double g;

  // localization
  double xMin, xMax, yMin, yMax, zMin, zMax;  

  // state
  double prev_time;
  
  // DOFs and fields used to compute state
  Kokkos::DynRankView<ScalarT,PHX::Device> n;
  Kokkos::DynRankView<ScalarT,PHX::Device> p;
  Kokkos::DynRankView<ScalarT,PHX::Device> T;
  Kokkos::DynRankView<ScalarT,PHX::Device> gamma_e;
  Kokkos::DynRankView<ScalarT,PHX::Device> gamma_h;
  Kokkos::DynRankView<ScalarT,PHX::Device> Nc;
  Kokkos::DynRankView<ScalarT,PHX::Device> Nv;
  Kokkos::DynRankView<ScalarT,PHX::Device> Eg;
  Kokkos::DynRankView<ScalarT,PHX::Device> Emag;

  // time-dependent level trap state 
  Kokkos::DynRankView<ScalarT,PHX::Device> cC;
  Kokkos::DynRankView<ScalarT,PHX::Device> cV;
  Kokkos::DynRankView<ScalarT,PHX::Device> eC;
  Kokkos::DynRankView<ScalarT,PHX::Device> eV;
  Kokkos::DynRankView<ScalarT,PHX::Device> prev_prob;
  Kokkos::DynRankView<ScalarT,PHX::Device> prob;

  // time-dependent continuous distribution states
  int nL; // no of discrete levels
  std::vector<ScalarT> enLevels; // discrete energy levels
  std::vector<ScalarT> norm_densities; // normalized trap densities at mid energies
  
  Kokkos::DynRankView<ScalarT,PHX::Device> cC_distr;
  Kokkos::DynRankView<ScalarT,PHX::Device> cV_distr;
  Kokkos::DynRankView<ScalarT,PHX::Device> eC_distr;
  Kokkos::DynRankView<ScalarT,PHX::Device> eV_distr;
  Kokkos::DynRankView<ScalarT,PHX::Device> prev_prob_distr;
  Kokkos::DynRankView<ScalarT,PHX::Device> prob_distr;

  void dicretizeContDistribution();
  

}; // Trap class





// Traps (trap container)
template<typename EvalT>
class DynamicTraps
{
public:
  using ScalarT = typename EvalT::ScalarT;
  
  DynamicTraps(Teuchos::RCP<std::vector<Teuchos::RCP<Trap<EvalT>>>> traps);

  // set-up the trap system state
  void initTrapsStateWithField(const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& field) const;

  void initTrapsState(const Kokkos::DynRankView<ScalarT,PHX::Device>& edens, 
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& hdens,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& lT,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& egamma,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& hgamma,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& e_effdos,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& h_effdos,
		      const Kokkos::DynRankView<ScalarT,PHX::Device>& eff_bg) const;

  // get the total number of valid traps as defined by user
  size_t GetTrapNo() const { return trap_entries->size(); } 
  // get a trap by index
  Teuchos::RCP<Trap<EvalT>> const GetTrap(size_t itrap) { return (*trap_entries)[itrap]; } 
  // check if traps need field to cross sections
  bool WithFieldDepXsec();
  // compute traps initial state
  void computeTrapsInitialState(double t0);
  // compute traps state
  void computeTrapsState(double t);
  // save trap state
  void saveTrapsState(double t);
  
private:
  Teuchos::RCP<std::vector<Teuchos::RCP<Trap<EvalT>>>> trap_entries;
}; // DynamicTraps class





//! obtain SRH recombination rate and trapped charge
template<typename EvalT, typename Traits>
class RecombRate_DynamicTraps
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
public:
  RecombRate_DynamicTraps(const Teuchos::ParameterList& p);

  void postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

  void evaluateFields(typename Traits::EvalData d);

private:
  using ScalarT = typename EvalT::ScalarT;

  double initial_time;
  double prev_time;

  // output
  PHX::MDField<ScalarT,Cell,Point> erecomb_rate;
  PHX::MDField<ScalarT,Cell,Point> hrecomb_rate;
  PHX::MDField<ScalarT,Cell,Point> etrapped_charge;
  PHX::MDField<ScalarT,Cell,Point> htrapped_charge;
  PHX::MDField<ScalarT,Cell,Point> trapped_charge;

  // input
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> e_gamma;
  PHX::MDField<const ScalarT,Cell,Point> h_gamma;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;
  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_field;

  // traps related
  Kokkos::DynRankView<ScalarT,PHX::Device> edens;
  Kokkos::DynRankView<ScalarT,PHX::Device> hdens;
  Kokkos::DynRankView<ScalarT,PHX::Device> lT;
  Kokkos::DynRankView<ScalarT,PHX::Device> egamma;
  Kokkos::DynRankView<ScalarT,PHX::Device> hgamma;
  Kokkos::DynRankView<ScalarT,PHX::Device> e_effdos;
  Kokkos::DynRankView<ScalarT,PHX::Device> h_effdos;
  Kokkos::DynRankView<ScalarT,PHX::Device> eff_bg;
  Kokkos::DynRankView<ScalarT,PHX::Device> field;

  // traps
  Teuchos::RCP<DynamicTraps<EvalT>> traps;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]
  double T0; // temperature scaling, [K]
  double E0; // electric field scaling, [V/cm]
  double C0; // concentration scaling, [cm^(-3)]
  double R0; // recomb./gen. scaling,[#/(cm^3.s)]

  int num_points, num_dims;
  int int_rule_degree;
  int num_nodes;
  std::size_t int_rule_index;
  std::size_t basis_index;
  std::string basis_name; 
  std::string driveForce;     
  bool isSGCVFEM;

  bool withField;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  /**
   * @brief Initialize trap parameters
   */
  void initDynamicTrapsParams(const std::string& matName, 
			      const Teuchos::ParameterList& trapsPL);

}; // end of class RecombRate_DynamicTraps


}

#endif
