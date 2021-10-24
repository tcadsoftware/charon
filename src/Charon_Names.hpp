
#ifndef CHARON_NAMES_HPP
#define CHARON_NAMES_HPP

#include <vector>
#include <string>

namespace charon {

  /*! \brief String names for phalanx fields

    To add a new equation set:
    1. Add new DOF name to DOF_Names struct
    2. Edit the method "setDOFNames" to add in the new string
    3. Add a test to make sure the name is working correctly
  */
  class Names {

  public:

    Names(int equation_dimension, const std::string prefix,
	  const std::string discfields, const std::string discsuffix, const std::string fd_suffix = "");

    Names() = default;
    
    ~Names();

    std::string concat(const std::string& s1, const std::string& s2,
                       const std::string& s3 = "")
      const;

    void applySuffixes(const std::string& fields, const std::string& suffix);

    const std::string& prefix() const;

    const std::string& FDsuffix() const;
    
  public:

    //! Unknown (Degree of Freedom) names
    struct DOF_Names {
      //! Prefix for this struct
      std::string prefix;
      //! Suffix for this struct, supporting frequency domain analysis
      std::string fd_suffix;
      //! Electric potential (scaled, no unit)
      std::string phi;
      //! Electron density (scaled, no unit)
      std::string edensity;
      //! Hole density (scaled, no unit)
      std::string hdensity;
      //! Lattice temperature (scaled, no unit)
      std::string latt_temp;
      //! Ion density (scaled, no unit)
      std::string iondensity;
    };

    //! Suffixes to be combined with the residual names
    struct Operators {
      //! Suffix for transient term entries
      std::string trans;
      //! Suffix for convection term entries
      std::string conv;
      //! Suffix for diffusion term entries
      std::string diff;
      //! Suffix for convection-diffusion term entries
      std::string conv_diff;
      //! Suffix for general source term entries
      std::string src;
       //! Suffix for divergence term entries
      std::string div;
      //! Suffix for general flux term
      std::string flux;
      //! Suffix for Laplacian term entry
      std::string laplacian;
      //! Suffix for SUPG convection stab.
      std::string supg_conv;
      //! Suffix for SUPG source stab.
      std::string supg_src;
      //! Suffix for symEFFPG stabilization residual
      std::string sym_effpg_stab;
    };

    //! Names of intermediate fields used in Charon.
    //! All intermediate fields are saved as scaled values in FieldManager
    //! except energies and scaling parameters.

    struct Fields {

      //! Net space charge
      std::string space_charge; 
      //! Net doping
      std::string doping;
      //! Acceptor doping
      std::string acceptor;
      //! Donor doping
      std::string donor;
      //! Net raw doping (fully ionized)
      std::string doping_raw;
      //! Acceptor raw doping (fully ionized)
      std::string acceptor_raw;
      //! Donor raw doping (fully ionized)
      std::string donor_raw;
      //! Relative permittivity [1]
      std::string rel_perm;
      //! Lattice temperature
      std::string latt_temp;
      //! Nonlinear Poisson source term (equilibrium)
      std::string nlprho;
      //! Poisson source term (together with continuity equations)
      std::string psrc;
      //! Potential flux
      std::string phi_flux;
      //! Equilibrium effective intrinsic concentration
      std::string intrin_conc;
      //! Negative potential gradient
      std::string grad_negpot_x;
      std::string grad_negpot_y;

      //! Electron effective DOS
      std::string elec_eff_dos;
      //! Electron mobility
      std::string elec_mobility;
      //! Electron Klaassen mobility
      std::string elec_klaassen_mobility;
      //! Electron Shirahata mobility
      std::string elec_shirahata_mobility;
      //! Electron Philips-Thomas mobility
      std::string elec_philips_thomas_mobility;
      //! Electron diffusion coefficient
      std::string elec_diff_coeff;
      //! Electron effective electric field
      std::string elec_efield;
      //! Electron current density
      std::string elec_curr_density;
      //! Electron current density at CV edge
      std::string elec_curr_dens_cvedge;
      //! Electron velocity
      std::string elec_velocity;
      //! Electron Peclet number (FEM)
      std::string elec_peclet;
      //! Stabilized residual for electron
      std::string R_e;
      //! Stabilization parameter for electron
      std::string tau_stab_e;
      //! Electron SRH lifetime
      std::string elec_lifetime;
      //! Electron quasi fermi potential gradient
      std::string elec_grad_qfp;
      //! Electron CVFEM-SG edge current density
      std::string elec_edge_currdens;
      //! Electron contact current density
      std::string elec_contact_currdens;
      //! Electron centroid negative potential gradient
      std::string elec_grad_negpot;
      //! Electron Fermi-Dirac degeneracy factor
      std::string elec_deg_factor;
      //! Electron effective velocity, include soret contribution
      std::string elec_eff_velocity;

      //! Hole effective DOS
      std::string hole_eff_dos;
      //! Hole mobility
      std::string hole_mobility;
      //! Hole Klaassen mobility
      std::string hole_klaassen_mobility;
      //! Hole Shirahata mobility
      std::string hole_shirahata_mobility;
      //! HOle Philips-Thomas mobility
      std::string hole_philips_thomas_mobility;
      //! Hole diffusion coefficient
      std::string hole_diff_coeff;
      //! Hole effective electric field
      std::string hole_efield;
      //! Hole current density
      std::string hole_curr_density;
      //! Hole current density at CV edge
      std::string hole_curr_dens_cvedge;
      //! Hole Velocity
      std::string hole_velocity;
      //! Hole Peclet number;
      std::string hole_peclet;
      //! Stabilized residual for hole
      std::string R_h;
      //! Stabilization parameter for hole
      std::string tau_stab_h;
      //! Hole SRH lifetime
      std::string hole_lifetime;
      //! Hole quasi fermi potential gradient
      std::string hole_grad_qfp;
      //! Hole CVFEM-SG edge current density
      std::string hole_edge_currdens;
      //! Hole contact current density
      std::string hole_contact_currdens;
      //! Hole centroid negative potential gradient
      std::string hole_grad_negpot;
      //! Hole Fermi-Dirac degeneracy factor;
      std::string hole_deg_factor;
      //! Hole effective velocity, include soret contribution
      std::string hole_eff_velocity;

      //! SRH recombination rate;
      std::string srh_recomb;
      //! SRH derivative
      std::string srh_deriv_e;
      std::string srh_deriv_h;
      //! Trap SRH recombination rate
      std::string trap_srh_recomb;
      //! Trap charge and traps are modeled as SRH
      std::string trap_srh_charge;
      //! Trao SRH derivative
      std::string trap_srh_deriv_e;
      std::string trap_srh_deriv_h;
      //! Radiative recombination rate
      std::string rad_recomb;
      //! Radiative derivative
      std::string rad_deriv_e;
      std::string rad_deriv_h;
      //! Auger recombination rate
      std::string auger_recomb;
      //! Auger derivative
      std::string auger_deriv_e;
      std::string auger_deriv_h;
      //! Avalanche generation (Impact ionization) rate
      std::string avalanche_rate;
      //! Avalanche derivative
      std::string ava_deriv_e;
      std::string ava_deriv_h;
      //! Defect cluster recombination rate
      std::string defect_cluster_recomb;
      //! Empirical Defect recombination rate
      std::string empirical_defect_recomb;
      //! Ionization Particle Strike Rate
      std::string ionization_particle_strike_rate;
      //! Optical generation rate
      std::string opt_gen;
      //! Total recombination rate (sum of all recomb. rates - avalanche gen.)
      std::string total_recomb;
      //! Total recombination derivative
      std::string recomb_deriv_e;
      std::string recomb_deriv_h;

      //! Band gap in [eV] (NO BGN)
      std::string band_gap;
      //! Electron affinity in [eV] (NO BGN)
      std::string affinity;
      //! Effective band gap in [eV] (include BGN)
      std::string eff_band_gap;
      //! Effective electron affinity in [eV] (include BGN)
      std::string eff_affinity;
      //! Intrinsic Fermi energy in [eV] (include BGN)
      std::string intrin_fermi;
      //! Conduction band energy in [eV] (include BGN)
      std::string cond_band;
      //! Valence band energy in [eV] (include BGN)
      std::string vale_band;
      //! Reference energy in [eV] needed for heterogeneous devices
      std::string ref_energy;

      //! Heat generation
      std::string heat_gen;
      //! Heat capacity
      std::string heat_cap;
      //! Thermal conductivity
      std::string kappa;
      //! Vacuum potential
      std::string vac_pot;
      //! Vacuum potential gradient
      std::string grad_vac_pot;

      //! Ion mobility
      std::string ion_mobility;
      //! Ion diffusion coefficient
      std::string ion_diff_coeff;
      //! Ion electric field
      std::string ion_efield;
      //! Ion current density
      std::string ion_curr_density;
      std::string ion_curr_dens_x;
      std::string ion_curr_dens_y;
      //! Ion velocity
      std::string ion_velocity;
      //! Ion Peclet number (FEM)
      std::string ion_peclet;
      //! Stabilized residual for ion
      std::string R_ion;
      //! Stabilization parameter for ion
      std::string tau_stab_ion;
      //! Ion Soret Coefficient
      std::string ion_soret_coeff;
      //! Ion effective velocity, include soret contribution
      std::string ion_eff_velocity;
      //! Ion thermodiffusion coefficient
      std::string ion_thermodiff_coeff;

      //! X Mole fraction [unitless]
      std::string mole_frac;
      
      //! Bulk fixed charge density; 
      std::string fixed_charge; 

      //! Insulator trapped hole charge density
      std::string ins_genpair_density;
      std::string ins_htrappedcharge;

      //! Polarization lattice constant
      std::string latt_const;
      //! Piezoelectric constant 33
      std::string e33;
      //! Piezoelectric constant 31
      std::string e31;
      //! Elastic constant 33
      std::string c33;
      //! Piezoelectric constant 13
      std::string c13;
      //! Spontaneous polarization
      std::string psp;

    };

    /** Keys for groups of parameters in closure models (use
        fields names if possible to avoid redundancy)
      */
    struct Closure_Model_Keys {
      // std::string doping;
      std::string material_name;
      std::string radiative_recombination;
      std::string auger_recombination;
      std::string trap_srh_recombination;
      std::string defect_cluster_recombination;
      std::string empirical_defect_recombination;
      std::string particle_strike;
      std::string incomplete_ionized_acceptor;
      std::string incomplete_ionized_donor;
      std::string tid;
    };

    struct Closure_Model_Values {
      // not used
    };

    //! DataLayout names used by all evaluators
    struct Default_DataLayouts {
      std::string basis_scalar;
      std::string basis_vector;
      std::string basis_matrix;
      std::string ip_scalar;
      std::string ip_vector;
      std::string ip_matrix;
    };

  public:

    std::vector<std::string> axes;

    //! DOF
    DOF_Names dof;

    //! DOF Gradients
    DOF_Names grad_dof;

    //! DOF Time derivatives
    DOF_Names dxdt;

    //! Residual
    DOF_Names res;

    //! Scatter
    DOF_Names scatter;

    //! Source
    DOF_Names src;

    //! Exact solution (use for comparison with analytic/mmf solutions)
    DOF_Names exact;

    //! Error in solution (use for verification studies)
    DOF_Names error;

    //! Operator suffixes - append these onto residual names
    Operators op;

    //! Internal fields
    Fields field;

    Closure_Model_Keys key;

    Closure_Model_Values value;

    Default_DataLayouts layouts;

  private:
    
    void setDOFNames(DOF_Names& names, const std::string& composite_prefix, const std::string& fd_suffix);
    
  private:

    int m_equation_dimension;

    std::string m_prefix;
    std::string m_fd_suffix;
    std::vector <std::string*> pdns;

    void register_pdn(std::string& name) {
      pdns.push_back(&name);
    }

  };

}

#endif
