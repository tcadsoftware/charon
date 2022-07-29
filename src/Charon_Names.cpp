
#include "Charon_Names.hpp"
#include "Teuchos_Assert.hpp"
#include "Panzer_String_Utilities.hpp"

// **********************************************************************
void charon::Names::
setDOFNames(DOF_Names& dof, const std::string& composite_prefix, const std::string& fd_suffix)
{
  //std::cout << "charon::Names has a fd_suffix value: " << fd_suffix << std::endl;
  dof.prefix = composite_prefix;
  dof.phi = composite_prefix+"ELECTRIC_POTENTIAL"+fd_suffix;
  register_pdn(dof.phi);
  dof.edensity = composite_prefix+"ELECTRON_DENSITY"+fd_suffix;
  register_pdn(dof.edensity);
  dof.hdensity = composite_prefix+"HOLE_DENSITY"+fd_suffix;
  register_pdn(dof.hdensity);
  dof.iondensity = composite_prefix+"ION_DENSITY"+fd_suffix;
  register_pdn(dof.iondensity);
  dof.elec_qpotential = composite_prefix+"ELECTRON_QUANTUM_POTENTIAL"+fd_suffix;
  register_pdn(dof.elec_qpotential);
  dof.hole_qpotential = composite_prefix+"HOLE_QUANTUM_POTENTIAL"+fd_suffix;
  register_pdn(dof.hole_qpotential);

  // this should be the same name as field.latt_temp, so that the temperature-dependent
  // evaluators can be used for both DD and DD+LatticeT formulations
  dof.latt_temp = composite_prefix+"Lattice Temperature"+fd_suffix;
}

// **********************************************************************
charon::Names::Names(int equation_dimension,
		     const std::string prefix,
		     const std::string discfields,
		     const std::string discsuffix,
                     const std::string fd_suffix /*= "" */) :
  m_equation_dimension(equation_dimension),
  m_prefix(prefix),
  m_discfields(discfields),
  m_discsuffix(discsuffix),
  m_fd_suffix(fd_suffix)
{
 {
  setDOFNames(dof, m_prefix, fd_suffix);
  setDOFNames(grad_dof, concat("GRAD_", m_prefix), fd_suffix);
  setDOFNames(dxdt, concat("DXDT_", m_prefix), fd_suffix);
  setDOFNames(res, concat("RESIDUAL_", m_prefix), fd_suffix);
  setDOFNames(scatter, concat("SCATTER_", m_prefix), fd_suffix);
  setDOFNames(src, concat("SOURCE_", m_prefix), fd_suffix);
  setDOFNames(exact, concat("EXACT_", m_prefix), fd_suffix);
  setDOFNames(error, concat("ERROR_", m_prefix), fd_suffix);
 }

 {
  op.trans = "_TRANSIENT_OP";
  op.conv = "_CONVECTION_OP";
  op.diff = "_DIFFUSION_OP";
  op.conv_diff = "_CONVECTION_DIFFUSION_OP";
  op.src = "_SOURCE_OP";
  op.div = "_DIV_OP";
  op.flux = "_FLUX_OP";
  op.laplacian = "_LAPLACIAN_OP";
  op.supg_conv = "_SUPG_CONVECTION_OP";
  op.supg_src = "_SUPG_SOURCE_OP";
  op.sym_effpg_stab = "_SymEFFPG_STABILIZATION_OP";
 }

 {
  field.space_charge = prefix+"Space Charge"+fd_suffix; 
  field.doping = prefix+"Doping"+fd_suffix;
  field.acceptor = prefix+"Acceptor Concentration"+fd_suffix;
  field.donor = prefix+"Donor Concentration"+fd_suffix;
  field.doping_raw = prefix+"Doping_Raw"+fd_suffix;
  field.acceptor_raw = prefix+"Acceptor_Raw Concentration"+fd_suffix;
  field.donor_raw = prefix+"Donor_Raw Concentration"+fd_suffix;
  field.rel_perm = prefix+"Relative Permittivity"+fd_suffix; // constant in CMF
  field.latt_temp = prefix+"Lattice Temperature"+fd_suffix;
  field.nlprho = prefix+"Nonlinear Poisson Source"+fd_suffix; // equilibrium
  field.psrc = prefix+"Poisson Source"+fd_suffix; // together with continuity equations
  field.phi_flux = prefix+"Potential Flux"+fd_suffix;
  field.intrin_conc = prefix+"Intrinsic Concentration"+fd_suffix;
  field.grad_negpot_x = prefix+"Negative Potential Gradient X"+fd_suffix;
  field.grad_negpot_y = prefix+"Negative Potential Gradient Y"+fd_suffix;

  field.elec_eff_dos = prefix+"Electron Effective DOS"+fd_suffix;
  field.elec_mobility = prefix+"Electron Mobility"+fd_suffix;
  field.elec_klaassen_mobility = prefix+"Klaassen Electron Mobility"+fd_suffix;
  field.elec_shirahata_mobility = prefix+"Shirahata Electron Mobility"+fd_suffix;
  field.elec_philips_thomas_mobility = prefix+"Philips-Thomas Electron Mobility"+fd_suffix;
  field.elec_arora_mobility = prefix+"Arora Electron Mobility"+fd_suffix; 
  field.elec_diff_coeff = prefix+"Electron Diffusion Coefficient"+fd_suffix;
  field.elec_efield = prefix+"Electron Electric Field"+fd_suffix;
  field.elec_grad_negpot = prefix+"Electron Centroid Negative Potential Gradient"+fd_suffix;
  field.elec_curr_density = prefix+"Electron Current Density"+fd_suffix;
  field.elec_curr_dens_cvedge = prefix+"Electron Current Density CV edge"+fd_suffix;
  field.elec_velocity = prefix+"Electron Velocity"+fd_suffix;
  field.elec_peclet = prefix+"Electron Peclet Number"+fd_suffix;
  field.R_e = prefix+"Electron PDE Residual"+fd_suffix;
  field.tau_stab_e = prefix+"Electron Stabilization Tau"+fd_suffix;
  field.elec_lifetime = prefix+"Electron Lifetime"+fd_suffix;
  field.elec_grad_qfp = prefix+"Electron GradQuasiFermiPotential"+fd_suffix;
  field.elec_edge_currdens = prefix+"Electron Edge Current Density" ;
  field.elec_contact_currdens = prefix+"Electron Contact Current Density"+fd_suffix;
  field.elec_deg_factor = prefix+"Electron Degeneracy Factor"+fd_suffix;
  field.elec_eff_velocity = prefix+"Electron Effective Velocity"+fd_suffix;

  field.hole_eff_dos = prefix+"Hole Effective DOS"+fd_suffix;
  field.hole_mobility = prefix+"Hole Mobility"+fd_suffix;
  field.hole_klaassen_mobility = prefix+"Klaassen Hole Mobility"+fd_suffix;
  field.hole_shirahata_mobility = prefix+"Shirahata Hole Mobility"+fd_suffix;
  field.hole_philips_thomas_mobility = prefix+"Philips-Thomas Hole Mobility"+fd_suffix;
  field.hole_arora_mobility = prefix+"Arora Hole Mobility"+fd_suffix; 
  field.hole_diff_coeff = prefix+"Hole Diffusion Coefficient"+fd_suffix;
  field.hole_efield = prefix+"Hole Electric Field"+fd_suffix;
  field.hole_grad_negpot = prefix+"Hole Centroid Negative Potential Gradient"+fd_suffix;
  field.hole_curr_density = prefix+"Hole Current Density"+fd_suffix;
  field.hole_curr_dens_cvedge = prefix+"Hole Current Density CV edge"+fd_suffix;
  field.hole_velocity = prefix+"Hole Velocity"+fd_suffix;
  field.hole_peclet = prefix+"Hole Peclet Number"+fd_suffix;
  field.R_h = prefix+"Hole PDE Residual"+fd_suffix;
  field.tau_stab_h = prefix+"Hole Stabilization Tau"+fd_suffix;
  field.hole_lifetime = prefix+"Hole Lifetime"+fd_suffix;
  field.hole_grad_qfp = prefix+"Hole GradQuasiFermiPotential"+fd_suffix;
  field.hole_edge_currdens = prefix+"Hole Edge Current Density"+fd_suffix;
  field.hole_contact_currdens = prefix+"Hole Contact Current Density"+fd_suffix;
  field.hole_deg_factor = prefix+"Hole Degeneracy Factor"+fd_suffix;
  field.hole_eff_velocity = prefix+"Hole Effective Velocity"+fd_suffix;

  field.displacement_curr_density = prefix+"Displacement Current Density"+fd_suffix;
  field.grad_phi_prev = prefix+"Potential Gradient Prev"+fd_suffix;
  field.cont_disp_curr_density = prefix+"Cont Disp Current Density"+fd_suffix;

  field.srh_recomb = prefix+"SRH Recombination"+fd_suffix;
  field.srh_deriv_e = prefix+"SRH Electron Derivative"+fd_suffix;
  field.srh_deriv_h = prefix+"SRH Hole Derivative"+fd_suffix;
  field.trap_srh_recomb = prefix+"Trap SRH Recombination"+fd_suffix;
  field.trap_srh_charge = prefix+"Trap SRH Charge"+fd_suffix;
  field.trap_srh_deriv_e = prefix+"Trap SRH Electron Derivative"+fd_suffix;
  field.trap_srh_deriv_h = prefix+"Trap SRH Hole Derivative"+fd_suffix;
  
  field.dynamic_traps_erecomb = prefix+"Dynamic Traps eRecombination"+fd_suffix;
  field.dynamic_traps_hrecomb = prefix+"Dynamic Traps hRecombination"+fd_suffix;
  field.etrapped_charge = prefix+"Electron Trapped Charge"+fd_suffix;
  field.htrapped_charge = prefix+"Hole Trapped Charge"+fd_suffix;
  field.trapped_charge = prefix+"Trapped Charge"+fd_suffix;
  field.eQF = prefix+"Electron Quasi Fermi Level"+fd_suffix;
  field.hQF = prefix+"Hole Quasi Fermi Level"+fd_suffix;

  field.rad_recomb = prefix+"Radiative Recombination"+fd_suffix;
  field.rad_deriv_e = prefix+"Radiative Electron Derivative"+fd_suffix;
  field.rad_deriv_h = prefix+"Radiative Hole Derivative"+fd_suffix;
  field.auger_recomb = prefix+"Auger Recombination"+fd_suffix;
  field.auger_deriv_e = prefix+"Auger Electron Derivative"+fd_suffix;
  field.auger_deriv_h = prefix+"Auger Hole Derivative"+fd_suffix;
  field.avalanche_rate = prefix+"Avalanche Generation"+fd_suffix;
  field.ava_deriv_e = prefix+"Avalanche Electron Derivative"+fd_suffix;
  field.ava_deriv_h = prefix+"Avalanche Hole Derivative"+fd_suffix;
  field.bbt_rate = prefix+"Band2Band Tunneling"+fd_suffix;
  field.defect_cluster_recomb = prefix+"Defect Cluster Rate"+fd_suffix;
  field.empirical_defect_recomb = prefix+"Empirical Defect Recombination"+fd_suffix;
  field.ionization_particle_strike_rate = prefix+"Ionization Particle Strike Rate"+fd_suffix;
  field.opt_gen = prefix+"Optical Generation"+fd_suffix;
  field.total_recomb = prefix+"Total Recombination"+fd_suffix;
  field.recomb_deriv_e = prefix+"Total Recombination Electron Derivative"+fd_suffix;
  field.recomb_deriv_h = prefix+"Total Recombination Hole Derivative"+fd_suffix;

  field.band_gap = prefix+"Band Gap"+fd_suffix;
  field.affinity = prefix+"Electron Affinity"+fd_suffix;
  field.eff_band_gap = prefix+"Effective Band Gap"+fd_suffix;
  field.eff_affinity = prefix+"Effective Electron Affinity"+fd_suffix;
  field.intrin_fermi = prefix+"Intrinsic Fermi Level"+fd_suffix;
  field.cond_band = prefix+"Conduction Band"+fd_suffix;
  field.vale_band = prefix+"Valence Band"+fd_suffix;
  field.ref_energy = prefix+"Reference Energy"+fd_suffix;

  field.heat_gen = prefix+"Heat Generation"+fd_suffix;
  field.heat_cap = prefix+"Heat Capacity"+fd_suffix;
  field.kappa = prefix+"Thermal Conductivity"+fd_suffix;
  field.vac_pot = prefix+"Vacuum Potential"+fd_suffix;
  field.grad_vac_pot = prefix+"Vacuum Potential Gradient"+fd_suffix;

  field.ion_mobility = prefix+"Ion Mobility"+fd_suffix;
  field.ion_diff_coeff = prefix+"Ion Diffusion Coefficient"+fd_suffix;
  field.ion_efield = prefix+"Ion Electric Field"+fd_suffix;
  field.ion_curr_density = prefix+"Ion Current Density"+fd_suffix;
  field.ion_curr_dens_x = prefix+"Ion Current Density X"+fd_suffix;
  field.ion_curr_dens_y = prefix+"Ion Current Density Y"+fd_suffix;
  field.ion_velocity = prefix+"Ion Velocity"+fd_suffix;
  field.ion_peclet = prefix+"Ion Peclet"+fd_suffix;
  field.R_ion = prefix+"Ion PDE Residual"+fd_suffix;
  field.tau_stab_ion = prefix+"Ion Stabilization Tau"+fd_suffix;
  field.ion_soret_coeff = prefix+"Ion Soret Coefficient"+fd_suffix;
  field.ion_eff_velocity = prefix+"Ion Effective Velocity"+fd_suffix;
  field.ion_thermodiff_coeff = prefix+"Ion Thermodiffusion Coefficient"+fd_suffix;

  field.mole_frac = prefix+"Mole Fraction"+fd_suffix;
  field.xMoleFrac = prefix+"xMoleFraction"+fd_suffix;
  field.yMoleFrac = prefix+"yMoleFraction"+fd_suffix;
  field.fixed_charge = prefix+"Fixed Charge"+fd_suffix; 
  field.latt_const = prefix+"Lattice Constant"+fd_suffix;
  field.e33 = prefix+"Piezoelectric Constant 33"+fd_suffix;
  field.e31 = prefix+"Piezoelectric Constant 31"+fd_suffix;
  field.c33 = prefix+"Elastic Constant 33"+fd_suffix;
  field.c13 = prefix+"Elastic Constant 13"+fd_suffix;
  field.psp = prefix+"Spontaneous Polarization"+fd_suffix;

  field.ins_genpair_density = prefix+"Ins Generated Pairs Density"+fd_suffix; 
  field.ins_htrappedcharge = prefix+"Ins Hole Trapped Charge"+fd_suffix; 

  // Fields needed for Quantum Correction
  field.e_qp_flux = prefix+"Electron Quantum Correction Potential Flux"+fd_suffix;
  field.e_qp_fieldmag = prefix+"Electron Quantum Correction Potential Field Magnitude"+fd_suffix;
  field.h_qp_flux = prefix+"Hole Quantum Correction Potential Flux"+fd_suffix;
  field.h_qp_fieldmag = prefix+"Hole Quantum Correction Potential Field Magnitude"+fd_suffix;

  // Save the initial ELECTRIC_POTENTIAL
  field.initial_phi = prefix+"Initial ELECTRIC_POTENTIAL"+fd_suffix;
  field.initial_grad_phi = prefix+"Initial GRAD_ELECTRIC_POTENTIAL"+fd_suffix; 

  register_pdn(field.doping);
  register_pdn(field.doping_raw);
  register_pdn(field.elec_eff_dos);
  register_pdn(field.hole_eff_dos);
  register_pdn(field.band_gap);
  register_pdn(field.affinity);
  register_pdn(field.eff_band_gap);
  register_pdn(field.eff_affinity);
  register_pdn(field.cond_band);
  register_pdn(field.vale_band);
 }

  // Closure Model Keys
  {
    key.material_name                  = "Material Name"+fd_suffix;
    key.radiative_recombination        = "Radiative Recombination"+fd_suffix;
    key.auger_recombination            = "Auger Recombination"+fd_suffix;
    key.trap_srh_recombination         = "Trap SRH Recombination"+fd_suffix;
    key.dynamic_traps_recombination    = "Dynamic Traps Recombination"+fd_suffix;
    key.defect_cluster_recombination   = "Defect Cluster Recombination"+fd_suffix;
    key.empirical_defect_recombination = "Empirical Defect Recombination"+fd_suffix;
    key.particle_strike                = "Particle Strike"+fd_suffix;
    key.incomplete_ionized_acceptor    = "Incomplete Ionized Acceptor"+fd_suffix;
    key.incomplete_ionized_donor       = "Incomplete Ionized Donor"+fd_suffix;
    key.tid                            = "TID"+fd_suffix;
  }
 
 // Closure Model Values
 {

 }

  // Do NOT use prefix on values
 
 { 
  layouts.basis_scalar = prefix+"BASIS SCALAR"+fd_suffix;
  layouts.basis_vector = prefix+"BASIS VECTOR"+fd_suffix;
  layouts.basis_matrix = prefix+"BASIS MATRIX"+fd_suffix;
  layouts.ip_scalar = prefix+"INTEGRATION POINT SCALAR"+fd_suffix;
  layouts.ip_vector = prefix+"INTEGRATION POINT VECTOR"+fd_suffix;
  layouts.ip_matrix = prefix+"INTEGRATION POINT MATRIX"+fd_suffix;
 }

 applySuffixes(discfields, discsuffix);

}

// **********************************************************************
charon::Names::~Names()
{
}

// **********************************************************************
void charon::Names::
applySuffixes(const std::string& fields, const std::string& suffix)
{
  // Convert the string of fields to a vector of strings
  std::vector<std::string> fields_vector;
  panzer::StringTokenizer(fields_vector, fields, ",");
  // Loop over fields
  for (std::vector<std::string>::iterator ifield=fields_vector.begin();
       ifield != fields_vector.end(); ++ifield) {
    // Loop over possibly discontinuous names to find a match
    for( std::vector<std::string*>::iterator it= pdns.begin();
         it != pdns.end(); ++it) {
      if (0 == (*it)->compare(concat(m_prefix,*ifield))
          || 0 == (*it)->compare(concat("GRAD_",m_prefix,*ifield))
          || 0 == (*it)->compare(concat("DXDT_",m_prefix,*ifield))
          || 0 == (*it)->compare(concat("RESIDUAL_",m_prefix,*ifield))
          || 0 == (*it)->compare(concat("SCATTER_",m_prefix,*ifield))
          || 0 == (*it)->compare(concat("EXACT_",m_prefix,*ifield))
          || 0 == (*it)->compare(concat("ERROR_",m_prefix,*ifield))
          ) {
        // Found a match, so append the suffix
        (*it)->append(suffix);
      }
    }
  }
}

// **********************************************************************
std::string charon::Names::concat(const std::string& s1,
                                  const std::string& s2,
                                  const std::string& s3) const
{
  std::string s = s1 + s2 + s3;
  return s;
}

// **********************************************************************
const std::string& charon::Names::prefix() const
{
  return m_prefix;
}

// **********************************************************************
const std::string& charon::Names::discfields() const
{
  return m_discfields;
}

// **********************************************************************
const std::string& charon::Names::discsuffix() const
{
  return m_discsuffix;
}

// **********************************************************************
const std::string& charon::Names::FDsuffix() const
{
  //return dof.fd_suffix;
  return m_fd_suffix;
}
