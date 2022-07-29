
#ifndef CHARON_AVALANCHE_CROWELLSZE_IMPL_HPP
#define CHARON_AVALANCHE_CROWELLSZE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"


/*
The Crowell-Sze avalanche model based on Baraff's curves aprroximations [1]. 
[1] Crowell, C. R., and S. M. Sze, “Temperature dependence of avalanche multiplication 
    in semiconductors,” Applied Physics Letters Vol. 9 (1966): 242–244

The ionization coefficients has the following form:

alpha_n = (1/lambda_n) * exp(C0(rn) + C1(rn)*xn + C2(rn)*xn*xn
alpha_p = (1/lambda_p) * exp(C0(rp) + C1(rp)*xp + C2(rp)*xp*xp

with 
C0(r) = -1.92 + 75.5*r - 757*r*r, 
C1(r) = 1.75e-2 - 11.9*r + 46*r*r,
C2(r) = 3.9e-4 - 1.17*r + 11.5*r*r, 

lambda_n = lambda300_e * ( tanh(E_opt_ph/2KT) / tanh(E_opt_ph/2KT300) )
lambda_p = lambda300_h * ( tanh(E_opt_ph/2KT) / tanh(E_opt_ph/2KT300) )

rn = E_opt_ph/Ei_e
rp = E_opt_ph/Ei_h

where E_opt_ph is the optical phonon energy, lambda300_e, lambda300_h are the mean 
free paths for electrons and holes at 300K, Ei_e, Ei_h are the electron and hole 
ionization energies. The E_opt_ph, lambda300_e, lambda300_h, Ei_e, Ei_h are model
parameters whose default values for Silicon are given in charon::Material_Properties.

F is the driving field, which can be GradQuasiFermi, EParallelJ, EParallelJtot,
EeffParallelJ, or EeffParallelJtot.

For numerical reasons, the field F is computed as n*Fn/(n+n0) for electrons,
and p*Fp/(p+p0) for holes. n0 and p0 are damping parameters and default to 0 and
can be set in the Avalanche Generation section through eDrForceRefDens and
hDrForceRefDens respectively. Using positive values for n0 and p0 can improve
convergence for problems where strong generation–recombination occurs in regions
with small density. 
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Avalanche_CrowellSze<EvalT, Traits>::
Avalanche_CrowellSze(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using Teuchos::ParameterList;
  using panzer::BasisIRLayout;

  RCP<ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  std::string eqset = p.get<string>("Equation Set Type");

  // get IR - It is FEM IP for SUPG-FEM and EFFPG-FEM, but is the SubCV centroid for SGCVFEM 
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  num_points = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // get Basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  num_nodes = basis->functional->dimension(1);
  basis_name = basis->name();
  RCP<DataLayout> basis_scalar = basis->functional;

  isSGCVFEM = false;  // default
  RCP<DataLayout> input_scalar = ip_scalar;
  if(eqset == "SGCVFEM Drift Diffusion") {
    RCP<DataLayout> input_scalar = basis_scalar; 
    isSGCVFEM = true; 
  } 

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  J0 = scaleParams->scale_params.J0;
  R0 = scaleParams->scale_params.R0;
  E0 = scaleParams->scale_params.E0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // initialize avalanche model parameters
  const string& materialName = p.get<string>("Material Name");
  const ParameterList& avaParamList = p.sublist("Avalanche ParameterList");
  initAvaParams(materialName, avaParamList);

  // Driving Force = GradQuasiFermi, GradPotentialParallelJ, GradPotentialParallelJtot,
  // EffectiveFieldParallelJ, or EffectiveFieldParallelJtot.
  driveForce = "EffectiveFieldParallelJ" ;   // default
  if (avaParamList.isParameter("Driving Force"))
    driveForce = avaParamList.get<string>("Driving Force");

  // Set driving force damping parameters
  eDrForceRefDens = 0.0;  // default, using positive value might improve convergence
  hDrForceRefDens = 0.0;
  if (avaParamList.isParameter("eDrForceRefDens"))
    eDrForceRefDens = avaParamList.get<double>("eDrForceRefDens");
  if (avaParamList.isParameter("hDrForceRefDens"))
    hDrForceRefDens = avaParamList.get<double>("hDrForceRefDens");

  // Define and add dependent fields according to driveForce
  if (driveForce == "GradQuasiFermi") {
    elec_drForce = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_qfp,ip_vector);
    hole_drForce = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_grad_qfp,ip_vector);
  } else if ( (driveForce == "GradPotentialParallelJ") ||
              (driveForce == "GradPotentialParallelJtot") ) {
    elec_drForce = !isSGCVFEM
      ? MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi,ip_vector)
      : MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_negpot,ip_vector);
    hole_drForce = !isSGCVFEM
      ? MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi,ip_vector)
      : MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_grad_negpot,ip_vector);
  } else if ( (driveForce == "EffectiveFieldParallelJ") ||
            (driveForce == "EffectiveFieldParallelJtot") ) {
    // if eqnset is not SGCVFEM Drift Diffusion, effective electric fields are always
    // computed because they are needed for derivatives calculations
    elec_drForce = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield,ip_vector);
    hole_drForce = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield,ip_vector);
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Error: invalid Driving Force, must be GradQuasiFermi, GradPotentialParallelJ, GradPotentialParallelJtot, EffectiveFieldParallelJ, or EffectiveFieldParallelJtot!");
  }
  this->addDependentField(elec_drForce);
  this->addDependentField(hole_drForce);

  // Always need current densities and carrier densities
  curr_dens_e = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_curr_density,ip_vector);
  curr_dens_h = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_curr_density,ip_vector);
  this->addDependentField(curr_dens_e);
  this->addDependentField(curr_dens_h);

  dens_e = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
  dens_h = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);
  this->addDependentField(dens_e);
  this->addDependentField(dens_h);

  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  this->addDependentField(latt_temp);
   
  // Add evaluated field
  avalanche_rate = MDField<ScalarT,Cell,Point>(n.field.avalanche_rate,ip_scalar);
  this->addEvaluatedField(avalanche_rate);

  std::string name = "CrowellSze_Avalanche";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Avalanche_CrowellSze<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if(isSGCVFEM) 
    basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Avalanche_CrowellSze<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Scale the damping parameters
  ScalarT scaled_eRefDens = eDrForceRefDens / C0;
  ScalarT scaled_hRefDens = hDrForceRefDens / C0;

  // Retrieve physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double eleQ = cpc.q;           // electron charge in [C]
  double kbBoltz = cpc.kb;       // Boltzmann constant in [eV/K]

  // Avalanche generation rate scaling
  ScalarT scaling = J0 / (eleQ * R0);

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    Kokkos::DynRankView<ScalarT,PHX::Device> dens_e_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> dens_h_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> latt_temp_cpts;
    if(isSGCVFEM) {
      // interpolate some fields at centroids from their values at basis
      const int num_ips = num_points;
      dens_e_cpts =
        createDynRankView(dens_e.get_static_view(),"dens_e_cpts",num_ips);
      Kokkos::deep_copy(dens_e_cpts,ScalarT(0.0));
      dens_h_cpts =
        createDynRankView(dens_h.get_static_view(),"dens_h_cpts",num_ips);
      Kokkos::deep_copy(dens_h_cpts,ScalarT(0.0));
      latt_temp_cpts =
        createDynRankView(latt_temp.get_static_view(),
                          "latt_temp_cpts",num_ips);
      Kokkos::deep_copy(latt_temp_cpts,ScalarT(0.0));

      for (int inode = 0; inode < num_nodes; ++inode)
        for (int ip = 0; ip < num_ips; ++ip) {
          dens_e_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) *
            dens_e(cell,inode);
          dens_h_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) *
            dens_h(cell,inode);
          latt_temp_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) *
            latt_temp(cell,inode);
        }
    }

    for (int point = 0; point < num_points; ++point)
    {
      ScalarT lattT, eden, hden;
      if(isSGCVFEM) {
	lattT = latt_temp_cpts(point)*T0; // [K]
	eden = dens_e_cpts(point); // scaled
        hden = dens_h_cpts(point); // scaled
      } else {
	lattT = latt_temp(cell,point)*T0; // [K]
	eden = dens_e(cell,point); // scaled
        hden = dens_h(cell,point); // scaled
      }

      // Compute 'lambda(T)' in the model [cm]
      ScalarT lambda_c = tanh(E_opt_ph/(2.*kbBoltz*lattT)) / 
	                 tanh(E_opt_ph/(2.*kbBoltz*300.0));
      ScalarT lambdae = lambda300_e * lambda_c;
      ScalarT lambdah = lambda300_h * lambda_c;
      
      ScalarT Fe = 0.0;        // driving fields
      ScalarT Fh = 0.0;
      if (driveForce == "GradQuasiFermi") {
        ScalarT egradqfp = 0.0;
        ScalarT hgradqfp = 0.0;
        for (int dim = 0; dim < num_dims; ++dim) {
          const ScalarT& gqfpe = elec_drForce(cell,point,dim);
          const ScalarT& gqfph = hole_drForce(cell,point,dim);
          egradqfp += gqfpe * gqfpe;
          hgradqfp += gqfph * gqfph;
        }
        // The density term is for damping to obtain better convergence in some cases
        Fe = std::sqrt(egradqfp) * E0 * eden/(eden + scaled_eRefDens);  // in [V/cm]
        Fh = std::sqrt(hgradqfp) * E0 * hden/(hden + scaled_hRefDens);
      } else {
        ScalarT normJe = 0.0;
        ScalarT normJh = 0.0;
        ScalarT normJtot = 0.0;
        ScalarT e_drForce_dot_Je = 0.0;
        ScalarT h_drForce_dot_Jh = 0.0;
        ScalarT e_drForce_dot_Jtot = 0.0;
        ScalarT h_drForce_dot_Jtot = 0.0;

        for (int dim = 0; dim < num_dims; ++dim) {
          const ScalarT& Je = curr_dens_e(cell,point,dim);
          const ScalarT& Jh = curr_dens_h(cell,point,dim);
          const ScalarT Jtot = Je + Jh;
          const ScalarT& elec_drF = elec_drForce(cell,point,dim);
          const ScalarT& hole_drF = hole_drForce(cell,point,dim);

          normJe += Je * Je;
          normJh += Jh * Jh;
          normJtot += Jtot * Jtot;

          e_drForce_dot_Je += elec_drF * Je;
          h_drForce_dot_Jh += hole_drF * Jh;
          e_drForce_dot_Jtot += elec_drF * Jtot;
          h_drForce_dot_Jtot += hole_drF * Jtot;

        }
        normJe = std::sqrt(normJe);
        normJh = std::sqrt(normJh);
        normJtot = std::sqrt(normJtot);

        ScalarT e_drForce_dot_J, h_drForce_dot_J;
        ScalarT e_normJ, h_normJ;
        if(driveForce == "GradPotentialParallelJ" || driveForce == "EffectiveFieldParallelJ") {
          e_drForce_dot_J = e_drForce_dot_Je;
          h_drForce_dot_J = h_drForce_dot_Jh;
          e_normJ = normJe;
          h_normJ = normJh;
        } else { // GradPotentialParallelJtot, EffectiveFieldParallelJtot
          e_drForce_dot_J = e_drForce_dot_Jtot;
          h_drForce_dot_J = h_drForce_dot_Jtot;
          e_normJ = h_normJ = normJtot;
        }

        if(e_normJ > 0.)
          Fe =  e_drForce_dot_J / e_normJ * E0 * eden/(eden + scaled_eRefDens);
        if(h_normJ > 0.)
          Fh =  h_drForce_dot_J / h_normJ * E0 * hden/(hden + scaled_hRefDens);
      }

      Fe = std::abs(Fe);
      Fh = std::abs(Fh);

      // Compute absolute values of the current densities (scaled)
      ScalarT normJe = 0.0;
      ScalarT normJh = 0.0;
      for (int dim = 0; dim < num_dims; ++dim) {
        const ScalarT& Je = curr_dens_e(cell,point,dim);
        const ScalarT& Jh = curr_dens_h(cell,point,dim);
        normJe += Je * Je;
        normJh += Jh * Jh;
      }
      normJe = std::sqrt(normJe);
      normJh = std::sqrt(normJh);

      // -----------------------------------------------------------------------
      // --- Ionization coefficients calculation depends on avalanche models ---
      // -----------------------------------------------------------------------

      // Compute impact ionization coefficients [1/cm]
      ScalarT iicoeffe=0.,iicoeffh=0.,r=0.,x=0.;
      const double C00 = -1.92; const double C01 = 75.5; const double C02 = -757.0;
      const double C10 = 1.75e-2; const double C11 = -11.9; const double C12 = 46.0;
      const double C20 = 3.9e-4; const double C21 = -1.17; const double C22 = 11.5;
      
      
      if (Fe > minField) {
	r = E_opt_ph / Ei_e;     // [1]
	x = Ei_e / lambdae / Fe; // [1]
        ScalarT C0 = C00 + C01*r + C02*r*r;
	ScalarT C1 = C10 + C11*r + C12*r*r; 
	ScalarT C2 = C20 + C21*r + C22*r*r;
	iicoeffe = exp(C0 + C1*x + C2*x*x) / lambdae;
      }
      if (Fh > minField) {
	r = E_opt_ph / Ei_h;     // [1]
	x = Ei_h / lambdah / Fh; // [1]
	ScalarT C0 = C00 + C01*r + C02*r*r;
	ScalarT C1 = C10 + C11*r + C12*r*r; 
	ScalarT C2 = C20 + C21*r + C22*r*r;
	iicoeffh = exp(C0 + C1*x + C2*x*x) / lambdah;
      }

      // Compute avalanche generation rate 
      ScalarT eAvaRate = 0.0;
      ScalarT hAvaRate = 0.0;
      if (normJe > 0.)  eAvaRate = iicoeffe*normJe;
      if (normJh > 0.)  hAvaRate = iicoeffh*normJh;

      avalanche_rate(cell,point) = (eAvaRate + hAvaRate)*scaling;

    }  // end of loop over point
  }  // end of loop over cell

}



template<typename EvalT, typename Traits>
void Avalanche_CrowellSze<EvalT, Traits>::initAvaParams(
const std::string& matName, const Teuchos::ParameterList& avaParamList) {
  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve Crowell-Sze avalanche model parameters
  E_opt_ph =  matProperty.getPropertyValue(matName, "Crowell-Sze E_opt_ph");

  lambda300_e = matProperty.getPropertyValue(matName, "Crowell-Sze Electron lambda300");
  Ei_e = matProperty.getPropertyValue(matName, "Crowell-Sze Electron Ei");

  lambda300_h = matProperty.getPropertyValue(matName, "Crowell-Sze Hole lambda300");
  Ei_h = matProperty.getPropertyValue(matName, "Crowell-Sze Hole Ei");

  // Overwrite parameters when specified by users
  if (avaParamList.isParameter("E_opt_ph"))
    E_opt_ph = avaParamList.get<double>("E_opt_ph");
  
  if (avaParamList.isParameter("lambda300_e"))
    lambda300_e = avaParamList.get<double>("lambda300_e");
  if (avaParamList.isParameter("Ei_e"))
    Ei_e = avaParamList.get<double>("Ei_e");

  if (avaParamList.isParameter("lambda300_h"))
    lambda300_h = avaParamList.get<double>("lambda300_h");
  if (avaParamList.isParameter("Ei_h"))
    Ei_h = avaParamList.get<double>("Ei_h");

  // Set minimum field below which avalanche generation = 0
  minField = 5e4;  // default, [V/cm]
  if (avaParamList.isParameter("Minimum Field"))
    minField = avaParamList.get<double>("Minimum Field");  
}



///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Avalanche_CrowellSze<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<std::string>("Equation Set Type", "?");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->sublist("Avalanche ParameterList", false, "");
  p->sublist("Avalanche ParameterList").set<std::string>("Value", "Crowell-Sze", "Use the Crowell-Sze avalanche model");
  p->sublist("Avalanche ParameterList").set<std::string>("Driving Force", "GradQuasiFermi", "Specify the driving force");
  p->sublist("Avalanche ParameterList").set<double>("eDrForceRefDens", 0., "[cm^-3]");
  p->sublist("Avalanche ParameterList").set<double>("hDrForceRefDens", 0., "[cm^-3]");
  p->sublist("Avalanche ParameterList").set<double>("Minimum Field", 0., "[V/cm]");

  p->sublist("Avalanche ParameterList").set<double>("E_opt_ph", 0., "eV");
  
  p->sublist("Avalanche ParameterList").set<double>("lambda300_e", 0., "[cm]");
  p->sublist("Avalanche ParameterList").set<double>("Ei_e", 0., "eV");

  p->sublist("Avalanche ParameterList").set<double>("lambda300_h", 0., "[cm]");
  p->sublist("Avalanche ParameterList").set<double>("Ei_h", 0., "eV");


  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
