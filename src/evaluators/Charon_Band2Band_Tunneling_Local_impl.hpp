
#ifndef CHARON_BAND2BAND_TUNNELING_LOCAL_IMPL_HPP
#define CHARON_BAND2BAND_TUNNELING_LOCAL_IMPL_HPP

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
Local band2band tunneling has three models: Kane, Hurkx, Schenk

The Kane's band2band tunneling model:
Gbbt = (D_Kane * A_Kane)/Eg(T)^(1/2) * F^(gamma_Kane) * exp(-B_Kane * Eg(T)^(3/2) / F),
where F is the magnitude of the electric field.

The A_Kane, B_Kane, alpha_Kane and gamma_Kane are model parameters whose
default values are given in charon::Material_Properties.
gamma_Kane = 2 (direct transition), 5/2 (indirect transition)

D is a statistical factor (D=1 (Default value)). It can be defined as:
    D = (nie^2 - n*p)/( (n+nie)*(p+nie) ) * ( (1-abs(alpha_Kane)) - alpha_Kane)
    alpha_Kane = 0 (original Hurkx model), 1 (recombination only,), -1 (generation only)

The Hurkx band2band tunneling model is given by
Gbbt = (D_Hurkx * A_Hurkx * F^(gamma_Hurkx) * exp(-B_Hurkx/F)
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////

template<typename EvalT, typename Traits>
Band2Band_Tunneling_Local<EvalT, Traits>::
Band2Band_Tunneling_Local(
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
  int_rule_degree = ir->cubature_degree;

  // get Basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  num_nodes = basis->functional->dimension(1);
  basis_name = basis->name();
  RCP<DataLayout> basis_scalar = basis->functional;

  isSGCVFEM = false;  // default
  RCP<DataLayout> input_scalar = ip_scalar;
  if(eqset == "SGCVFEM Drift Diffusion") 
  {
     RCP<DataLayout> input_scalar = basis_scalar; 
     isSGCVFEM = true; 
  } 

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  R0 = scaleParams->scale_params.R0;
  E0 = scaleParams->scale_params.E0;
  C0 = scaleParams->scale_params.C0;

  // initialize band2band tunneling model parameters
  const string& materialName = p.get<string>("Material Name");
  const ParameterList& bbtParamList = p.sublist("Band2Band Tunneling ParameterList");
  initBBTParams(materialName, bbtParamList);

  // Driving Force = GradQuasiFermi, GradPotential, EffectiveField, 
  // GradPotentialParallelJ, GradPotentialParallelJtot,
  // EffectiveFieldParallelJ, or EffectiveFieldParallelJtot.
  driveForce = "EffectiveFieldParallelJ";   // default
  if (bbtParamList.isParameter("Driving Force"))
      driveForce = bbtParamList.get<string>("Driving Force");

  minField = 1e3;
  if (bbtParamList.isParameter("Minimum Field"))
      minField = bbtParamList.get<double>("Minimum Field");

  // Define and add dependent fields according to driveForce
  if (driveForce == "GradQuasiFermi") 
  {
      elec_drForce = MDField<const ScalarT, Cell, Point, Dim>(n.field.elec_grad_qfp, ip_vector);
      hole_drForce = MDField<const ScalarT, Cell, Point, Dim>(n.field.hole_grad_qfp, ip_vector);
  }
  else if ((driveForce == "GradPotentialParallelJ")    ||
           (driveForce == "GradPotentialParallelJtot") ||
           (driveForce == "GradPotential")) 
  {
      elec_drForce = !isSGCVFEM
          ? MDField<const ScalarT, Cell, Point, Dim>(n.grad_dof.phi, ip_vector)
          : MDField<const ScalarT, Cell, Point, Dim>(n.field.elec_grad_negpot, ip_vector);
      hole_drForce = !isSGCVFEM
          ? MDField<const ScalarT, Cell, Point, Dim>(n.grad_dof.phi, ip_vector)
          : MDField<const ScalarT, Cell, Point, Dim>(n.field.hole_grad_negpot, ip_vector);
  }
  else if ((driveForce == "EffectiveFieldParallelJ")    ||
           (driveForce == "EffectiveFieldParallelJtot") ||
           (driveForce == "EffectiveField")) 
  {
      // if eqnset is not SGCVFEM Drift Diffusion, effective electric fields are always
      // computed because they are needed for derivatives calculations
      elec_drForce = MDField<const ScalarT, Cell, Point, Dim>(n.field.elec_efield, ip_vector);
      hole_drForce = MDField<const ScalarT, Cell, Point, Dim>(n.field.hole_efield, ip_vector);
  }
  else 
  {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error: invalid Driving Force, must be GradQuasiFermi, GradPotential, EffectiveField, GradPotentialParallelJ, GradPotentialParallelJtot, EffectiveFieldParallelJ, or EffectiveFieldParallelJtot!");
  }
  this->addDependentField(elec_drForce);
  this->addDependentField(hole_drForce);

  // add the initial electric field 
  if (bAddFactor)
  {
    initial_grad_phi = MDField<const ScalarT, Cell, Point, Dim>(n.field.initial_grad_phi, ip_vector);  
    this->addDependentField(initial_grad_phi);
  }

  // Always need current densities and carrier densities
  curr_dens_e = MDField<const ScalarT, Cell, Point, Dim>(n.field.elec_curr_density, ip_vector);
  curr_dens_h = MDField<const ScalarT, Cell, Point, Dim>(n.field.hole_curr_density, ip_vector);
  this->addDependentField(curr_dens_e);
  this->addDependentField(curr_dens_h);

  // Add electron and hole densities
  dens_e = MDField<const ScalarT, Cell, Point>(n.dof.edensity, input_scalar);
  dens_h = MDField<const ScalarT, Cell, Point>(n.dof.hdensity, input_scalar);
  this->addDependentField(dens_e);
  this->addDependentField(dens_h);

  // Add effective intrinsic concentration
  intrin_conc = MDField<const ScalarT, Cell, Point>(n.field.intrin_conc, input_scalar);
  this->addDependentField(intrin_conc);

  // Add effective band gap
  eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,input_scalar);
  this->addDependentField(eff_bandgap);

  // Add additional fields for computing the threshold field
  net_doping = MDField<const ScalarT,Cell,Point>(n.field.doping,input_scalar);
  rel_perm = MDField<const ScalarT,Cell,Point>(n.field.rel_perm,input_scalar);
  this->addDependentField(net_doping);
  this->addDependentField(rel_perm);

  // Add evaluated field
  bbt_rate = MDField<ScalarT,Cell,Point>(n.field.bbt_rate,ip_scalar);
  this->addEvaluatedField(bbt_rate);

  std::string name = "Band2Band_Tunneling";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////

template<typename EvalT, typename Traits>
void
Band2Band_Tunneling_Local<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);

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
Band2Band_Tunneling_Local<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Band2Band Tunneling generation rate scaling
  ScalarT scaling = 1.0 / R0; // 1/[1/(cm^3 s)] = [cm^3 s]
  
  // use a single tolerance value for all the if statements
  double tol = 1.e-10; 

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {

    Kokkos::DynRankView<ScalarT, PHX::Device> dens_e_cpts;
    Kokkos::DynRankView<ScalarT, PHX::Device> dens_h_cpts;
    Kokkos::DynRankView<ScalarT, PHX::Device> intrin_conc_cpts;
    Kokkos::DynRankView<ScalarT, PHX::Device> eff_bandgap_cpts;

    if(isSGCVFEM) 
    {
      // interpolate some fields at centroids from their values at basis
      const int num_ips = num_points;

      dens_e_cpts = createDynRankView(dens_e.get_static_view(), "dens_e_cpts", num_ips);
      Kokkos::deep_copy(dens_e_cpts, ScalarT(0.0));
      dens_h_cpts = createDynRankView(dens_h.get_static_view(), "dens_h_cpts", num_ips);
      Kokkos::deep_copy(dens_h_cpts, ScalarT(0.0));
      intrin_conc_cpts = createDynRankView(intrin_conc.get_static_view(), "intrin_conc_cpts", num_ips);
      Kokkos::deep_copy(intrin_conc_cpts, ScalarT(0.0));
      eff_bandgap_cpts = createDynRankView(eff_bandgap.get_static_view(),"eff_bandgap_cpts",num_ips);
      Kokkos::deep_copy(eff_bandgap_cpts,ScalarT(0.0));

      for (int inode = 0; inode < num_nodes; ++inode)
      {
         for (int ip = 0; ip < num_ips; ++ip) 
         {
            dens_e_cpts(ip)      += (workset.bases[basis_index])->basis_scalar(cell, inode, ip) * dens_e(cell, inode);
            dens_h_cpts(ip)      += (workset.bases[basis_index])->basis_scalar(cell, inode, ip) * dens_h(cell, inode);
            intrin_conc_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell, inode, ip) * intrin_conc(cell, inode);
            eff_bandgap_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell, inode, ip) * eff_bandgap(cell,inode);
         }
      }
    }

    for (int point = 0; point < num_points; ++point)
    {
     double x = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
     double y = 0.0, z = 0.0;
     if (num_dims == 2)
       y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
     if (num_dims == 3)
     {
       y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
       z = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,2);
     }
  
     if ((x >= Xmin) && (x <= Xmax) && (y >= Ymin) && (y <= Ymax) && (z >= Zmin) && (z <= Zmax))
     {  // compute the band2band tunneling rate
      ScalarT effEg, eden, hden, nie;
      if (isSGCVFEM) 
      {
          eden  = dens_e_cpts(point);       // scaled
          hden  = dens_h_cpts(point);       // scaled
          nie   = intrin_conc_cpts(point);  // scaled
          effEg = eff_bandgap_cpts(point);  // [eV]
      }
      else 
      {
          eden  = dens_e(cell, point);       // scaled
          hden  = dens_h(cell, point);       // scaled
          nie   = intrin_conc(cell, point);  // scaled
          effEg = eff_bandgap(cell, point);  // eV
      }

      ScalarT Fe = 0.0;        // driving fields
      ScalarT Fh = 0.0;
      ScalarT F = 0.0;
      if (driveForce == "GradQuasiFermi" || driveForce == "GradPotential" || driveForce == "EffectiveField") 
      {
          ScalarT egradqfp = 0.0;
          ScalarT hgradqfp = 0.0;
          for (int dim = 0; dim < num_dims; ++dim) 
          {
              const ScalarT& gqfpe = elec_drForce(cell, point, dim);
              const ScalarT& gqfph = hole_drForce(cell, point, dim);
              egradqfp += gqfpe * gqfpe;
              hgradqfp += gqfph * gqfph;
          }

          if (Sacado::ScalarValue<ScalarT>::eval(egradqfp) > tol)
              Fe = std::sqrt(egradqfp);  // scaled
          if (Sacado::ScalarValue<ScalarT>::eval(hgradqfp) > tol)
              Fh = std::sqrt(hgradqfp);  // scaled
      }
      else 
      {
          ScalarT normJe = 0.0;
          ScalarT normJh = 0.0;
          ScalarT normJtot = 0.0;
          ScalarT e_drForce_dot_Je = 0.0;
          ScalarT h_drForce_dot_Jh = 0.0;
          ScalarT e_drForce_dot_Jtot = 0.0;
          ScalarT h_drForce_dot_Jtot = 0.0;

          for (int dim = 0; dim < num_dims; ++dim) 
          {
              const ScalarT& Je = curr_dens_e(cell, point, dim);
              const ScalarT& Jh = curr_dens_h(cell, point, dim);
              const ScalarT Jtot = Je + Jh;
              const ScalarT& elec_drF = elec_drForce(cell, point, dim);
              const ScalarT& hole_drF = hole_drForce(cell, point, dim);

              normJe += Je * Je;
              normJh += Jh * Jh;
              normJtot += Jtot * Jtot;

              e_drForce_dot_Je += elec_drF * Je;
              h_drForce_dot_Jh += hole_drF * Jh;
              e_drForce_dot_Jtot += elec_drF * Jtot;
              h_drForce_dot_Jtot += hole_drF * Jtot;
          }

          if (Sacado::ScalarValue<ScalarT>::eval(normJe) > tol)
              normJe = std::sqrt(normJe);
          else
              normJe = 0.0; 
          if (Sacado::ScalarValue<ScalarT>::eval(normJh) > tol)
              normJh = std::sqrt(normJh);
          else
              normJh = 0.0; 
          if (Sacado::ScalarValue<ScalarT>::eval(normJtot) > tol)
              normJtot = std::sqrt(normJtot);
          else
              normJtot = 0.0; 

          ScalarT e_drForce_dot_J, h_drForce_dot_J;
          ScalarT e_normJ, h_normJ;
          if (driveForce == "GradPotentialParallelJ" || driveForce == "EffectiveFieldParallelJ") 
          {
              e_drForce_dot_J = e_drForce_dot_Je;
              h_drForce_dot_J = h_drForce_dot_Jh;
              e_normJ = normJe;
              h_normJ = normJh;
          }
          else  // GradPotentialParallelJtot, EffectiveFieldParallelJtot
          {
              e_drForce_dot_J = e_drForce_dot_Jtot;
              h_drForce_dot_J = h_drForce_dot_Jtot;
              e_normJ = normJtot;
              h_normJ = normJtot;
          }

          if (Sacado::ScalarValue<ScalarT>::eval(e_normJ) > tol) 
              Fe = e_drForce_dot_J / e_normJ; // scaled
          if (Sacado::ScalarValue<ScalarT>::eval(h_normJ) > tol) 
              Fh = h_drForce_dot_J / h_normJ; // scaled
      }

      Fe = std::abs(Fe); // scaled
      Fh = std::abs(Fh); // scaled
      F = (Fe + Fh) * E0 * 0.5; // in [V/cm]
      
      if (bbtModel == "Kane") 
      {
        // Kane model
        ScalarT D = 1.0; // Default value
        ScalarT ni2 = nie * nie; //scaled
        ScalarT numer = ni2 - eden * hden;
        ScalarT denom = (eden + nie) * (hden + nie);
        D = (numer/denom)*( 1.0-std::abs(alpha_Kane) ) - alpha_Kane ;
        // alpha_Kane=0: original Hurkx model; alpha_Kane=1: only recombination; alpha_Kane=-1: only generation

        ScalarT value = 0.0;
        // when F is very small it can lead to numerical issues
        if (Sacado::ScalarValue<ScalarT>::eval(F) > minField)
          value = std::exp(-B_Kane*std::pow(effEg,1.5)/F); 

        // add additional field factor to reduce unphysial tunneling current at zero drain voltage 
        ScalarT addFactor = 1.0;
        if (bAddFactor)  // do the following only when user wants to
        { 
          // compute the initial electric field magnitude
          ScalarT Finit = 0.0; 
          for (int dim = 0; dim < num_dims; ++dim)
            Finit += initial_grad_phi(cell,point,dim) * initial_grad_phi(cell,point,dim);  

          if (Sacado::ScalarValue<ScalarT>::eval(Finit) > tol)
            Finit = std::sqrt(Finit) * E0;
          
          double Fth = Sacado::ScalarValue<ScalarT>::eval(Finit);  
          // if ((Sacado::ScalarValue<ScalarT>::eval(std::abs(F-Fth)) > tol) &&
          if (Sacado::ScalarValue<ScalarT>::eval(Fth) > tol)
            addFactor = std::pow(std::abs(F-Fth)/Fth, beta_Kane);
        }

        // Kane band-to-band tunneling rate (scaled)
        bbt_rate(cell,point) = D*A_Kane*std::pow(F,gamma_Kane)/std::sqrt(effEg)* value * scaling * addFactor;
      }

      else if (bbtModel == "Hurkx")  // Hurkx model
      {
        // Hurkx model
        ScalarT D = 1.0; // Default value
        ScalarT ni2 = nie * nie; //scaled
        ScalarT numer = ni2 - eden * hden;
        ScalarT denom = (eden + nie) * (hden + nie);
        D = (numer / denom) * (1.0 - std::abs(alpha_Hurkx)) - alpha_Hurkx;
        // alpha_Hurkx=0: original Hurkx model; alpha_Hurkx=1: only recombination; alpha_Hurkx=-1: only generation
      
        ScalarT value = 0.0;
        // when F is very small it can lead to numerical issues
        if (Sacado::ScalarValue<ScalarT>::eval(F) > minField)
          value = std::exp(-B_Hurkx / F); 

        // add an additional field factor to reduce unphysical tunneling current at zero drain voltage
        ScalarT addFactor = 1.0;
        if (bAddFactor)  // do the following only when user wants to
        { 
          ScalarT Finit = 0.0;
          for (int dim = 0; dim < num_dims; ++dim)
            Finit += initial_grad_phi(cell,point,dim) * initial_grad_phi(cell,point,dim);

          if (Sacado::ScalarValue<ScalarT>::eval(Finit) > tol)
            Finit = std::sqrt(Finit) * E0;

          double Fth = Sacado::ScalarValue<ScalarT>::eval(Finit);
          // if ( (Sacado::ScalarValue<ScalarT>::eval(std::abs(F-Fth)) > tol) &&
          if (Sacado::ScalarValue<ScalarT>::eval(Fth) > tol)
            addFactor = std::pow(std::abs(F-Fth)/Fth, beta_Hurkx); 
        }

        // Hurkx band-to-band tunneling rate (scaled)
        bbt_rate(cell, point) = D * A_Hurkx * std::pow(F, gamma_Hurkx) * value * scaling * addFactor;
      }

      else if (bbtModel == "Schenk") 
      {
        //Schenk model
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Error: not implemented Schenk Evaluation Field!");
      }
     }  // end of the calculation if within a spatial window
   
     else  // set the tunneling rate to 0 outside a spatial window
      bbt_rate(cell, point) = 0.0; 

    }  // end of loop over point
  
  }  // end of loop over cell

}

///////////////////////////////////////////////////////////////////////////////
//
//  initBBTParams()
//
///////////////////////////////////////////////////////////////////////////////

template<typename EvalT, typename Traits>
void Band2Band_Tunneling_Local<EvalT, Traits>::initBBTParams(
  const std::string& matName, const Teuchos::ParameterList& bbtParamList) 
{
  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve band2band tunneling model
  bbtModel = "Kane"; //use Kane model by default
  if (bbtParamList.isParameter("Model"))
      bbtModel = bbtParamList.get<std::string>("Model");

  if (bbtModel == "Kane")
  {
      // Retrieve band2band tunneling model parameters
      A_Kane = matProperty.getPropertyValue(matName, "Kane A");
      B_Kane = matProperty.getPropertyValue(matName, "Kane B");
      alpha_Kane = matProperty.getPropertyValue(matName, "Kane alpha");
      gamma_Kane = matProperty.getPropertyValue(matName, "Kane gamma");
      beta_Kane = matProperty.getPropertyValue(matName, "Kane beta");

      // Overwrite parameters when specified by users
      if (bbtParamList.isParameter("Kane A"))
          A_Kane = bbtParamList.get<double>("Kane A");
      if (bbtParamList.isParameter("Kane B"))
          B_Kane = bbtParamList.get<double>("Kane B");
      if (bbtParamList.isParameter("Kane alpha"))
          alpha_Kane = bbtParamList.get<double>("Kane alpha");
      if (bbtParamList.isParameter("Kane gamma"))
          gamma_Kane = bbtParamList.get<double>("Kane gamma");
      if (bbtParamList.isParameter("Kane beta"))
      {
          beta_Kane = bbtParamList.get<double>("Kane beta");
          bAddFactor = true; 
      }

  } //end Kane init BBT parameters
  
  else if (bbtModel == "Hurkx")
  {
      // Retrieve band2band tunneling model parameters
      A_Hurkx = matProperty.getPropertyValue(matName, "Hurkx A");
      B_Hurkx = matProperty.getPropertyValue(matName, "Hurkx B");
      alpha_Hurkx = matProperty.getPropertyValue(matName, "Hurkx alpha");
      gamma_Hurkx = matProperty.getPropertyValue(matName, "Hurkx gamma");
      beta_Hurkx = matProperty.getPropertyValue(matName, "Hurkx beta");

      // Overwrite parameters when specified by users
      if (bbtParamList.isParameter("Hurkx A"))
          A_Hurkx = bbtParamList.get<double>("Hurkx A");
      if (bbtParamList.isParameter("Hurkx B"))
          B_Hurkx = bbtParamList.get<double>("Hurkx B");
      if (bbtParamList.isParameter("Hurkx alpha"))
          alpha_Hurkx = bbtParamList.get<double>("Hurkx alpha");
      if (bbtParamList.isParameter("Hurkx gamma"))
          gamma_Hurkx = bbtParamList.get<double>("Hurkx gamma");
      if (bbtParamList.isParameter("Hurkx beta"))
      {
          beta_Hurkx = bbtParamList.get<double>("Hurkx beta");
          bAddFactor = true; 
      }

  } //end Hurkx init BBT parameters
  
  else if (bbtModel == "Schenk")
  {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Error: invalid Schenk init BBT parameters, not implemented yet!");

  } //end Schenk init BBT parameters

  // set the min and max coordinate values
  Xmin = -1e100; Xmax = 1e100;
  Ymin = -1e100; Ymax = 1e100;
  Zmin = -1e100; Zmax = 1e100;

  if (bbtParamList.isParameter("Xmin"))  Xmin = bbtParamList.get<double>("Xmin");
  if (bbtParamList.isParameter("Xmax"))  Xmax = bbtParamList.get<double>("Xmax");
  if (bbtParamList.isParameter("Ymin"))  Ymin = bbtParamList.get<double>("Ymin");
  if (bbtParamList.isParameter("Ymax"))  Ymax = bbtParamList.get<double>("Ymax");
  if (bbtParamList.isParameter("Zmin"))  Zmin = bbtParamList.get<double>("Zmin");
  if (bbtParamList.isParameter("Zmax"))  Zmax = bbtParamList.get<double>("Zmax");

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////

template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Band2Band_Tunneling_Local<EvalT, Traits>::getValidParameters() const
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

  p->sublist("Band2Band Tunneling ParameterList", false, "");
  p->sublist("Band2Band Tunneling ParameterList").set<std::string>("Value", "Local", "Use local Band2Band Tunneling model"); // Value = Local
  p->sublist("Band2Band Tunneling ParameterList").set<std::string>("Model", "Kane", "Use the Kane local Band2Band Tunneling model"); // Model = Kane, Hurkx, Schenk
  p->sublist("Band2Band Tunneling ParameterList").set<std::string>("Driving Force", "GradQuasiFermi", "Specify the driving force"); // Driving Force = GradQuasiFermi, GradPotential, EffectiveField, 
                                                                                                                                       //                 GradPotentialParallelJ, GradPotentialParallelJtot,
                                                                                                                                       //                 EffectiveFieldParallelJ, or EffectiveFieldParallelJtot
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Minimum Field", 1.0, "Minimum value of the electric field from which the BTB model is turn on");
  
  // Parameter List for Kane model
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Kane A",      0.,        "[eV^(1/2)/(cm.s.V^2)]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Kane B",      0.,        "[V/(cm.eV^(3/2))]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Kane gamma",  0.,        "[1]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Kane alpha",  0.,        "[1]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Kane beta",  0.,        "[1]");

  // Parameter List for Hurkx model
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Hurkx A",      0.,        "[eV^(1/2)/(cm.s.V^2)]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Hurkx B",      0.,        "[V/(cm.eV^(3/2))]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Hurkx gamma",  0.,        "[1]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Hurkx alpha",  0.,        "[1]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Hurkx beta",  0.,        "[1]");

  // Parameter List for Schenk model
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Schenk A",  0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Schenk B",  0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Schenk HW", 0., "[]");

  // Spatial coordinates 
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Xmin", 0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Xmax", 0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Ymin", 0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Ymax", 0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Zmin", 0., "[]");
  p->sublist("Band2Band Tunneling ParameterList").set<double>("Zmax", 0., "[]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
