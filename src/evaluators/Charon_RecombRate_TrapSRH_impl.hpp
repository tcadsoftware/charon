
#ifndef CHARON_RECOMBRATE_TRAPSRH_IMPL_HPP
#define CHARON_RECOMBRATE_TRAPSRH_IMPL_HPP

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/math/special_functions/airy.hpp>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Material_Properties.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
// #include <chrono>
// #include <ctime> 

const int MAX_NUM_TRAPS = 50;
const int MAX_NUM_ENERGIES = 200;
const int ENERGY_GRID_RATIO = 4;
const double MIN_ELECTRIC_FIELD = 1e5;  // [V/m]
const double MAX_BETA = 30;
const double MAX_DISTANCE = 20E-9;  // [m]
const double MAX_ENHANCE_FACTOR = 1E8;

const std::vector<double> abscissas20_dos = {
 0.005718877300142 ,0.0227446680053324,0.0506878420247731,
 0.0889090923694288,0.1365339617565427,0.1924728491834011,
 0.255445938746293 ,0.3240124803676208,0.3966037525316293,
 0.4715589528970582,0.5471631956847183,0.6216867465723728,
 0.6934245975978842,0.7607354769985799,0.8220794024218905,
 0.8760529211255577,0.9214212407914146,0.9571465601852601,
 0.9824122393491017,0.9966459942144552
};
const std::vector<double> weights20_dos = {
 0.0008633073867834,0.0034039498722146,0.0074766474904512,
 0.0128476655801293,0.0192069589304947,0.0261842900158401,
 0.0333684924530623,0.0403288889874457,0.046637762368297 ,
 0.0518927231195774,0.0557378232382089,0.0578823290177204,
 0.0581161862319689,0.0563213807476829,0.0524786088535414,
 0.0466689146143901,0.0390702202307166,0.0299489979627076,
 0.0196481244884538,0.0085833950769805
};
const int order20_dos = static_cast<int>(abscissas20_dos.size());


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_TrapSRH<EvalT, Traits>::
RecombRate_TrapSRH(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using PHX::DataLayout;
  using PHX::MDField;

  // validate the parameter list
  RCP<ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // retrieve parameters
  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));
  const string& matName = p.get<string>("Material Name");
  const string& eqnSetType = p.get<string>("Equation Set Type");
  driveForce = p.get<string>("Driving Force");

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  T0 = scaleParams->scale_params.T0;
  E0 = scaleParams->scale_params.E0;
  X0 = scaleParams->scale_params.X0;
  C0 = scaleParams->scale_params.C0;

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

  // Get and validate the Trap SRH ParameterList
  //RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  //trapSRHParamList.validateParameters(*valid_params);

  // initialize the trap parameters
  const ParameterList& trapSRHParamList = p.sublist("Trap SRH ParameterList");
  initTrapSRHParams(matName, trapSRHParamList);

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

  // evaluated fields
  trap_srh_rate = MDField<ScalarT,Cell,Point>(n.field.trap_srh_recomb,output_scalar);
  trap_srh_charge = MDField<ScalarT,Cell,Point>(n.field.trap_srh_charge,output_scalar);
  trap_srh_deriv_e = MDField<ScalarT,Cell,Point>(n.field.trap_srh_deriv_e,output_scalar);
  trap_srh_deriv_h = MDField<ScalarT,Cell,Point>(n.field.trap_srh_deriv_h,output_scalar);
  this->addEvaluatedField(trap_srh_rate);
  this->addEvaluatedField(trap_srh_charge);
  this->addEvaluatedField(trap_srh_deriv_e);
  this->addEvaluatedField(trap_srh_deriv_h);

  // dependent fields
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos, input_scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos, input_scalar);
  eff_bandgap =  MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, input_scalar);
  latt_temp =  MDField<const ScalarT,Cell,Point>(n.field.latt_temp, input_scalar);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(eff_bandgap);
  this->addDependentField(latt_temp);

  // determine the dependent driving force
  if (driveForce == "GradQuasiFermi")
  {
    elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_qfp, input_vector);
    hole_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_grad_qfp, input_vector);
    fieldSign = 1.0;
  }
  else if (driveForce == "EffectiveField")
  {
    elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield, input_vector);
    hole_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield, input_vector);
    fieldSign = 1.0;
  }
  else if (driveForce == "GradPotential")
  {
    if (isSGCVFEM)
    {
      elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_negpot, input_vector);
      hole_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_grad_negpot, input_vector);
      fieldSign = 1.0;
    }
    else  // for FEM-based discretization
    {
      elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi, input_vector);
      hole_field = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.phi, input_vector);
      fieldSign = -1.0;
    }
  }
  this->addDependentField(elec_field);
  this->addDependentField(hole_field);

  // universal constants
  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  hbar = phyConst.hbar;  // reduced Planck constant [J.s]
  q = phyConst.q;        // elemental charge [C]
  m0 = phyConst.m0;      // free electron mass [kg]
  pi = phyConst.pi;
  kb = phyConst.kb;      // Boltzmann constant in [eV/K]

  std::string name = "RecombRate_TrapSRH";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_TrapSRH<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_TrapSRH<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using panzer::index_t;
  using std::string;

  int numofTraps = trapDensity.size();

  for (int cell = 0; cell < workset.num_cells; ++cell)
  {
    Kokkos::DynRankView<ScalarT,PHX::Device> edens_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> hdens_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> elec_effdos_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> hole_effdos_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> eff_bandgap_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> latt_temp_cpts;

    // interpolate scalar fields at nodes to the centroids of subcv for SGCVFEM
    if (isSGCVFEM)
    {
      // initialization
      const int num_ips = num_points;
      edens_cpts = createDynRankView(edensity.get_static_view(),"edens_cpts",num_ips);
      Kokkos::deep_copy(edens_cpts,ScalarT(0.0));
      hdens_cpts = createDynRankView(hdensity.get_static_view(),"hdens_cpts",num_ips);
      Kokkos::deep_copy(hdens_cpts,ScalarT(0.0));
      elec_effdos_cpts = createDynRankView(elec_effdos.get_static_view(),"elec_effdos_cpts",num_ips);
      Kokkos::deep_copy(elec_effdos_cpts,ScalarT(0.0));
      hole_effdos_cpts = createDynRankView(hole_effdos.get_static_view(),"hole_effdos_cpts",num_ips);
      Kokkos::deep_copy(hole_effdos_cpts,ScalarT(0.0));
      eff_bandgap_cpts = createDynRankView(eff_bandgap.get_static_view(),"eff_bandgap_cpts",num_ips);
      Kokkos::deep_copy(eff_bandgap_cpts,ScalarT(0.0));
      latt_temp_cpts = createDynRankView(latt_temp.get_static_view(),"latt_temp_cpts",num_ips);
      Kokkos::deep_copy(latt_temp_cpts,ScalarT(0.0));

      // interpolation
      for (int inode = 0; inode < num_nodes; ++inode)
      {
        for (int ip = 0; ip < num_ips; ++ip)
        {
          edens_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * edensity(cell,inode);
          hdens_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * hdensity(cell,inode);
          elec_effdos_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * elec_effdos(cell,inode);
          hole_effdos_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * hole_effdos(cell,inode);
          eff_bandgap_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * eff_bandgap(cell,inode);
          latt_temp_cpts(ip) += (workset.bases[basis_index])->basis_scalar(cell,inode,ip) * latt_temp(cell,inode);
        }
      }
    }  // end of the node-to-centroid interpolation for SGCVFEM

    // compute trap_srh_rate at IPs (FEM IPs for SUPG-FEM and EFFPG-FEM, but SubCV Centroids for SGCVFEM)
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT n, p, Nc, Nv, Eg, kbT, latt;
      if (isSGCVFEM)  // use values at the centroids of SubCV for SGCVFEM
      {
        n = edens_cpts(point);  // scaled
        p = hdens_cpts(point);
        Nc = elec_effdos_cpts(point);
        Nv = hole_effdos_cpts(point);
        Eg = eff_bandgap_cpts(point);  // [eV]
        kbT = latt_temp_cpts(point) * T0 * kb;  // [eV]
        latt = latt_temp_cpts(point) * T0;      // [K]
      }
      else  // use values at the IPs of finite elements for SUPG-FEM and EFFPG-FEM
      {
        n = edensity(cell,point);  // scaled
        p = hdensity(cell,point);
        Nc = elec_effdos(cell,point);
        Nv = hole_effdos(cell,point);
        Eg = eff_bandgap(cell,point);  // [eV]
        kbT = latt_temp(cell,point) * T0 * kb; // [eV]
        latt = latt_temp(cell,point) * T0;     // [K]
      }
/*
      if ( (Sacado::ScalarValue<ScalarT>::eval(n) <= 0.0) ||  // or, not and
           (Sacado::ScalarValue<ScalarT>::eval(p) <= 0.0) )
      {
        trap_srh_rate(cell,point) = 0.0;
        trap_srh_charge(cell,point) = 0.0;
        trap_srh_deriv_e(cell,point) = 0.0;
        trap_srh_deriv_h(cell,point) = 0.0;
        std::cout << "n=" << n << ", p=" << p << ", trap_srh_rate=" << trap_srh_rate(cell,point) << std::endl;
        break;
      }
*/
      if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
          (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
      {
        // calculate trap_srh_rate when n and p are both positive
        trap_srh_rate(cell,point) = 0.0;  // zero out trap_srh_recomb
        trap_srh_charge(cell,point) = 0.0;  // zero out 
        trap_srh_deriv_e(cell,point) = 0.0;  // zero out 
        trap_srh_deriv_h(cell,point) = 0.0;  // zero out 

        // compute the electric field amplitude
        ScalarT Fn = 0.0, Fp = 0.0;
        for (int dim = 0; dim < num_dims; ++dim)  // fieldSign is not necessary here since |Fn| and |Fp| are needed here
        {
           Fn += elec_field(cell,point,dim) * elec_field(cell,point,dim) * fieldSign * fieldSign; 
           Fp += hole_field(cell,point,dim) * hole_field(cell,point,dim) * fieldSign * fieldSign;
        }
        Fn = std::sqrt(Fn) * E0 * 100.0;  // [V/m], 100 is used to convert from [V/cm] to [V/m]
        Fp = std::sqrt(Fp) * E0 * 100.0;

        // obtain the spatial coordinates at IPs (FEM IPs for SUPG-FEM and EFFPG-FEM, but subcv centroid for SGCVFEM)
        double xcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
        double ycoord = 0.0, zcoord = 0.0;
        if (num_dims == 2)
          ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
        if (num_dims == 3)
        {
          ycoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
          zcoord = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,2);
        }
      
        // loop over the number of traps
        for (int itrap = 0; itrap < numofTraps; ++itrap)
        {
          ScalarT taun0 = eLifeTime[itrap]; // [s]
          ScalarT taup0 = hLifeTime[itrap];

          if (!eTimeGiven[itrap])  // compute taun0 from electron cross section
            taun0 *= std::sqrt(300.0 / latt);
          if (!hTimeGiven[itrap])  // compute taup0 from hole cross section
            taup0 *= std::sqrt(300.0 / latt);

          // obtain the field-independent lifetimes (scaled)
          ScalarT taun = taun0 / t0;
          ScalarT taup = taup0 / t0;

          // obtain the tunneling model name
          string eTunModel = eTunnelModel[itrap];
          string hTunModel = hTunnelModel[itrap];

          // initialize the field enhancement factor
          double gn = 1.0, gp = 1.0;

          // compute the field enhancement factor for the electron lifetime
          if ( (eTunModel.find("Schenk") != std::string::npos) )  // find "Schenk"
          {
            if ((eTunModel != "Schenk NewDOS") && (Sacado::ScalarValue<ScalarT>::eval(Fn) > MIN_ELECTRIC_FIELD) )
              gn = evalSchenkFieldFactor(Fn, kbT, Eg, itrap, "Electron");

            else if (eTunModel == "Schenk NewDOS")
            {
              string tunDir = eTunnelDir[itrap];
              double spcLoc = xcoord;  // default to "X"
              if (tunDir == "Y") spcLoc = ycoord;
              if (tunDir == "Z") spcLoc = zcoord;
              gn = evalFieldFactorWithNewDOS(Fn, kbT, Eg, itrap, "Electron", spcLoc);
            }
          }

          // compute the field enhancement factor for the hole lifetime
          if ( (hTunModel.find("Schenk") != std::string::npos) )  // find "Schenk"
          {
            if ((hTunModel != "Schenk NewDOS") && (Sacado::ScalarValue<ScalarT>::eval(Fp) > MIN_ELECTRIC_FIELD) )
              gp = evalSchenkFieldFactor(Fp, kbT, Eg, itrap, "Hole");

            else if (hTunModel == "Schenk NewDOS")
            {
              string tunDir = hTunnelDir[itrap];
              double spcLoc = xcoord;  // default to "X"
              if (tunDir == "Y") spcLoc = ycoord;
              if (tunDir == "Z") spcLoc = zcoord;
              gp = evalFieldFactorWithNewDOS(Fp, kbT, Eg, itrap, "Hole", spcLoc);
            }
          }

          // if (gn > MAX_ENHANCE_FACTOR) gn = MAX_ENHANCE_FACTOR;
          // if (gp > MAX_ENHANCE_FACTOR) gp = MAX_ENHANCE_FACTOR;
          // std::cout << "Trap " << itrap << ", Fn = " << Fn << ", Fp = " << Fp << ", gn = " << gn << ", gp = " << gp << std::endl;

          // add the field enhancement factor to the lifetime
          taun = taun / gn;
          taup = taup / gp;

          double Et = energyLevel[itrap];
          ScalarT ni2 = Nc * Nv * std::exp(-Eg / kbT);  // scaled
          ScalarT n1 = Nc * std::exp(-Et / kbT);
          ScalarT p1 = Nv * std::exp((Et - Eg) / kbT);

          // scaled and unscaled Rsrh have the same expression
          ScalarT numer = n*p - ni2;
          ScalarT denom = taup*(n+n1) + taun*(p+p1);
          ScalarT Rsrh = numer / denom;       // scaled
          ScalarT deriv_e = p / denom - numer*taup / (denom*denom);
          ScalarT deriv_h = n / denom - numer*taun / (denom*denom);

          // sum up the contribution from all specified traps
          trap_srh_rate(cell,point) += Rsrh;
          trap_srh_deriv_e(cell,point) += deriv_e;
          trap_srh_deriv_h(cell,point) += deriv_h;

          // compute and sum up the trap charge (scaled)
          ScalarT ft = 0.0;  // trap occupation within [0,1]
          ScalarT Qt = 0.0;  // trap charge (scaled)
          if (trapType[itrap] == "Acceptor")      // electron capture (close to Ec)
          { 
            ft = (taup*n + taun*p1) / denom; 
            Qt = - trapDensity[itrap] * ft / C0;  // carry -e when occupied
          }
          else if (trapType[itrap] == "Donor")    // hole capture (close to Ev)
          {
            ft = (taup*n1 + taun*p) / denom; 
            Qt = trapDensity[itrap] * ft / C0;    // carry +e when occupied
          }

          trap_srh_charge(cell,point) += Qt;  // will be added to the RHS of the Poisson equation
       }  // end of loop over traps

     }  // end of if (n > 0. && p > 0.)
     else
     {
       trap_srh_rate(cell,point) = 0.0;
       trap_srh_charge(cell,point) = 0.0;
       trap_srh_deriv_e(cell,point) = 0.0;
       trap_srh_deriv_h(cell,point) = 0.0;
     } 
      
    }  // end of loop over points
  }  // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSchenkFieldFactor()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::evalSchenkFieldFactor(const ScalarT& fieldAmp,
  const ScalarT& kbT, const ScalarT& bandGap, const int& itrap, const std::string& carrType)
{
  double effMass = 0.0;
  double Et = 0.0;
  std::string tunModel;

  if (carrType == "Electron")
  {
    effMass = eEffMass[itrap];
    Et = energyLevel[itrap];
    tunModel = eTunnelModel[itrap];
  }
  else if (carrType == "Hole")
  {
    effMass = hEffMass[itrap];
    Et = Sacado::ScalarValue<ScalarT>::eval(bandGap) - energyLevel[itrap];
    tunModel = hTunnelModel[itrap];
  }

  // convert from ScalarT to double
  double field = Sacado::ScalarValue<ScalarT>::eval(fieldAmp);
  double dkbT = Sacado::ScalarValue<ScalarT>::eval(kbT);

  // compute the field-dependent electro-optical energy in [eV]
  double hbarTheta = hbar/q * std::pow(q*q*field*field/(2.0*hbar*m0*effMass), 1.0/3.0);

  // initialize the field enhancement factor
  double fEnhFactor = 0.0;

  // compute the field factor using the temperature approximate expressions in Schenk's paper
  if ((tunModel == "Schenk HighTemp") || (tunModel == "Schenk LowTemp"))
    fEnhFactor = schenkTemperatureApprox(hbarTheta, dkbT, Et, itrap, tunModel);

  // compute the field factor for Tunneling Model = Schenk AsymConstFDOS and Schenk ConstFDOS
  else
  {
    // obtain the phonon energy and related parameters
    double Eph = phononEnergy[itrap];
    double fB = 1.0 / (std::exp(Eph/dkbT) - 1.0);
    double z = 2.0 * hRhysFactor[itrap] * std::sqrt(fB * (fB + 1.0));

    double numerator = schenkFieldFactorNumerator(hbarTheta, dkbT, Et, Eph, z, tunModel);
    double denominator = schenkFieldFactorDenominator(dkbT, Et, Eph, z);

    // compute the field enhancement factor (dosPrefac is canceled due to the division)
    fEnhFactor = numerator / denominator;
  }

  // add the factor of 1
  fEnhFactor = fEnhFactor + 1.0;
  return fEnhFactor;
}


///////////////////////////////////////////////////////////////////////////////
//
//  schenkTemperatureApprox()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::schenkTemperatureApprox(
  const double& hbarTheta, const double& dkbT, const double& Et,
  const int& itrap, const std::string& tunModel)
{
  // the lattice relaxation energy in [eV]
  double ER = hRhysFactor[itrap] * phononEnergy[itrap];

  // the intermediate energy in [eV]
  double EF = std::pow(2.0*ER*dkbT, 2.0) / std::pow(hbarTheta, 3.0);

  // the transition energy in [eV]
  double E0 = 2.0 * std::sqrt(EF*(EF + Et + ER)) - 2.0*EF - ER;

  // std::sqrt(Et - E0) requires a non-negative value of Et - E0
  double Et0 = Et - E0;
  if (Et0 < 0.0)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Cannot use HighTemp Approx or LowTemp Approx when Et-E0 < 0.0 ! Et-E0=" << Et0 << std::endl);

  // std::sqrt(Et*E0) requires a non-negative value of E0 (since Et always > 0)
  if ( (E0 < 0.0) && (tunModel == "Schenk LowTemp") )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Cannot use LowTemp Approx when E0 < 0.0 !");

  double fEnhFactor = 0.0;
  if (tunModel == "Schenk HighTemp")  // Eq.(40) in Schenk's paper
  {
    // the zero-field activation energy
    double Eact0 = (Et - ER) * (Et - ER) / (4.0 * ER);

    // the field-dependent activation energy
    double Eact = (E0 - ER) * (E0 - ER) / (4.0 * ER);

    // temporary variables
    double tmp1 = std::pow(1.0 + 2. * ER * dkbT / (std::pow(hbarTheta, 3./2.) * std::sqrt(Et - E0)), -1./2.);
    double tmp2 = (Eact0 + Et)/dkbT * std::pow(hbarTheta / (Et + ER), 3./2.);
    double tmp3 = std::exp((Eact0 - Eact + Et - E0) / dkbT);
    double tmp4 = std::exp(-4./3. * std::pow((Et - E0) / hbarTheta, 3./2.));

    // the field enhancement factor
    fEnhFactor = tmp1 * tmp2 * tmp3 * tmp4;
  }

  else if (tunModel == "Schenk LowTemp")  // Eq. (49) in Schenk's paper
  {
    // temporary variables
    double Eph = phononEnergy[itrap];
    double tmp1 = std::pow(1.0 + std::pow(hbarTheta, 3./2.) * std::sqrt(Et - E0) / (E0 * Eph), -1./2.);
    double tmp2 = std::pow(hbarTheta, 3./4.) * std::pow(Et-E0, 1./4.) * std::pow(hbarTheta/dkbT, 3./2.) / (2.0*std::sqrt(Et * E0));
    double tmp3 = std::exp(-(Et-E0)/Eph + (Eph-dkbT)/2./Eph + (Et+dkbT/2.)/Eph*std::log(Et/ER) - E0/Eph*std::log(E0/ER));
    double tmp4 = std::exp((Et-E0)/dkbT - 4./3.*std::pow((Et-E0)/hbarTheta, 3./2.));

    // the field enhancement factor
    fEnhFactor = tmp1 * tmp2 * tmp3 * tmp4;
  }

  // otherwise, throw out an error message
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Invalid Tunneling Model name of " << tunModel << std::endl);

  return fEnhFactor;
}


///////////////////////////////////////////////////////////////////////////////
//
//  schenkFieldFactorDenominator()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::schenkFieldFactorDenominator(
  const double& dkbT, const double& Et, const double& Eph, const double& z)
{
  // smaller energy grid for the denominator due to sharp profile in energy
  double dE = Et / double(MAX_NUM_ENERGIES * ENERGY_GRID_RATIO);

  double denIntegrand[MAX_NUM_ENERGIES * ENERGY_GRID_RATIO + 1];
  int numE = MAX_NUM_ENERGIES * ENERGY_GRID_RATIO + 1;

  // compute and store the denominator integrand for E within [Et, 2*Et];
  // for E > 2*Et, the integrand is negligibly small in most cases
  for (int iE = 0; iE < numE; ++iE)
  {
    double E = Et + iE * dE;

    // double dosPrefac = 1/(2*pi^2)*(2*me*q/(hbar*hbar))^1.5; % [eV^(-1.5).m^(-3)];
    // need to be multiplied by dosPrefac to get DOS in [eV^(-1).m^(-3)]
    double zeroFDOS = std::sqrt(E - Et);

    // the number of phonons
    double Nph = E / Eph;

    // compute the modified Bessel function (use its asymptotic form)
    double Iph = 1./std::sqrt(2.*pi) * std::pow(Nph*Nph + z*z, -1./4.) * std::exp(std::sqrt(Nph*Nph + z*z))
               * std::exp(-Nph * std::log(Nph/z + std::sqrt(1.+Nph*Nph/(z*z)) ) );

    // compute the denominator integrand
    denIntegrand[iE] = zeroFDOS * Iph * std::exp(-E/(2.*dkbT));
  }

  // sum up the denominator integrand using the trapezoidal rule
  double denominator = 0.0;
  for (int iE = 0; iE < numE-1; ++iE)
    denominator = denominator + (denIntegrand[iE] + denIntegrand[iE+1]) / 2.0;

  // multiply the sum by the energy spacing
  denominator = denominator * dE;

  return denominator;
}


///////////////////////////////////////////////////////////////////////////////
//
//  schenkFieldFactorNumerator()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::schenkFieldFactorNumerator(
  const double& hbarTheta, const double& dkbT, const double& Et,
  const double& Eph, const double& z, const std::string& tunModel)
{
  double dE = Et / double(MAX_NUM_ENERGIES); // energy grid for the numerator
  double numIntegrand[MAX_NUM_ENERGIES + 1];

  // compute and store the numerator integrand
  for (int iE = 0; iE < MAX_NUM_ENERGIES+1; ++iE)
  {
    double E = iE * dE;
    double y = (Et-E) / hbarTheta;
    if (y < 0.0) y = 0.0;  // y < 0.0 can occur when iE = MAX_NUM_ENERGIES due to numerical round error

    // compute the electro-optical function that appears in the constant field DOS expression,
    double elop_func = 0.0;
    if (tunModel == "Schenk AsymConstFDOS")  // Eq.(27) in Schenk's paper
      // approximate Ai(x) by its asymptotic form at large positive x and approximate Et-E by Et/2
      elop_func = 1./(8.*pi) * (2.*hbarTheta/Et) * std::exp(-4./3.*std::pow(y, 1.5));

    else if (tunModel == "Schenk ConstFDOS")
    {
      // the relative error b.t.w. Airy and asymtotic form is less than 1% for y > 5.0
      if (y > 5.0)
        elop_func = 1./(8.*pi*y) * std::exp(-4./3. * std::pow(y, 1.5));
      else
        elop_func = std::pow(std::fabs(boost::math::airy_ai_prime(y)), 2.0)
                  - y * std::pow(std::fabs(boost::math::airy_ai(y)), 2.0);
    }

    // otherwise, throw out an error message
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Invalid Tunneling Model name of " << tunModel << std::endl);

    // compute the constant field DOS (multiplied by dosPrefac to obtain the correct DOS)
    // double dosPrefac = 1/(2*pi^2)*(2*me*q/(hbar*hbar))^1.5; // [eV^(-1.5).m^(-3)];
    double constFDOS = pi * std::sqrt(hbarTheta) * elop_func;

    // the number of phonons
    double Nph = E / Eph;

    // compute the modified Bessel function (use its asymptotic form)
    double Iph = 1./std::sqrt(2.*pi) * std::pow(Nph*Nph + z*z, -1./4.) * std::exp(std::sqrt(Nph*Nph + z*z))
               * std::exp(-Nph * std::log(Nph/z + std::sqrt(1. + Nph*Nph/(z*z)) ) );

    // compute the numerator integrand
    numIntegrand[iE] = constFDOS * Iph * std::exp(-E/(2.*dkbT));
  }

  // sum up the numerator integrand using the trapzoidal rule
  double numerator = 0.0;
  for (int iE = 0; iE < MAX_NUM_ENERGIES; ++iE)
    numerator = numerator + (numIntegrand[iE] + numIntegrand[iE+1]) / 2.0;

  // multiply the sum by the energy spacing
  numerator = numerator * dE;
  return numerator;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalFieldFactorWithNewDOS()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::evalFieldFactorWithNewDOS(
  const ScalarT& fieldAmp, const ScalarT& kbT, const ScalarT& bandGap,
  const int& itrap, const std::string& carrType, const double& spcLoc)
{
  // obtain the location of heterojunction in the tunneling direction
  double Et = 0.0, hjLoc = 0.0;
  if (carrType == "Electron")
  {
    Et = energyLevel[itrap];       // in unit of [eV]
    hjLoc = eHJLocation[itrap];    // in unit of [X0]
  }
  else if (carrType == "Hole")
  {
    Et = Sacado::ScalarValue<ScalarT>::eval(bandGap) - energyLevel[itrap];
    hjLoc = hHJLocation[itrap];
  }
  else // otherwise, throw out an error message
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Invalid carrier type of " << carrType << std::endl);

  // compute the distance to the heterojunction (HJ) in unit of [meters]
  double xloc = std::fabs(spcLoc - hjLoc) * X0 * 0.01;

  // convert from ScalarT to double
  double field = Sacado::ScalarValue<ScalarT>::eval(fieldAmp);  // [V/m]
  double dkbT = Sacado::ScalarValue<ScalarT>::eval(kbT);        // [eV]

  // do not compute the field factor when the following condition holds
  if ((xloc >= MAX_DISTANCE) && (field <= MIN_ELECTRIC_FIELD) )
  {
    double fEnhFactor = 1.0;
    return fEnhFactor;
  }

  // set the parameters for adaptive trapezoidal integration
  // these parameters are used in fieldFactorIntegrand(...)
  field4Adaptive = field;  // [V/m]
  kbT4Adaptive = dkbT;     // [eV];
  bandGap4Adaptive = Sacado::ScalarValue<ScalarT>::eval(bandGap);
  spcLoc4Adaptive = spcLoc;
  itrap4Adaptive = itrap;
  carrType4Adaptive = carrType;
/*
 {
  // compute the numerator using the uniform trapezoidal rule
  auto time1 = std::chrono::system_clock::now();
  double numerator = fieldFactorWithNewDOSNumerator(fieldAmp, kbT, bandGap, itrap, carrType, spcLoc);
  auto time2 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds1 = time2-time1;
  std::cout << "Numerator uniform trapezoidal: " << std::scientific << numerator;
  std::cout << ", time =  " << std::to_string(elapsed_seconds1.count()) << "s" << std::endl;

  // obtain the phonon energy and related parameters
  double Eph = phononEnergy[itrap];
  double fB = 1.0 / (std::exp(Eph/dkbT) - 1.0);
  double z = 2.0 * hRhysFactor[itrap] * std::sqrt(fB * (fB + 1.0));

  // compute the denominator using the uniform trapezoidal rule
  auto time3 = std::chrono::system_clock::now();
  double denominator = schenkFieldFactorDenominator(dkbT, Et, Eph, z);
  auto time4 = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds2 = time4-time3;
  std::cout << "Denominator uniform trapezoidal: " << std::scientific << denominator;
  std::cout << ", time =  " << std::to_string(elapsed_seconds2.count()) << "s" << std::endl;
 }
*/

  // compute the numerator using the adaptive quadrature rule from github
  double tolerance = 1e-4;
  isNewDOSNum = true;
  double numerator = adaptiveIntegrate(0, Et, tolerance);

  // compute the denominator using the adaptive quadrature rule from github
  isNewDOSNum = false;
  double denominator = adaptiveIntegrate(Et, 2.0*Et, tolerance);

  // compute the field enhancement factor (dosPrefac is canceled due to the division)
  double fEnhFactor = numerator / denominator;

  // add the factor of 1
  fEnhFactor = fEnhFactor + 1.0;

  return fEnhFactor;
}


///////////////////////////////////////////////////////////////////////////////
//
//  fieldFactorWithNewDOSNumerator()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::fieldFactorWithNewDOSNumerator(
  const ScalarT& fieldAmp, const ScalarT& kbT, const ScalarT& bandGap,
  const int& itrap, const std::string& carrType, const double& spcLoc)
{
  // obtain parameters
  double effMass = 0.0, Et = 0.0, offset = 0.0, hjLoc = 0.0;
  if (carrType == "Electron")
  {
    effMass = eEffMass[itrap];     // in unit of [m0]
    Et = energyLevel[itrap];       // in unit of [eV]
    offset = eHJBandOffset[itrap]; // in unit of [eV]
    hjLoc = eHJLocation[itrap];    // in unit of [X0]
  }
  else if (carrType == "Hole")
  {
    effMass = hEffMass[itrap];
    Et = Sacado::ScalarValue<ScalarT>::eval(bandGap) - energyLevel[itrap];
    offset = hHJBandOffset[itrap];
    hjLoc = hHJLocation[itrap];
  }
  // otherwise, throw out an error message
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Invalid carrier type of " << carrType << std::endl);

  // convert from ScalarT to double
  double field = Sacado::ScalarValue<ScalarT>::eval(fieldAmp);  // [V/m]
  double dkbT = Sacado::ScalarValue<ScalarT>::eval(kbT);  // [eV]

  // compute the distance to the heterojunction (HJ) in unit of [meters]
  double xloc = std::fabs(spcLoc - hjLoc) * X0 * 0.01;

  // compute the field-dependent electro-optical energy in [eV]
  double hbarTheta = hbar/q * std::pow(q*q*field*field/(2.0*hbar*m0*effMass), 1.0/3.0);

  // compute relevant parameters used in the Schenk NewDOS model
  double Eloc = offset + field * xloc;   // [eV]
  double alpha = std::pow(2.0*effMass*m0*q*field/(hbar*hbar), 1./3.);  // [1/m]

  // obtain the phonon energy and related parameters
  double Eph = phononEnergy[itrap];
  double fB = 1.0 / (std::exp(Eph/dkbT) - 1.0);
  double z = 2.0 * hRhysFactor[itrap] * std::sqrt(fB * (fB + 1.0));

  // energy grid for the numerator
  double dE = Et / double(MAX_NUM_ENERGIES);  // [eV]
  double numIntegrand[MAX_NUM_ENERGIES + 1];

  // compute and store the numerator integrand
  for (int iE = 0; iE < MAX_NUM_ENERGIES; ++iE)
  {
    double E = iE * dE;  // 0 <= E <= Et, [eV]

    // compute the tunneling dos for Tunneling Model = Schenk NewDOS
    // need to be multiplied by dosPrefac to obtain the DOS in unit of [eV^(-1).m^(-3)]
    // double dosPrefac = 1/(2*pi^2)*(2*me*q/(hbar*hbar))^1.5; // [eV^(-1.5).m^(-3)];
    double tunDOS = calcTunnelDOSForSchenkNewModel(hbarTheta, offset, Et, Eloc, E, dE, xloc, alpha, effMass);

    // the number of phonons
    double Nph = E / Eph;

    // compute the modified Bessel function (use its asymptotic form)
    double Iph = 1./std::sqrt(2.*pi) * std::pow(Nph*Nph + z*z, -1./4.) * std::exp(std::sqrt(Nph*Nph + z*z))
               * std::exp(-Nph * std::log(Nph/z + std::sqrt(1. + Nph*Nph/(z*z)) ) );

    // compute the numerator integrand
    numIntegrand[iE] = tunDOS * Iph * std::exp(-E/(2.*dkbT));
  } // end of loop over energies

  // compute the numerator integrand at E = Et, since E = MAX_NUM_ENERGIES * Et / double(MAX_NUM_ENERGIES)
  // can be slightly larger than Et due to numerical round error, which causes y = (Et-E) / hbarTheta < 0;
  {
    int iE = MAX_NUM_ENERGIES;
    double E = Et;
    double tunDOS = calcTunnelDOSForSchenkNewModel(hbarTheta, offset, Et, Eloc, E, dE, xloc, alpha, effMass);
    double Nph = E / Eph;
    double Iph = 1./std::sqrt(2.*pi) * std::pow(Nph*Nph + z*z, -1./4.) * std::exp(std::sqrt(Nph*Nph + z*z))
               * std::exp(-Nph * std::log(Nph/z + std::sqrt(1. + Nph*Nph/(z*z)) ) );
    numIntegrand[iE] = tunDOS * Iph * std::exp(-E/(2.*dkbT));
  }

  // sum up the numerator integrand using the trapzoidal rule
  double numerator = 0.0;
  for (int iE = 0; iE < MAX_NUM_ENERGIES; ++iE)
    numerator = numerator + (numIntegrand[iE] + numIntegrand[iE+1]) / 2.0;

  // multiply the sum by the energy spacing
  numerator = numerator * dE;
  return numerator;
}


///////////////////////////////////////////////////////////////////////////////
//
//  adaptiveIntegrate() - The adaptive quadrature (a.k.a. trapezoid) rule
//  comes from https://github.com/nblinov/Adaptive-Integrator  which implements
//  https://link.springer.com/article/10.1023%2FA%3A1022318402393,
//  also found in Numerical Recipes. The code from github passes a function that
//  describes the integrand as an argument to adaptiveIntegrate(...), but the
//  function cannot be passed as an argument here for some reason. Now the
//  integrand is directly by the fieldFactorIntegrand(...) routine.
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::adaptiveIntegrate(const double a, const double b, const double tol_)
{
    const double alpha = std::sqrt(2./3.);
    const double beta = 1./std::sqrt(5.);
    const double x1 = .94288241569547971905635175843185720232;
    const double x2 = .64185334234578130578123554132903188354;
    const double x3 = .23638319966214988028222377349205292599;

    double tol, eps;
    eps = std::numeric_limits<double>::epsilon();
    tol = (tol_ < eps) ?  eps : tol_;

    double m, h;
    m = (a+b)/2.; h = (b-a)/2.;

    double y[13] = {this->fieldFactorIntegrand(a), this->fieldFactorIntegrand(m-x1*h),
                    this->fieldFactorIntegrand(m-alpha*h), this->fieldFactorIntegrand(m-x2*h),
                    this->fieldFactorIntegrand(m-beta*h), this->fieldFactorIntegrand(m-x3*h),
                    this->fieldFactorIntegrand(m), this->fieldFactorIntegrand(m+x3*h),
                    this->fieldFactorIntegrand(m+beta*h), this->fieldFactorIntegrand(m+x2*h),
                    this->fieldFactorIntegrand(m+alpha*h), this->fieldFactorIntegrand(m+x1*h),
                    this->fieldFactorIntegrand(b)};

    double fa, fb;
    fa = y[0]; fb = y[12];

    double i1, i2, is;

    i2 = (h/6.)*(y[0] + y[12] + 5.*(y[4] + y[8]));
    i1 = (h/1470.)*(77.*(y[0]+y[12]) + 432.*(y[2]+y[10]) + 625.*(y[4]+y[8]) + 672.*y[6]);
    is = h*(.0158271919734802*(y[0]+y[12]) + .0942738402188500*(y[1]+y[11])
           + .155071987336585*(y[2]+y[10]) + .188821573960182*(y[3]+y[9])
           + .199773405226859*(y[4]+y[8]) + .224926465333340*(y[5]+y[7])
           + .242611071901408*y[6]);

    double erri1, erri2, R;
    erri1 = fabs(i1 - is); erri2 = fabs(i2 - is);
    R = (erri2 != 0.) ? erri1/erri2 : 1.;

    tol = (R > 0. and R < 1.) ? tol_/R : tol_;
    is = fabs(is)*tol/eps;
    if (is == 0.) is = b-a;

    return adaptlobstp(a, b, fa, fb, is);
}


///////////////////////////////////////////////////////////////////////////////
//
//  adaptlobstp()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::adaptlobstp(const double a,
  const double b, const double fa, const double fb, const double is)
{
    const double alpha = std::sqrt(2./3.);
    const double beta = 1./std::sqrt(5.);
    bool terminated = false;

    double m, h;
    m = (a+b)/2.; h = (b-a)/2.;

    double mll, ml, mr, mrr;
    mll = m - alpha*h; ml = m - beta*h; mr = m + beta*h; mrr = m + alpha*h;

    double fmll, fml, fm, fmr, fmrr;
    fmll = this->fieldFactorIntegrand(mll);
    fml = this->fieldFactorIntegrand(ml);
    fm = this->fieldFactorIntegrand(m);
    fmr = this->fieldFactorIntegrand(mr);
    fmrr = this->fieldFactorIntegrand(mrr);

    double i1, i2;
    i2 = (h/6.)*(fa + fb + 5.*(fml+fmr));
    i1 = (h/1470.)*(77*(fa+fb) + 432.*(fmll + fmrr) + 625.*(fml + fmr) + 672.*fm);

    if (is + (i1-i2) == is or mll <= a or b <= mrr)
    {
        if ( (m <= a or b <= m) and !terminated)
        {
            std::cout << "m=" << m << ", a=" << a << ", b=" << b << std::endl;
            std::cerr << "No machine number in the interval. Requested tolerance may not be met.\n";
            terminated = true;
        }
        return i1;
    }
    else
    {
        return adaptlobstp(a,mll,fa,fmll,is)
             + adaptlobstp(mll,ml,fmll,fml,is)
             + adaptlobstp(ml,m,fml,fm,is)
             + adaptlobstp(m,mr,fm,fmr,is)
             + adaptlobstp(mr,mrr,fmr,fmrr,is)
             + adaptlobstp(mrr,b,fmrr,fb,is);
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//  fieldFactorIntegrand() - compute the numerator and denominator integrands
//  of the field factor for Tunneling Model = Schenk NewDOS
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::fieldFactorIntegrand(const double E)
{
  // obtain parameters
  double effMass = 0.0, Et = 0.0, offset = 0.0, hjLoc = 0.0;

  if (carrType4Adaptive == "Electron")
  {
    effMass = eEffMass[itrap4Adaptive];     // in unit of [m0]
    Et = energyLevel[itrap4Adaptive];       // in unit of [eV]
    offset = eHJBandOffset[itrap4Adaptive]; // in unit of [eV]
    hjLoc = eHJLocation[itrap4Adaptive];    // in unit of [X0]
  }
  else if (carrType4Adaptive == "Hole")
  {
    effMass = hEffMass[itrap4Adaptive];
    Et = bandGap4Adaptive - energyLevel[itrap4Adaptive];
    offset = hHJBandOffset[itrap4Adaptive];
    hjLoc = hHJLocation[itrap4Adaptive];
  }
  // otherwise, throw out an error message
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Invalid carrier type of " << carrType4Adaptive << std::endl);

  double tunDOS = 0.0;
  if (isNewDOSNum)  // compute the new dos for the numerator of field factor
  {
    double field = field4Adaptive;  // [V/m]

    // compute the distance to the heterojunction (HJ) in unit of [meters]
    double xloc = std::fabs(spcLoc4Adaptive - hjLoc) * X0 * 0.01;

    // compute the field-dependent electro-optical energy in [eV]
    double hbarTheta = hbar/q * std::pow(q*q*field*field/(2.0*hbar*m0*effMass), 1.0/3.0);

    // compute relevant parameters used in the Schenk NewDOS model
    double Eloc = offset + field * xloc;   // [eV]
    double alpha = std::pow(2.0*effMass*m0*q*field/(hbar*hbar), 1./3.);  // [1/m]

    // energy grid for the numerator
    double dE = Et / double(MAX_NUM_ENERGIES);  // [eV]

    // call the function to compute the new dos
    tunDOS = calcTunnelDOSForSchenkNewModel(hbarTheta, offset, Et, Eloc, E, dE, xloc, alpha, effMass);
  }
  else  // compute the zero-field dos for the denominator of field factor
    tunDOS = std::sqrt(E - Et);

  double dkbT = kbT4Adaptive;  // [eV]

  // obtain the phonon energy and related parameters
  double Eph = phononEnergy[itrap4Adaptive];
  double fB = 1.0 / (std::exp(Eph/dkbT) - 1.0);
  double z = 2.0 * hRhysFactor[itrap4Adaptive] * std::sqrt(fB * (fB + 1.0));
  double Nph = E / Eph;  // the number of phonons

  // compute the modified Bessel function (use its asymptotic form)
  double Iph = 1./std::sqrt(2.*pi) * std::pow(Nph*Nph + z*z, -1./4.) * std::exp(std::sqrt(Nph*Nph + z*z))
               * std::exp(-Nph * std::log(Nph/z + std::sqrt(1. + Nph*Nph/(z*z)) ) );

  // compute the integrand
  double integrand = tunDOS * Iph * std::exp(-E/(2.*dkbT));
  return integrand;
}


///////////////////////////////////////////////////////////////////////////////
//
//  calcTunnelDOSForSchenkNewModel()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::calcTunnelDOSForSchenkNewModel(
const double& hbarTheta, const double& bOffset, const double& Et, const double& Eloc,
const double& E, const double& /* dEx */, const double& xloc, const double& alpha, const double& me)
{
  double tunDOS = 0.0;

  // use the asymptotic constant-field DOS since E << Et (far from band edge)
  if (E <= Et - Eloc)  // this part executes only when Et > Eloc
  {
    double y = (Et - E) / hbarTheta;
    double elop_func = 1./(8.*pi) * (2.*hbarTheta/Et) * std::exp(-4./3.*std::pow(y, 1.5));
    tunDOS = pi * std::sqrt(hbarTheta) * elop_func;

    // maybe, should set tunDOS = 0.0 here since there is no carrier for E < Et-Eloc
  }
  else
  {
    double Enew = Eloc - Et + E;  // energy upper bound, and the lower bound is 0
    double minBeta = (bOffset - Enew) / hbarTheta;  // unitless
    double maxBeta = bOffset / hbarTheta;  // unitless

    // this condition may occur for small fields
    if ((minBeta > MAX_BETA) || (maxBeta > MAX_BETA))
    {
      if (xloc <= MAX_DISTANCE)  // use the step-barrier DOS
      {
        Enew = bOffset -  Et + E;  // reset Enew for step barrier
        if (Enew > 0.0)
        {
          //tunDOS = calcDOSForStepBarrier(Enew, dEx, bOffset, xloc, me);
          tunDOS = calcDOSForStepBarrierGaussQR(Enew, bOffset, xloc, me);
        }
        else
          tunDOS = 0.0;
      }
      else  // use the asymptotic constant-field DOS (use the asymptotic form
      {     // of Airy(x) at large positive x and approximate Et-E by Et/2)
        double  y = (Et-E) / hbarTheta;
        double elop_func = 1./(8.*pi) * (2.*hbarTheta/Et) * std::exp(-4./3.*std::pow(y, 1.5));
        tunDOS = pi * std::sqrt(hbarTheta) * elop_func;
      }
    }

    else  // when all beta > MAX_BETA
    {
      if (xloc <= MAX_DISTANCE)  // use the new DOS calculation
      {
        //tunDOS = calcDOSForLinPotWithOffset(hbarTheta, bOffset, Enew, dEx, xloc, alpha);
        tunDOS = calcDOSForLinPotWithOffsetGaussQR(hbarTheta, bOffset, Enew, xloc, alpha);
      }
      else  // use the asymptotic constant-field DOS for distances far from the HJ
      {
        double y = (Et-E) / hbarTheta;  // note y is always >= 0
        double elop_func = 1./(8.*pi) * (2.*hbarTheta/Et) * std::exp(-4./3.*std::pow(y, 1.5));
        tunDOS = pi * std::sqrt(hbarTheta) * elop_func;
      }
    }
  }
  return tunDOS;  // tunDOS in the unit of [eV^0.5]
}


///////////////////////////////////////////////////////////////////////////////
//
//  calcDOSForStepBarrier()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::calcDOSForStepBarrier(const double& Enew,
const double& dEx, const double& bOffset, const double& xloc, const double& me)
{
  // This new implementation of the uniform trapezoid rule does not reduce Charon's
  // simulation due to the negligible ime for the DOS calculation, when compared to
  // the original implementation where the integrand is computed and stored first
  // and summed up later.
  const int nEx = (int)std::floor(Enew/dEx);
  double tunDOS = 0.0;

  double firstIntegrand = 0.0; // equal to zero for Ex = 0;
  double lastIntegrand = 0.0;  // default
  double Ex = 0.0, nv = 0.0;

  // compute the last integrand
  if (nEx > 0)
  {
    Ex = nEx * dEx;  // [eV]
    nv = std::sqrt(2.0*me*m0*(bOffset-Ex)*q / (hbar*hbar));  // [1/m]
    lastIntegrand = std::sqrt(Ex) * std::exp(-2.*nv*xloc);
  }

  // compute and sum the integrand for iEx = 1 to (nEx-1) using trapzoidal rule
  for (int iEx = 1; iEx < nEx; iEx++)
  {
    Ex = iEx * dEx;  // [eV]
    nv = std::sqrt(2.0*me*m0*(bOffset-Ex)*q / (hbar*hbar));  // [1/m]
    tunDOS += std::sqrt(Ex) * std::exp(-2.*nv*xloc);
  }

  // sum in the first and last integrands and multiply dEx
  tunDOS = (tunDOS + (firstIntegrand + lastIntegrand) * 0.5) * dEx;

  // compute the DOS integrand for Ex = Enew
  double extraIntegrand = 0.0;
  {
    Ex = Enew;
    nv = std::sqrt(2.0*me*m0*(bOffset-Ex)*q / (hbar*hbar));  // [1/m]
    extraIntegrand = std::sqrt(Ex) * std::exp(-2.*nv*xloc);
  }

  // add the last interval contribution if any
  double lastdEx = Enew - nEx*dEx;  // note: lastdEx could be 0 if Enew = nEx*dEx
  tunDOS = tunDOS + (lastIntegrand + extraIntegrand) * lastdEx * 0.5;

  // divide by bOffset
  tunDOS = tunDOS / bOffset;   // tunDOS in the unit of [eV^0.5]
  return tunDOS;
}


///////////////////////////////////////////////////////////////////////////////
//
//  calcDOSForLinPotWithOffset()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::calcDOSForLinPotWithOffset(
const double& hbarTheta, const double& bOffset, const double& Enew,
const double& dEx, const double& xloc, const double& alpha)
{
  const int nEx = (int)std::floor(Enew/dEx);
  double tunDOS = 0.0;

  double firstIntegrand = 0.0; // equal to zero for Ex = 0;
  double lastIntegrand = 0.0;  // default
  double Ex = 0.0, beta = 0.0, y = 0.0, airy1 = 0.0, airy2 = 0.0, airyp = 0.0;

  // compute the last integrand
  if (nEx > 0)
  {
    Ex = nEx * dEx;  // [eV]
    beta = (bOffset - Ex) / hbarTheta;
    y = alpha * xloc + beta;
    airy1 = std::pow(std::fabs(boost::math::airy_ai(y)), 2.0);
    airy2 = std::pow(std::fabs(boost::math::airy_ai(beta)), 2.0);
    airyp = std::pow(std::fabs(boost::math::airy_ai_prime(beta)), 2.0);
    lastIntegrand = 1./std::sqrt(Ex) * airy1 / (airy2 + hbarTheta/Ex * airyp);
  }

  // compute and sum the integrand for iEx = 1 to (nEx-1) using trapzoidal rule
  for (int iEx = 1; iEx < nEx; iEx++)
  {
    Ex = iEx * dEx;  // [eV]
    beta = (bOffset - Ex) / hbarTheta;
    y = alpha * xloc + beta;
    airy1 = std::pow(std::fabs(boost::math::airy_ai(y)), 2.0);
    airy2 = std::pow(std::fabs(boost::math::airy_ai(beta)), 2.0);
    airyp = std::pow(std::fabs(boost::math::airy_ai_prime(beta)), 2.0);
    tunDOS += 1./std::sqrt(Ex) * airy1 / (airy2 + hbarTheta/Ex * airyp);
  }

  // sum in the first and last integrands and multiply dEx
  tunDOS = (tunDOS + (firstIntegrand + lastIntegrand) * 0.5) * dEx;

  // compute the DOS integrand for Ex = Enew
  double extraIntegrand = 0.0;
  {
    Ex = Enew;
    beta = (bOffset - Ex) / hbarTheta;
    y = alpha * xloc + beta;
    airy1 = std::pow(std::fabs(boost::math::airy_ai(y)), 2.0);
    airy2 = std::pow(std::fabs(boost::math::airy_ai(beta)), 2.0);
    airyp = std::pow(std::fabs(boost::math::airy_ai_prime(beta)), 2.0);
    extraIntegrand = 1./std::sqrt(Ex) * airy1 / (airy2 + hbarTheta/Ex * airyp);
  }

  // add the last interval contribution if any
  double lastdEx = Enew - nEx*dEx;  // note: lastdEx could be 0 if Enew = nEx*dEx
  tunDOS = tunDOS + (lastIntegrand + extraIntegrand) * lastdEx * 0.5;
  return tunDOS;
}


///////////////////////////////////////////////////////////////////////////////
//
//  calcDOSForStepBarrierGaussQR()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::calcDOSForStepBarrierGaussQR(const double& Enew,
const double& bOffset, const double& xloc, const double& me)
{
  double tunDOS = 0.0;
  double p1 = 2.*xloc *std::sqrt(2.0*me*m0*Enew*q /(hbar*hbar)); // unitless
  double p2 = bOffset / Enew;  // unitless
  int order = order20_dos;

  // perform integration using the Gauss-Jacobi quadrature rule with alpha=0.5, beta=0 for the Jacobi polynomials
  for (int ii = 0; ii < order; ii++)
  {
    double xi = abscissas20_dos[ii];
    double wi = weights20_dos[ii];
    double fi = std::exp(-p1*std::sqrt(p2-xi));
    tunDOS += wi * fi;
  }

  // multiply the necessary factor due to change of variables
  tunDOS = tunDOS * std::pow(Enew,1.5) / bOffset;   // tunDOS in the unit of [eV^0.5]
  return tunDOS;
}


///////////////////////////////////////////////////////////////////////////////
//
//  calcDOSForLinPotWithOffsetGaussQR()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double RecombRate_TrapSRH<EvalT, Traits>::calcDOSForLinPotWithOffsetGaussQR(
const double& hbarTheta, const double& bOffset, const double& Enew,
const double& xloc, const double& alpha)
{
  double tunDOS = 0.0;
  double p1 = alpha * xloc;     // unitless
  double p2 = Enew / hbarTheta; // unitless
  double p3 = bOffset / Enew;   // unitless
  int order = order20_dos;

  // perform integration using the Gauss-Jacobi quadrature rule with alpha=0.5, beta=0 for the Jacobi polynomials
  for (int ii = 0; ii < order; ii++)
  {
    double xi = abscissas20_dos[ii];
    double wi = weights20_dos[ii];
    double arg = p2 * (p3-xi);
    double airy1 = std::pow(std::fabs(boost::math::airy_ai(p1+arg)), 2.0);
    double airy2 = std::pow(std::fabs(boost::math::airy_ai(arg)), 2.0);
    double airyp = std::pow(std::fabs(boost::math::airy_ai_prime(arg)), 2.0);
    double fi = airy1 / (xi*airy2 + airyp/p2);
    tunDOS += wi * fi;
  }

  // multiply the necessary factor due to change of variables
  tunDOS = tunDOS * std::sqrt(Enew);  // tunDOS in the unit of [eV^0.5]
  return tunDOS;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
RecombRate_TrapSRH<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<ParameterList> p = Teuchos::rcp(new ParameterList);
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

  ParameterList& trapPL = p->sublist("Trap SRH ParameterList", false, "Sublist defining the Trap SRH ParameterList");

  for (int i = 0; i < MAX_NUM_TRAPS; i++)
  {
    std::stringstream ss;
    ss << i;
    string subListName("Trap " + ss.str());
    trapPL.sublist(subListName, false, "Sublist defining the parameters for one type of trap");

    // trap related parameters
    trapPL.sublist(subListName).set<double>("Energy Level", 0.0, "Trap energy level measured from the conduction band in [eV]");
    trapPL.sublist(subListName).set<double>("Trap Density", 0.0, "Trap density in [cm^-3]");
    trapPL.sublist(subListName).set<string>("Trap Type", "", "Either Acceptor (0 if unoccupied and -1 if occupied, electron capture) or Donor (0 if unoccupied and +1 if occupied, hole capture)");
    trapPL.sublist(subListName).set<string>("Spatial Profile", "Uniform", "Spatial profile for the trap distribution in space, currently only Uniform");
    trapPL.sublist(subListName).set<double>("X Min", 0., "X min for the spatial box containing trap i");
    trapPL.sublist(subListName).set<double>("X Max", 0., "X max for the spatial box containing trap i");
    trapPL.sublist(subListName).set<double>("Y Min", 0., "Y min for the spatial box containing trap i");
    trapPL.sublist(subListName).set<double>("Y Max", 0., "Y max for the spatial box containing trap i");
    trapPL.sublist(subListName).set<double>("Z Min", 0., "Z min for the spatial box containing trap i");
    trapPL.sublist(subListName).set<double>("Z Max", 0., "Z max for the spatial box containing trap i");

    // phonon parameters
    trapPL.sublist(subListName).set<double>("Phonon Energy", 0.0, "Phonon energy in [eV] for trap i");
    trapPL.sublist(subListName).set<double>("Huang-Rhys Factor", 0.0, "Huang-Rhys Factor [unitless] for trap i");

    // field-independent lifetime parameters
    trapPL.sublist(subListName).set<double>("Electron Lifetime", 0., "Field-independent electron lifetime in [s] for trap i");
    trapPL.sublist(subListName).set<double>("Electron Cross Section", 0.0, "Electron capture cross section in [cm^2]");
    trapPL.sublist(subListName).set<double>("Hole Lifetime", 0., "Field-independent hole lifetime in [s] for trap i");
    trapPL.sublist(subListName).set<double>("Hole Cross Section", 0.0, "Hole capture cross section in [cm^2]");

    // Schenk band-to-trap tunneling related parameters
    trapPL.sublist(subListName).set<string>("Electron Tunneling Model", "None", "None, Schenk HighTemp, Schenk LowTemp, Schenk AsymConstFDOS, Schenk ConstFDOS, or Schenk NewDOS");
    trapPL.sublist(subListName).set<string>("Electron Tunneling Direction", "X", "Required for the Schenk NewDOS model: X, Y, or Z");
    trapPL.sublist(subListName).set<double>("Electron HJ Location", 0.0, "Required for the Schenk NewDOS model, in [um]");
    trapPL.sublist(subListName).set<double>("Electron HJ Band Offset", 0.0, "Required for the Schenk NewDOS model, in [eV]");
    trapPL.sublist(subListName).set<double>("Electron Effective Mass", 0.0, "Electron effective mass for the Schenk's model in [m0]");

    trapPL.sublist(subListName).set<string>("Hole Tunneling Model", "None", "None, Schenk HighTemp, Schenk LowTemp, Schenk AsymConstFDOS, Schenk ConstFDOS, or Schenk NewDOS");
    trapPL.sublist(subListName).set<string>("Hole Tunneling Direction", "X", "Required for the Schenk NewDOS model: X, Y, or Z");
    trapPL.sublist(subListName).set<double>("Hole HJ Location", 0.0, "Required for the Schenk NewDOS model, in [um]");
    trapPL.sublist(subListName).set<double>("Hole HJ Band Offset", 0.0, "Required for the Schenk NewDOS model, in [eV]");
    trapPL.sublist(subListName).set<double>("Hole Effective Mass", 0.0, "Better to use the light hole effective mass for the Schenk's model in [m0]");
  }

  return p;
}


///////////////////////////////////////////////////////////////////////////////
//
//  initTrapSRHParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void RecombRate_TrapSRH<EvalT, Traits>::initTrapSRHParams(
const std::string& matName, const Teuchos::ParameterList& trapSRHParamList)
{
  using Teuchos::ParameterList;
  using std::string;

  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();
  double eMass = matProperty.getPropertyValue(matName, "Electron Effective Mass");
  double hMass = matProperty.getPropertyValue(matName, "Hole Effective Mass");

  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  double q = phyConst.q;    // elemental charge in [C]
  double m0 = phyConst.m0;  // electron mass in [kg]
  double velPrefactor = std::sqrt(3.0*kb*q*300./m0) * 100.;  // [cm/s]

  double xmin = -1e20, xmax = 1e20;
  double ymin = -1e20, ymax = 1e20;
  double zmin = -1e20, zmax = 1e20;
  string spaceProfile = "Uniform";

  for (ParameterList::ConstIterator it = trapSRHParamList.begin(); it != trapSRHParamList.end(); ++it)
  {
    const string key = it->first;
    const Teuchos::ParameterEntry& entry = it->second;
    const ParameterList& trapPList = Teuchos::getValue<ParameterList>(entry);

    // trap related parameters (required)
    energyLevel.push_back(trapPList.get<double>("Energy Level"));  // [eV]
    double Nt = trapPList.get<double>("Trap Density");  // [cm^-3]
    trapDensity.push_back(Nt);

    string tType = trapPList.get<string>("Trap Type");
    if ((tType != "Acceptor") && (tType != "Donor"))  // either Acceptor or Donor
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Either Acceptor or Donor must be specified for Trap Type!");
    trapType.push_back(tType);

    // trap related parameters (optional)
    if (trapPList.isParameter("Spatial Profile"))
      spaceProfile = trapPList.get<string>("Spatial Profile");
    trapSpcProfile.push_back(spaceProfile);
    if (trapPList.isParameter("X Min"))  xmin = trapPList.get<double>("X Min");
    if (trapPList.isParameter("X Max"))  xmax = trapPList.get<double>("X Max");
    if (trapPList.isParameter("Y Min"))  ymin = trapPList.get<double>("Y Min");
    if (trapPList.isParameter("Y Max"))  ymax = trapPList.get<double>("Y Max");
    if (trapPList.isParameter("Z Min"))  zmin = trapPList.get<double>("Z Min");
    if (trapPList.isParameter("Z Max"))  zmax = trapPList.get<double>("Z Max");
    minXCoord.push_back(xmin);
    maxXCoord.push_back(xmax);
    minYCoord.push_back(ymin);
    maxYCoord.push_back(ymax);
    minZCoord.push_back(zmin);
    maxZCoord.push_back(zmax);

    // field-independent lifetime parameters (electrons)
    bool eTime = false;
    double taun0 = 0.0;
    if (trapPList.isParameter("Electron Lifetime"))
    {
      taun0 = trapPList.get<double>("Electron Lifetime");  // [s]
      eTime = true;
    }
    else if (trapPList.isParameter("Electron Cross Section"))
    {
      double eXsec = trapPList.get<double>("Electron Cross Section");  // [cm^2]
      double eVel = velPrefactor / std::sqrt(eMass);  // [cm/s]
      taun0 = 1.0 / (eXsec * eVel * Nt);  // [s]
    }
    else if (trapPList.isParameter("Electron Lifetime") && trapPList.isParameter("Electron Cross Section"))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Electron Lifetime and Electron Cross Section cannot be specified simultaneously !");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Either Electron Lifetime or Electron Cross Section must be specified !");
    eLifeTime.push_back(taun0);
    eTimeGiven.push_back(eTime);

    // field-independent lifetime parameters (holes)
    bool hTime = false;
    double taup0 = 0.0;
    if (trapPList.isParameter("Hole Lifetime"))
    {
      taup0 = trapPList.get<double>("Hole Lifetime");  // [s]
      hTime = true;
    }
    else if (trapPList.isParameter("Hole Cross Section"))
    {
      double hXsec = trapPList.get<double>("Hole Cross Section");  // [cm^2]
      double hVel = velPrefactor / std::sqrt(hMass);  // [cm/s]
      taup0 = 1.0 / (hXsec * hVel * Nt);  // [s]
    }
    else if (trapPList.isParameter("Hole Lifetime") && trapPList.isParameter("Hole Cross Section"))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Hole Lifetime and Hole Cross Section cannot be specified simultaneously !");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! Either Hole Lifetime or Hole Cross Section must be specified !");
    hLifeTime.push_back(taup0);
    hTimeGiven.push_back(hTime);

    // default tunneling model
    string eTunModel = "None", hTunModel = "None";

    // obtain specified tunneling model
    if (trapPList.isParameter("Electron Tunneling Model"))
      eTunModel = trapPList.get<string>("Electron Tunneling Model");
    if (trapPList.isParameter("Hole Tunneling Model"))
      hTunModel = trapPList.get<string>("Hole Tunneling Model");
    eTunnelModel.push_back(eTunModel);
    hTunnelModel.push_back(hTunModel);

    // get phonon parameters (must be specified for all Schenk models)
    double phE = 0.0, hRhys = 0.0;
    if ( (eTunModel.find("Schenk") != std::string::npos) ||
         (hTunModel.find("Schenk") != std::string::npos) )  // find "Schenk"
    {
      phE = trapPList.get<double>("Phonon Energy");  // [eV]
      hRhys = trapPList.get<double>("Huang-Rhys Factor");  // [1]
    }
    phononEnergy.push_back(phE);
    hRhysFactor.push_back(hRhys);

    // get the effective mass (needed for all Schenk models)
    if (trapPList.isParameter("Electron Effective Mass"))
      eMass = trapPList.get<double>("Electron Effective Mass");
    eEffMass.push_back(eMass);
    if (trapPList.isParameter("Hole Effective Mass"))
      hMass = trapPList.get<double>("Hole Effective Mass");
    hEffMass.push_back(hMass);

    // get the Schenk NewDOS model parameters (electrons)
    double eHJLoc = 0.0, eHJBand = 0.0;
    string eTunDir = "X";
    if (eTunModel == "Schenk NewDOS")
    {
      eHJLoc = trapPList.get<double>("Electron HJ Location");  // [um]
      eHJBand = trapPList.get<double>("Electron HJ Band Offset");  // [eV]
      eTunDir = trapPList.get<string>("Electron Tunneling Direction");
    }
    eHJLocation.push_back(eHJLoc);
    eHJBandOffset.push_back(eHJBand);
    eTunnelDir.push_back(eTunDir);

    // get the Schenk NewDOS model parameters (holes)
    double hHJLoc = 0.0, hHJBand = 0.0;
    string hTunDir = "X";
    if (hTunModel == "Schenk NewDOS")
    {
      hHJLoc = trapPList.get<double>("Hole HJ Location");  // [um]
      hHJBand = trapPList.get<double>("Hole HJ Band Offset");  // [eV]
      hTunDir = trapPList.get<string>("Hole Tunneling Direction");
    }
    hHJLocation.push_back(hHJLoc);
    hHJBandOffset.push_back(hHJBand);
    hTunnelDir.push_back(hTunDir);
  }  // end of for loop
}

}

#endif
