
#ifndef CHARON_DDLATTICEBC_OHMICCONTACT_IMPL_HPP
#define CHARON_DDLATTICEBC_OHMICCONTACT_IMPL_HPP

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
#include "Charon_Util.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

/*
The Ohmic BC calculations here for electrons and holes are applicable to
the DD+Lattice, DD+Ion, and DD+Lattice+Ion formulations. The potential
here is the intrinsic Fermi potential, i.e., -q*\phi = Ei.
*/

namespace charon {

// select a physical solution for carrier density at contact
double select_physical_sol(std::vector<double>& sol, double dop) {
  std::vector<double> phy_sol;
  // only accept positive roots
  for(size_t i=0; i<sol.size(); i++)
    if(sol[i] > 0) phy_sol.push_back(sol[i]);
  // only accept roots below dop level
  for(size_t i=0; i<phy_sol.size(); i++)
    if(phy_sol[i] > dop) phy_sol.erase(phy_sol.begin() + i);
  // select root closest to dop
  double phy_solution = 0;
  double d = 1e30;
  for(size_t i=0; i<phy_sol.size(); i++) {
    if((dop - phy_sol[i]) < d) {
      phy_solution = phy_sol[i];
      d = dop - phy_sol[i];
    }
  }
  assert(phy_solution > 0);
  return phy_solution;
}

double evaluateIonizEnFromFile(
               const std::vector<double>& dop_table,
               const std::map<double,double>& IonizEn,
               double dop) {
  // dop_table is assumed to be ordered from low to high doping
  // find the two adjacent values in the dop_table bounding dop
  std::pair<double,double> dop_bounds;
  charon::findDopingPoints(dop_table,dop,dop_bounds);
  // compute ionization energy
  return charon::interpolateIonizEn(IonizEn,dop_bounds,dop);
}


// compute carrier density with MB Statistics at an Ohmic contact
// when incomplete ionization model is active
double compute_MB_carrier_dens (
     int q, const Teuchos::ParameterList& incmpl_ioniz_param,
     double kbT, double Nc, double Nv, double ni, double Na,
     double Nd, double Nion,
     double dens_sc) {
  double phy_sol = 0;
  assert(q != 0);

  if(q < 0) { // electron density in n-type semic
    assert(incmpl_ioniz_param.sublist("Donor").numParams() > 0);
    const Teuchos::ParameterList& params =
      incmpl_ioniz_param.sublist("Donor");
    std::string approx_type = params.get<std::string>("Approximation");
    double Ed = 0.;
    if(params.isSublist("DonIncmplIonizData")) {
      Teuchos::RCP<std::vector<double>> conc =
        params.sublist("DonIncmplIonizData").
        get<Teuchos::RCP<std::vector<double> > >("donConc");
      Teuchos::RCP<std::map<double,double> > ioniz_en =
        params.sublist("DonIncmplIonizData").
        get<Teuchos::RCP<std::map<double,double> > >("donIonizEn");
      Ed = evaluateIonizEnFromFile(*conc,*ioniz_en,Nd*dens_sc);
    } else {
      Ed = params.get<double>("Ionization Energy");
    }
    double gD = params.get<double>("Degeneracy Factor");
    double n1 = Nc*exp(-Ed/kbT);

    if(approx_type == "I") { // acceptor fully ionized
      double a3 = gD;
      double a2 = n1 + gD*(Na-Nion);
      double a1 = n1*Na - n1*Nd - ni*ni*gD - n1*Nion;
      double a0 = -n1*ni*ni;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol2;
      sol2.push_back(r1_re); sol2.push_back(r2_re); sol2.push_back(r3_re);
      phy_sol = select_physical_sol(sol2, Nd);
    } else if(approx_type == "II") { // acceptor neglected
      double a3 = gD;
      double a2 = n1 - gD*Nion;
      double a1 = -(ni*ni*gD + n1*Nd + n1*Nion);
      double a0 = -n1*ni*ni;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol1;
      sol1.push_back(r1_re); sol1.push_back(r2_re); sol1.push_back(r3_re);
      phy_sol = select_physical_sol(sol1, Nd);
    } else if(approx_type == "III") { // full model
      // make sure that Acceptor Incomplete Ionization model is active
      if(incmpl_ioniz_param.sublist("Acceptor").numParams() == 0)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
             "Error! Incomplete Ionization Model for Acceptor MUST be defined \
             when Approximation type III is used for Donor Incomplete Ionization");
      const Teuchos::ParameterList& params1 =
        incmpl_ioniz_param.sublist("Acceptor");
      double Ea = 0.;
      if(params1.isSublist("AccIncmplIonizData")) {
        Teuchos::RCP<std::vector<double>> conc1 =
          params1.sublist("AccIncmplIonizData").
          get<Teuchos::RCP<std::vector<double>>>("accConc");
        Teuchos::RCP<std::map<double,double>> ioniz_en1 =
          params1.sublist("AccIncmplIonizData").
          get<Teuchos::RCP<std::map<double,double>>>("accIonizEn");
        Ea = evaluateIonizEnFromFile(*conc1,*ioniz_en1,Na*dens_sc);
      } else {
        Ea = params1.get<double>("Ionization Energy");
      }
      double gA = params1.get<double>("Degeneracy Factor");
      double p1 = Nv*exp(-Ea/kbT);
      double a4 = gD*p1;
      double a3 = gA*gD*ni*ni + n1*p1 + gD*Na*p1 - gD*Nion*p1;
      double a2 = gA*n1*ni*ni +  n1*p1*Na - n1*p1*Nd -
                  gD*p1*ni*ni - gA*gD*ni*ni*Nion - n1*Nion*p1;
      double a1 = -gA*n1*Nd*ni*ni - gA*gD*ni*ni*ni*ni -
                  n1*p1*ni*ni - gA*n1*ni*ni*Nion;
      double a0 = -gA*n1*ni*ni*ni*ni;
      // solve for solution analytically
      double r1=0, r2=0, r3=0, r4=0;
      charon::quarticsolve_salzer(a4,a3,a2,a1,a0,r1,r2,r3,r4);
      std::vector<double> sol3;
      sol3.push_back(r1); sol3.push_back(r2);
      sol3.push_back(r3); sol3.push_back(r4);
      phy_sol = select_physical_sol(sol3, Nd);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
           "Error! Incomplete Ionization Donor Model MUST specify an \
           Approximation of type I, II or III");
    }
  } else { // hole density in p-type semic
    assert(incmpl_ioniz_param.sublist("Acceptor").numParams() > 0);
    const Teuchos::ParameterList& params =
      incmpl_ioniz_param.sublist("Acceptor");
    std::string approx_type = params.get<std::string>("Approximation");
    double Ea = 0.;
    if(params.isSublist("AccIncmplIonizData")) {
      Teuchos::RCP<std::vector<double> > conc =
        params.sublist("AccIncmplIonizData").
        get<Teuchos::RCP<std::vector<double> > >("accConc");
      Teuchos::RCP<std::map<double,double> > ioniz_en =
        params.sublist("AccIncmplIonizData").
        get<Teuchos::RCP<std::map<double,double> > >("accIonizEn");
      Ea = evaluateIonizEnFromFile(*conc,*ioniz_en,Na*dens_sc);
    } else {
      Ea = params.get<double>("Ionization Energy");
    }
    double gA = params.get<double>("Degeneracy Factor");
    double p1 = Nv*exp(-Ea/kbT);

    if(approx_type == "I") { // donor fully ionized
      double a3 = gA;
      double a2 = p1 + gA*Nd + gA*Nion;
      double a1 = -gA*ni*ni - p1*(Na-Nd-Nion);
      double a0 = -p1*ni*ni;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol2;
      sol2.push_back(r1_re); sol2.push_back(r2_re); sol2.push_back(r3_re);
      phy_sol = select_physical_sol(sol2, Na);
    } else if(approx_type == "II") { // donor neglected
      double a3 = gA;
      double a2 = p1 + gA*Nion;
      double a1 = -gA*ni*ni - p1*(Na-Nion);
      double a0 = -p1*ni*ni;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol1;
      sol1.push_back(r1_re); sol1.push_back(r2_re); sol1.push_back(r3_re);
      phy_sol = select_physical_sol(sol1, Na);
    } else if(approx_type == "III") { // full model
      // make sure that Donor Incomplete Ionization model is active
      if(incmpl_ioniz_param.sublist("Donor").numParams() == 0)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Error! Incomplete Ionization Model for Donor MUST be defined \
             when Approximation type III is used for Acceptor Incomplete Ionization");
      const Teuchos::ParameterList& params1 =
        incmpl_ioniz_param.sublist("Donor");
      double Ed = 0.;
      if(params1.isSublist("DonIncmplIonizData")) {
        Teuchos::RCP<std::vector<double>> conc1 =
          params1.sublist("DonIncmplIonizData").
          get<Teuchos::RCP<std::vector<double> > >("donConc");
        Teuchos::RCP<std::map<double,double> > ioniz_en1 =
          params1.sublist("DonIncmplIonizData").
          get<Teuchos::RCP<std::map<double,double> > >("donIonizEn");
        Ed = evaluateIonizEnFromFile(*conc1,*ioniz_en1,Nd*dens_sc);
      } else {
        Ed = params.get<double>("Ionization Energy");
      }
      double gD = params1.get<double>("Degeneracy Factor");
      double n1 = Nc*exp(-Ed/kbT);
      double a4 = gA*n1;
      double a3 = gA*gD*ni*ni + n1*p1 + gA*n1*(Nd+Nion);
      double a2 = n1*ni*ni*(gD-gA) - n1*p1*(Na-Nd-Nion);
      double a1 = -gA*gD*ni*ni*ni*ni -n1*p1*ni*ni -
                  gD*p1*ni*ni*(Na-Nion);
      double a0 = -gD*p1*ni*ni*ni*ni;
      // solve for solution analytically
      double r1=0, r2=0, r3=0, r4=0;
      charon::quarticsolve_salzer(a4,a3,a2,a1,a0,r1,r2,r3,r4);
      std::vector<double> sol3;
      sol3.push_back(r1); sol3.push_back(r2);
      sol3.push_back(r3); sol3.push_back(r4);
      phy_sol = select_physical_sol(sol3, Na);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
           "Error! Incomplete Ionization Acceptor Model MUST specify an \
           Approximation of type I, II or III");
    }
  }
  return phy_sol;
}



// compute carrier density with FD Statistics at an Ohmic contact
// when incomplete ionization model is active
double compute_FD_carrier_dens (
                int q, const Teuchos::ParameterList& incmpl_ioniz_param,
                double kbT, double Nc, double Nv, double ni, double Na,
                double Nd, double Nion, double gamma_n, double gamma_p,
                double dens_sc) {
  double phy_sol = 0;
  assert(q != 0);

  if(q < 0) { // electron density in n-type semic
    assert(incmpl_ioniz_param.sublist("Donor").numParams() > 0);
    const Teuchos::ParameterList& params =
      incmpl_ioniz_param.sublist("Donor");
    std::string approx_type = params.get<std::string>("Approximation");
    double Ed = 0.;
    if(params.isSublist("DonIncmplIoniz File")) {
      Teuchos::RCP<std::vector<double>> conc =
        params.sublist("DonIncmplIonizData").
        get<Teuchos::RCP<std::vector<double>>>("donConc");
      Teuchos::RCP<std::map<double,double>> ioniz_en =
        params.sublist("DonIncmplIonizData").
        get<Teuchos::RCP<std::map<double,double>>>("donIonizEn");
      Ed = evaluateIonizEnFromFile(*conc,*ioniz_en,Nd*dens_sc);
    } else {
      Ed = params.get<double>("Ionization Energy");
    }
    double gD = params.get<double>("Degeneracy Factor");
    double n1 = Nc*exp(-Ed/kbT);

    if(approx_type == "I") { // acceptor fully ionized
      const double tt1 = gD*(Na-Nion) + gamma_n*n1;
      const double tt2 = tt1*tt1 + 4*gD*gamma_n*n1*(Nd-Na+Nion);
      assert(tt2 >= 0);
      phy_sol = -tt1 + sqrt(tt2);
      assert(phy_sol >= 0);
      phy_sol /= 2*gD;
    } else if(approx_type == "II") { // acceptor neglected
      const double tt1 = -gD*Nion + gamma_n*n1;
      const double tt2 = tt1*tt1 + 4*gD*gamma_n*n1*(Nd+Nion);
      assert(tt2 >= 0);
      phy_sol = -tt1 + sqrt(tt2);
      assert(phy_sol >= 0);
      phy_sol /= 2*gD;
    } else if(approx_type == "III") { // full model
      // make sure that Acceptor Incomplete Ionization model is active
      if(incmpl_ioniz_param.sublist("Acceptor").numParams() == 0)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
             "Error! Incomplete Ionization Model for Acceptor MUST be defined \
             when Approximation type III is used for Donor Incomplete Ionization");
      const Teuchos::ParameterList& params1 =
        incmpl_ioniz_param.sublist("Acceptor");
      double Ea = 0.;
      if(params1.isSublist("AccIncmplIonizData")) {
        Teuchos::RCP<std::vector<double>> conc1 =
          params1.sublist("AccIncmplIonizData").
          get<Teuchos::RCP<std::vector<double>>>("accConc");
        Teuchos::RCP<std::map<double,double>> ioniz_en1 =
          params1.sublist("AccIncmplIonizData").
          get<Teuchos::RCP<std::map<double,double>>>("accIonizEn");
        Ea = evaluateIonizEnFromFile(*conc1,*ioniz_en1,Na*dens_sc);
      } else {
        Ea = params1.get<double>("Ionization Energy");
      }
      double gA = params1.get<double>("Degeneracy Factor");
      double p1 = Nv*exp(-Ea/kbT);
      double a3 = gD*p1;
      double a2 = gA*gD*gamma_n*ni*ni + gamma_n*n1*p1 + gD*p1*(Na+Nion);
      double a1 = gA*gamma_n*gamma_n*n1*ni*ni +
                  gamma_n*n1*p1*(Na-Nd+Nion) +
                  gA*gD*gamma_n*ni*ni*Nion;
      double a0 = -gA*gamma_n*gamma_n*n1*ni*ni*Nd +
                  gA*gamma_n*gamma_n*n1*ni*ni*Nion;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol1;
      sol1.push_back(r1_re); sol1.push_back(r2_re); sol1.push_back(r3_re);
      phy_sol = select_physical_sol(sol1, Nd);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
           "Error! Incomplete Ionization Donor Model MUST specify an \
             Approximation of type I, II or III");
    }
  } else { // hole density in p-type semic
    assert(incmpl_ioniz_param.sublist("Acceptor").numParams() > 0);
    const Teuchos::ParameterList& params =
      incmpl_ioniz_param.sublist("Acceptor");
    std::string approx_type = params.get<std::string>("Approximation");
    double Ea = 0.;
    if(params.isSublist("AccIncmplIonizData")) {
      Teuchos::RCP<std::vector<double> > conc =
        params.sublist("AccIncmplIonizData").
        get<Teuchos::RCP<std::vector<double> > >("accConc");
      Teuchos::RCP<std::map<double,double> > ioniz_en =
        params.sublist("AccIncmplIonizData").
        get<Teuchos::RCP<std::map<double,double> > >("accIonizEn");
      Ea = evaluateIonizEnFromFile(*conc,*ioniz_en,Na*dens_sc);
    } else {
      Ea = params.get<double>("Ionization Energy");
    }
    double gA = params.get<double>("Degeneracy Factor");
    double p1 = Nv*exp(-Ea/kbT);

    if(approx_type == "I") { // donor fully ionized
      double tt1 = gA*(Nd+Nion) + gamma_p*p1;
      const double tt2 = tt1*tt1 + 4*gA*gamma_p*p1*(Na-Nd-Nion);
      assert(tt2 >= 0);
      phy_sol = -tt1 + sqrt(tt2);
      assert(phy_sol >= 0);
      phy_sol /= 2*gA;
    } else if(approx_type == "II") { // donor neglected
      double tt1 = gA*Nion + gamma_p*p1;
      double tt2 = tt1*tt1 + 4*gA*gamma_p*p1*(Na-Nion);
      assert(tt2 >= 0);
      phy_sol = -tt1 + sqrt(tt2);
      assert(phy_sol >= 0);
      phy_sol /= 2*gA;
    } else if(approx_type == "III") { // full model
      // make sure that Donor Incomplete Ionization model is active
      if(incmpl_ioniz_param.sublist("Donor").numParams() == 0)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
             "Error! Incomplete Ionization Model for Donor MUST be defined \
             when Approximation type III is used for Acceptor Incomplete Ionization");
      const Teuchos::ParameterList& params1 =
        incmpl_ioniz_param.sublist("Donor");
      double Ed = 0.;
      if(params1.isSublist("DonIncmplIonizData")) {
        Teuchos::RCP<std::vector<double>> conc1 =
          params1.sublist("DonIncmplIonizData").
          get<Teuchos::RCP<std::vector<double> > >("donConc");
        Teuchos::RCP<std::map<double,double> > ioniz_en1 =
          params1.sublist("DonIncmplIonizData").
          get<Teuchos::RCP<std::map<double,double> > >("donIonizEn");
        Ed = evaluateIonizEnFromFile(*conc1,*ioniz_en1,Nd*dens_sc);
      } else {
        Ed = params.get<double>("Ionization Energy");
      }
      double gD = params1.get<double>("Degeneracy Factor");
      double n1 = Nc*exp(-Ed/kbT);
      double a3 = gA*n1;
      double a2 = gA*gD*gamma_p*ni*ni + gamma_p*n1*p1 + gA*Nd*n1;
      double a1 = gD*p1*gamma_p*gamma_p*ni*ni -
        gamma_p*n1*p1*(Na - Nd);
      double a0 = -gD*gamma_p*gamma_p*p1*ni*ni*Na;
      // solve for solution analytically
      double r1_re,r1_im,r2_re,r2_im,r3_re,r3_im;
      charon::cubicsolve(a3,a2,a1,a0,r1_re,r1_im,r2_re,r2_im,r3_re,r3_im);
      std::vector<double> sol1;
      sol1.push_back(r1_re); sol1.push_back(r2_re); sol1.push_back(r3_re);
      phy_sol = select_physical_sol(sol1, Na);
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
           "Error! Incomplete Ionization Acceptor Model MUST specify an \
            Approximation of type I, II or III");
    }
  }

  return phy_sol;
}



///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLatticeBC_OhmicContact<EvalT, Traits>::
DDLatticeBC_OhmicContact(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names")) ;
  std::string prefix = p.get<std::string>("Prefix");

  // basis
  RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary =
    p.get<RCP<const panzer::FieldLibraryBase> >("Field Library");

  RCP<const panzer::PureBasis> basis = fieldLayoutLibrary->lookupBasis(names.dof.phi);
  RCP<PHX::DataLayout> data_layout = basis->functional;
  num_basis = data_layout->dimension(1);

  // read in user-given values
  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(0);
  if (p.isType<double>("Voltage"))
    user_value->setRealValue(p.get<double>("Voltage"));
  else if (p.isType<std::string>("Varying Voltage"))
  {
    if (p.get<std::string>("Varying Voltage") == "Parameter")
      user_value =
        panzer::createAndRegisterScalarParameter<EvalT>(
        std::string("Varying Voltage"),
        *p.get<RCP<panzer::ParamLib> >("ParamLib"));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "DDLatticeBC_OhmicContact():  Error:  Expecting Varying Voltage "     \
        "value of \"Parameter\"; received \""
        << p.get<std::string>("Varying Voltage") << "\".")
  }
  bUseFD = p.get<bool>("Fermi Dirac");

  incmpl_ioniz = p.sublist("Incomplete Ionization");
  expandIonizEnParams(incmpl_ioniz);

  bSolveIon = p.get<bool>("Solve Ion");
  ion_charge = p.get<int>("Ion Charge");
  bUseFermiPin = p.get<bool>("Fermi Level Pinning");

  // evaluated fields
  potential = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.phi, data_layout);
  edensity = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.edensity, data_layout);
  hdensity = MDField<ScalarT,Cell,BASIS>(prefix+names.dof.hdensity, data_layout);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // add evaluated fields
  this->addEvaluatedField(potential);
  this->addEvaluatedField(edensity);
  this->addEvaluatedField(hdensity);

  // dependent fields
  doping = MDField<const ScalarT,Cell,BASIS>(names.field.doping_raw, data_layout);
  acceptor = MDField<const ScalarT,Cell,BASIS>(names.field.acceptor_raw, data_layout);
  donor = MDField<const ScalarT,Cell,BASIS>(names.field.donor_raw, data_layout);
  gamma_e = MDField<const ScalarT,Cell,BASIS>(names.field.elec_deg_factor, data_layout);
  gamma_h = MDField<const ScalarT,Cell,BASIS>(names.field.hole_deg_factor, data_layout);
  intrin_conc = MDField<const ScalarT,Cell,BASIS>(names.field.intrin_conc, data_layout);
  elec_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.elec_eff_dos, data_layout);
  hole_effdos = MDField<const ScalarT,Cell,BASIS>(names.field.hole_eff_dos, data_layout);
  eff_affinity =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_affinity, data_layout);
  eff_bandgap =  MDField<const ScalarT,Cell,BASIS>(names.field.eff_band_gap, data_layout);
  latt_temp =  MDField<const ScalarT,Cell,BASIS>(names.field.latt_temp, data_layout);

  // add dependent fields
  this->addDependentField(doping);
  this->addDependentField(acceptor);
  this->addDependentField(donor);
  this->addDependentField(gamma_e);
  this->addDependentField(gamma_h);
  this->addDependentField(intrin_conc);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(latt_temp);

  if (bSolveIon)
  {
    ion_density = MDField<const ScalarT,Cell,BASIS>(names.dof.iondensity, data_layout);
    this->addDependentField(ion_density);
  }

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral =
    Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

  std::string n = "Ohmic Contact for DDLattice";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLatticeBC_OhmicContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  // typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  // size_type num_basis = potential.dimension(1);
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& phyConst = charon::PhysicalConstants::Instance();
  double kbBoltz = phyConst.kb;   // Boltzmann constant in [eV/K]

  const bool withAccIncmplIoniz =
    (incmpl_ioniz.sublist("Acceptor").numParams() == 0) ? false : true;
  const bool withDonIncmplIoniz =
    (incmpl_ioniz.sublist("Donor").numParams() == 0) ? false : true;

 // use the Fermi-Dirac (FD) statistics. To compare FD with MB, we should choose
 // the right intrinsic conc. model so that nie = sqrt(Nc*Nv)*exp(-Egeff/2kbT) for MB
 if (bUseFD)
 {

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      // obtain concentrations (scaled)
      const ScalarT& dop = doping(cell,basis);
      const ScalarT& Nc = elec_effdos(cell,basis);
      const ScalarT& Nv = hole_effdos(cell,basis);
      const ScalarT& Na = acceptor(cell,basis);
      const ScalarT& Nd = donor(cell,basis);
      const ScalarT& nie = intrin_conc(cell,basis);
      const ScalarT& e_gamma = gamma_e(cell,basis);
      const ScalarT& h_gamma = gamma_h(cell,basis);

      // obtain energies in [eV]
      const ScalarT& effChi = eff_affinity(cell,basis);
      const ScalarT& effEg = eff_bandgap(cell,basis);

      // obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,basis)*T0;
      //if (lattT <= 0)
      //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
      //      "Error: Lattice temperature must be > 0 ! But its value is " << lattT << " ! \n");

       // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      ScalarT kbT = kbBoltz * lattT;  // [eV]
      ScalarT qtheta = effChi + 0.5*effEg + 0.5*kbT*log(Nc/Nv);  // [eV]

      // set Neff
      ScalarT Neff = dop;
      ScalarT Nion = 0.;
      if (bSolveIon)
      {
        //ScalarT Nion;
        if (bUseFermiPin)
          Nion = 5.e20 / C0;  // hard-code 5.e20 for now ???
        else
          Nion = ion_density(cell,basis);
        Neff = dop + ion_charge * Nion;
      }

      // compute carrier density and potential (scaled) at ohmic contact
      if (Neff >= 0.0)  // n-type contact
      {
        ScalarT n0 = 0.;
        if(!withDonIncmplIoniz) {
          n0 = Neff;  // assume n0 = Neff
        } else {
          double ND_crit = incmpl_ioniz.sublist("Donor").
            get<double>("Critical Doping Value");
          if(Nd > ND_crit/C0)
            // dopant fully ionized
            n0 = Neff;
          else
            // with incomplete ionized donor
            n0 = compute_FD_carrier_dens(-1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(e_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(h_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        ScalarT Ef_minus_Ec = kbT * (*inverseFermiIntegral)(n0/Nc);  // [eV]
        ScalarT Ef_minus_Ev = Ef_minus_Ec + effEg;  // [eV]
        ScalarT p0 = Nv * exp(-Ef_minus_Ev/kbT);
        ScalarT potBC = qtheta - effChi + Ef_minus_Ec +
          user_value->getValue();  // [eV]=[V]

        edensity(cell,basis) = n0;
        hdensity(cell,basis) = p0;
        potential(cell,basis) = potBC / V0;
      }
      else  // p-type contact
      {
        ScalarT p0 = 0.;
        if(!withAccIncmplIoniz) {
          p0 = -Neff;  // assume p0 = -Neff
        } else {
          double NA_crit = incmpl_ioniz.sublist("Acceptor").
            get<double>("Critical Doping Value");
          if(Na > NA_crit/C0)
            // dopant fully ionized
            p0 = -Neff;
          else
            // with incomplete ionized acceptor
            p0 = compute_FD_carrier_dens(1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(e_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(h_gamma),
                       Sacado::ScalarValue<ScalarT>::eval(C0));
        }
        ScalarT Ev_minus_Ef = kbT * (*inverseFermiIntegral)(p0/Nv);  // [eV]
        ScalarT Ec_minus_Ef = Ev_minus_Ef + effEg;  // [eV]
        ScalarT n0 = Nc * exp(-Ec_minus_Ef/kbT);
        ScalarT potBC = qtheta - effChi -effEg - Ev_minus_Ef +
          user_value->getValue();  // [eV]=[V]

        edensity(cell,basis) = n0;
        hdensity(cell,basis) = p0;
        potential(cell,basis) = potBC / V0;
      }

    }
  }
 }  // end of the if (bUseFD) block


 // use the Maxwell-Boltzmann (MB) statistics. The MB implementation assumes
 // nie = sqrt(Nc*Nv)*exp(-Egeff/2kbT).
 else
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      // obtain concentrations (scaled)
      const ScalarT& dop = doping(cell,basis);
      const ScalarT& Na = acceptor(cell,basis);
      const ScalarT& Nd = donor(cell,basis);
      const ScalarT& nie = intrin_conc(cell,basis);
      const ScalarT& Nc = elec_effdos(cell,basis);
      const ScalarT& Nv = hole_effdos(cell,basis);

      // obtain energies in [eV]
      const ScalarT& effChi = eff_affinity(cell,basis);
      const ScalarT& effEg = eff_bandgap(cell,basis);

      // obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,basis)*T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      ScalarT kbT = kbBoltz * lattT;  // [eV]
      ScalarT qtheta = effChi + 0.5*effEg + 0.5*kbT*log(Nc/Nv);  // [eV]

      // set Neff
      ScalarT Nion = 0.;
      ScalarT Neff = dop;
      if (bSolveIon)
      {
        if (bUseFermiPin)
          Nion = 5.e20 / C0;  // hard-code 5.e20 for now ???
        else
          Nion = ion_density(cell,basis);
        Neff = dop + ion_charge * Nion;
      }

      ScalarT y, potBC;

      // compute carrier density and potential (scaled) at ohmic contact
      if (Neff >= 0.0)  // n-type contact
      {
        ScalarT n0 = 0.;
        if(!withDonIncmplIoniz) {
          n0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) + 0.5*Neff;
        } else {
          double ND_crit = incmpl_ioniz.sublist("Donor").
            get<double>("Critical Doping Value");
          if(Nd > ND_crit/C0)
            // dopant fully ionized
            n0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) + 0.5*Neff;
          else
            // donor incomplete ionization active
            // acceptor incomplete ionization can be active or not
            // depending on the user input
            n0 = compute_MB_carrier_dens(-1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(C0) );
        }
        edensity(cell,basis) = n0;
        hdensity(cell,basis) = nie*nie/n0;

        y = 0.5*Neff/Nc + sqrt(pow(0.5*Neff/Nc,2.0) + Nv/Nc*exp(-effEg/kbT));
        potBC = qtheta - effChi + kbT*log(y) + user_value->getValue();  // [V]
        potential(cell,basis) = potBC / V0;  // scaled
      }

      else  // p-type contact
      {
        ScalarT p0 = 0.;
        if(!withAccIncmplIoniz) {
          p0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) - 0.5*Neff;
        } else {
          double NA_crit = incmpl_ioniz.sublist("Acceptor").
            get<double>("Critical Doping Value");
          if(Na > NA_crit/C0)
            p0 = sqrt(pow(0.5*Neff,2.0) + pow(nie,2.0)) - 0.5*Neff;
          else
            // acceptor incomplete ionization active
            // donor incomplete ionization can be active or not
            // depending on the user input
            p0 = compute_MB_carrier_dens(1, incmpl_ioniz,
                       Sacado::ScalarValue<ScalarT>::eval(kbT),
                       Sacado::ScalarValue<ScalarT>::eval(Nc),
                       Sacado::ScalarValue<ScalarT>::eval(Nv),
                       Sacado::ScalarValue<ScalarT>::eval(nie),
                       Sacado::ScalarValue<ScalarT>::eval(Na),
                       Sacado::ScalarValue<ScalarT>::eval(Nd),
                       Sacado::ScalarValue<ScalarT>::eval(ion_charge*Nion),
                       Sacado::ScalarValue<ScalarT>::eval(C0) );
        }
        hdensity(cell,basis) = p0;
        edensity(cell,basis) = pow(nie,2.0)/p0;

        y = -0.5*Neff/Nv + sqrt(pow(0.5*Neff/Nv,2.0) + Nc/Nv*exp(-effEg/kbT));
        potBC = qtheta - effChi - effEg - kbT*log(y) +
          user_value->getValue();  // [V]
        potential(cell,basis) = potBC / V0;  // scaled
      }

    }
  }
 }  // end of the else block

}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DDLatticeBC_OhmicContact<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Prefix", "");

  Teuchos::RCP<const panzer::FieldLibraryBase> fieldLayoutLibrary;
  p->set("Field Library", fieldLayoutLibrary);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);
  p->set<bool>("Frequency Domain", false);

  p->set<double>("Voltage", 0.0);
  p->set<std::string>("Varying Voltage", "Parameter");
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
    Teuchos::rcp(new panzer::ParamLib));
  p->set<bool>("Fermi Dirac", false);
  p->set<bool>("Acceptor Incomplete Ionization", false);
  p->set<bool>("Donor Incomplete Ionization", false);
  p->set<bool>("Solve Ion", false);
  p->set<int>("Ion Charge", 1);
  p->set<bool>("Fermi Level Pinning", false);
  p->sublist("Incomplete Ionization");
  p->sublist("Incomplete Ionization").sublist("Acceptor");
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Critical Doping Value", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Degeneracy Factor", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<double>("Ionization Energy", 0.0);
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<std::string>("AccIncmplIoniz File", "");
  p->sublist("Incomplete Ionization").sublist("Acceptor").
    set<std::string>("Approximation", "None");
  p->sublist("Incomplete Ionization").sublist("Donor");
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Critical Doping Value", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Degeneracy Factor", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<double>("Ionization Energy", 0.0);
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<std::string>("DonIncmplIoniz File", "");
  p->sublist("Incomplete Ionization").sublist("Donor").
    set<std::string>("Approximation", "None");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  Teuchos::RCP<charon::EmpiricalDamage_Data> dmgdata;
  p->set("empirical damage data", dmgdata);

  return p;
}

}

#endif

