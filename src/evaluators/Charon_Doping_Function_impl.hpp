
#ifndef CHARON_DOPING_FUNCTION_IMPL_HPP
#define CHARON_DOPING_FUNCTION_IMPL_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Util.hpp"

namespace charon {


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Doping_Function<EvalT, Traits>::
Doping_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  basis_name = basis->name();

  with_IonizAcc = p.isSublist("IncmplIonizAcc Doping ParameterList");
  with_IonizDon = p.isSublist("IncmplIonizDon Doping ParameterList");
  if(with_IonizAcc)
    ionizacc_dopParamList = p.sublist("IncmplIonizAcc Doping ParameterList");
  if(with_IonizDon)
    ionizdon_dopParamList = p.sublist("IncmplIonizDon Doping ParameterList");

  // initialize model parameters
  initParam(ionizacc_dopParamList, ionizdon_dopParamList);

  // output fields
  doping = MDField<ScalarT,Cell,IP>(n.field.doping,scalar);
  acceptor = MDField<ScalarT,Cell,IP>(n.field.acceptor,scalar);
  donor = MDField<ScalarT,Cell,IP>(n.field.donor,scalar);

  doping_basis = MDField<ScalarT,Cell,BASIS>(n.field.doping,data_layout);
  acceptor_basis = MDField<ScalarT,Cell,BASIS>(n.field.acceptor,data_layout);
  donor_basis = MDField<ScalarT,Cell,BASIS>(n.field.donor,data_layout);

  this->addEvaluatedField(doping);
  this->addEvaluatedField(acceptor);
  this->addEvaluatedField(donor);

  this->addEvaluatedField(doping_basis);
  this->addEvaluatedField(acceptor_basis);
  this->addEvaluatedField(donor_basis);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;

  // input
  doping_raw = MDField<const ScalarT,Cell,IP>(n.field.doping_raw,scalar);
  acceptor_raw = MDField<const ScalarT,Cell,IP>(n.field.acceptor_raw,scalar);
  donor_raw = MDField<const ScalarT,Cell,IP>(n.field.donor_raw,scalar);
  doping_raw_basis = MDField<const ScalarT,Cell,BASIS>(n.field.doping_raw,data_layout);
  acceptor_raw_basis = MDField<const ScalarT,Cell,BASIS>(n.field.acceptor_raw,data_layout);
  donor_raw_basis = MDField<const ScalarT,Cell,BASIS>(n.field.donor_raw,data_layout);
  gamma_e = MDField<const ScalarT,Cell,IP>(n.field.elec_deg_factor,scalar);
  gamma_h = MDField<const ScalarT,Cell,IP>(n.field.hole_deg_factor,scalar);
  gamma_e_basis =
    MDField<const ScalarT,Cell,BASIS>(n.field.elec_deg_factor,data_layout);
  gamma_h_basis =
    MDField<const ScalarT,Cell,BASIS>(n.field.hole_deg_factor,data_layout);
  // dependencies for Incomplete Ionized Acceptor
  if(with_IonizAcc) {
    dens_h = MDField<const ScalarT,Cell,IP>(n.dof.hdensity,scalar);
    dens_h_basis =
      MDField<const ScalarT,Cell,BASIS>(n.dof.hdensity,data_layout);
    hole_effdos =
      MDField<const ScalarT,Cell,IP>(n.field.hole_eff_dos,scalar);
    hole_effdos_basis =
      MDField<const ScalarT,Cell,BASIS>(n.field.hole_eff_dos,data_layout);
  }
  // dependencies for Incomplete Ionized Donor
  if(with_IonizDon) {
    dens_e = MDField<const ScalarT,Cell,IP>(n.dof.edensity,scalar);
    dens_e_basis =
      MDField<const ScalarT,Cell,BASIS>(n.dof.edensity,data_layout);
    elec_effdos =
      MDField<const ScalarT,Cell,IP>(n.field.elec_eff_dos,scalar);
    elec_effdos_basis =
      MDField<const ScalarT,Cell,BASIS>(n.field.elec_eff_dos,data_layout);
  }
  if(with_IonizAcc || with_IonizDon) {
    T0 = scaleParams->scale_params.T0;
    latt_temp = MDField<const ScalarT,Cell,IP>(n.field.latt_temp,scalar);
    latt_temp_basis =
      MDField<const ScalarT,Cell,BASIS>(n.field.latt_temp,data_layout);
  }

  this->addDependentField(doping_raw);
  this->addDependentField(acceptor_raw);
  this->addDependentField(donor_raw);
  this->addDependentField(doping_raw_basis);
  this->addDependentField(acceptor_raw_basis);
  this->addDependentField(donor_raw_basis);
  this->addDependentField(gamma_e);
  this->addDependentField(gamma_h);
  this->addDependentField(gamma_e_basis);
  this->addDependentField(gamma_h_basis);
  if(with_IonizAcc) {
    this->addDependentField(dens_h);
    this->addDependentField(dens_h_basis);
    this->addDependentField(hole_effdos);
    this->addDependentField(hole_effdos_basis);
  }
  if(with_IonizDon) {
    this->addDependentField(dens_e);
    this->addDependentField(dens_e_basis);
    this->addDependentField(elec_effdos);
    this->addDependentField(elec_effdos_basis);
  }
  if(with_IonizAcc || with_IonizDon) {
    this->addDependentField(latt_temp);
    this->addDependentField(latt_temp_basis);
  }

  std::string name = "Doping_Function";
  this->setName(name);
}



///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Doping_Function<EvalT, Traits>::
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
Doping_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ using std::cout;
  using std::endl;
  using std::ofstream;
  using std::ios;
  using panzer::index_t;

  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  size_type num_basis = doping_basis.dimension(1);

  // Boltzmann constant in [eV/K]
  double kbBoltz = charon::PhysicalConstants::Instance().kb;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // doping at IPs
    for (int ip = 0; ip < num_ip; ++ip)
    {
      // evaluate the acceptor and donor doping
      acceptor(cell,ip) = acceptor_raw(cell,ip);
      donor(cell,ip) = donor_raw(cell,ip);
      if(with_IonizAcc) {
        if(acceptor(cell,ip) <= NA_crit/C0) {
          double en_ionz = WithAccEnFromFile
            ? evaluateIonizEnFromFile(accConc,accIonizEn,
                  Sacado::ScalarValue<ScalarT>::eval(acceptor(cell,ip)*C0))
            : en_ionz_acc;
          ScalarT kbT = kbBoltz * latt_temp(cell,ip) * T0;  // [eV]
          ScalarT p1 = hole_effdos(cell,ip) * exp(-en_ionz/kbT);
          // partially ionized acceptor conc
          acceptor(cell,ip) = acceptor(cell,ip) /
            ( 1. + gD_acc*dens_h(cell,ip)/gamma_h(cell,ip)/p1 );
        }
      }
      if(with_IonizDon) {
        if(donor(cell,ip) <= ND_crit/C0) {
          double en_ionz = WithDonEnFromFile
            ? evaluateIonizEnFromFile(donConc,donIonizEn,
                  Sacado::ScalarValue<ScalarT>::eval(donor(cell,ip)*C0))
            : en_ionz_don;
          ScalarT kbT = kbBoltz * latt_temp(cell,ip) * T0;  // [eV]
          ScalarT n1 = elec_effdos(cell,ip) * exp(-en_ionz/kbT);
          // partially ionized donor conc
          donor(cell,ip) = donor(cell,ip) /
            ( 1. + gD_don*dens_e(cell,ip)/gamma_e(cell,ip)/n1 );
        }
      }
      doping(cell,ip) = donor(cell,ip) - acceptor(cell,ip);
    }

    // doping at basis points
    for (size_type basis = 0; basis < num_basis; ++basis)
    {
      acceptor_basis(cell,basis) = acceptor_raw_basis(cell,basis);
      donor_basis(cell,basis) = donor_raw_basis(cell,basis);
      if(with_IonizAcc) {
        if(acceptor_basis(cell,basis) <= NA_crit/C0) {
          double en_ionz = WithAccEnFromFile
            ? evaluateIonizEnFromFile(accConc,accIonizEn,
                  Sacado::ScalarValue<ScalarT>::eval(acceptor(cell,basis)*C0))
            : en_ionz_acc;
          ScalarT kbT = kbBoltz * latt_temp_basis(cell,basis) * T0;  // [eV]
          ScalarT p1 = hole_effdos_basis(cell,basis) * exp(-en_ionz/kbT);
          // partially ionized acceptor conc
          acceptor_basis(cell,basis) = acceptor_basis(cell,basis) /
            ( 1. + gD_acc*dens_h_basis(cell,basis)/gamma_h_basis(cell,basis)/p1 );
        }
      }
      if(with_IonizDon) {
        if(donor_basis(cell,basis) <= ND_crit/C0) {
          double en_ionz = WithDonEnFromFile
            ? evaluateIonizEnFromFile(donConc,donIonizEn,
                  Sacado::ScalarValue<ScalarT>::eval(donor(cell,basis)*C0))
            : en_ionz_don;
          ScalarT kbT = kbBoltz * latt_temp_basis(cell,basis) * T0;  // [eV]
          ScalarT n1 = elec_effdos_basis(cell,basis) * exp(-en_ionz/kbT);
          // partially ionized donor conc
          donor_basis(cell,basis) = donor_basis(cell,basis) /
            ( 1. + gD_don*dens_e_basis(cell,basis)/gamma_e_basis(cell,basis)/n1 );
        }
      }
      doping_basis(cell,basis) = donor_basis(cell,basis) - acceptor_basis(cell,basis);
    }

  } // end of loop over cells
}


template<typename EvalT, typename Traits>
void Doping_Function<EvalT, Traits>::initParam(
          const Teuchos::ParameterList& ionizacc_dopParamList,
          const Teuchos::ParameterList& ionizdon_dopParamList) {
  gD_acc = 0.; en_ionz_acc = 0.; NA_crit = 0.;
  gD_don = 0.; en_ionz_don = 0.; ND_crit = 0.;
  WithAccEnFromFile = false; WithDonEnFromFile = false;

  if(with_IonizAcc) {
    gD_acc = ionizacc_dopParamList.get<double>("Degeneracy Factor");
    NA_crit = ionizacc_dopParamList.get<double>("Critical Doping Value");

    if(ionizacc_dopParamList.isParameter("Ionization Energy") &&
       ionizacc_dopParamList.isParameter("AccIncmplIoniz File") )
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
         << "Either a value or a file must be specified for Acceptor Ionization Energy, not both" << std::endl);

    if(ionizacc_dopParamList.isParameter("Ionization Energy")) {
      en_ionz_acc = ionizacc_dopParamList.get<double>("Ionization Energy");
    } else {
      // read ionization energy from file
      std::string ionizEnFile =
        ionizacc_dopParamList.get<std::string>("AccIncmplIoniz File");
      if(ionizEnFile == "")
       TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "AccIncmplIoniz File name cannot be empty !" << std::endl);
      std::ifstream info(ionizEnFile.c_str());
      if(!info)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "Cannot read AccIncmplIoniz File '" << ionizEnFile << "'" << std::endl);
      WithAccEnFromFile =  true;
      std::string line;
      while(std::getline(info,line)) {
        std::istringstream iss(line);
        double dop=0, ioniz_en=0;
        if(!(iss >> dop >> ioniz_en)) { break; }
        if(dop <= 0) dop = 1e-100;
        accConc.push_back(dop);
        accIonizEn.insert(std::pair<double,double>(dop,ioniz_en));
      }
      std::sort(accConc.begin(),accConc.end());
    }
  }

  if(with_IonizDon) {
    gD_don = ionizdon_dopParamList.get<double>("Degeneracy Factor");
    ND_crit = ionizdon_dopParamList.get<double>("Critical Doping Value");

    if(ionizdon_dopParamList.isParameter("Ionization Energy") &&
       ionizdon_dopParamList.isParameter("DonIncmplIoniz File") )
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
         << "Either a value or a file must be specified for Donnor Ionization Energy, not both" << std::endl);

    if(ionizdon_dopParamList.isParameter("Ionization Energy")) {
      en_ionz_don = ionizdon_dopParamList.get<double>("Ionization Energy");
    } else {
      // read ionization energy from file
      std::string ionizEnFile =
        ionizdon_dopParamList.get<std::string>("DonIncmplIoniz File");
      if(ionizEnFile == "")
       TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "DonIncmplIoniz File name cannot be empty !" << std::endl);
      std::ifstream info(ionizEnFile.c_str());
      if(!info)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "Cannot read DonIncmplIoniz File '" << ionizEnFile << "'" << std::endl);
      WithDonEnFromFile =  true;
      std::string line;
      while(std::getline(info,line)) {
        std::istringstream iss(line);
        double dop=0, ioniz_en=0;
        if(!(iss >> dop >> ioniz_en)) { break; }
        if(dop <= 0) dop = 1e-100;
        donConc.push_back(dop);
        donIonizEn.insert(std::pair<double,double>(dop,ioniz_en));
      }
      std::sort(donConc.begin(),donConc.end());
    }
  }
}


template<typename EvalT, typename Traits>
double Doping_Function<EvalT, Traits>::evaluateIonizEnFromFile(
                     const std::vector<double>& dop_table,
                     const std::map<double,double>& IonizEn,
                     double dop) {
  // dop_table is assumed to be ordered from low to high doping
  // find the two adjacent values in the dop_table bounding dop
  std::pair<double,double> dop_bounds;
  findDopingPoints(dop_table,dop,dop_bounds);
  // compute ionization energy
  return interpolateIonizEn(IonizEn,dop_bounds,dop);
}


} // namespace charon

#endif
