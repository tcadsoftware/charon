
#ifndef CHARON_DOPING_FUNCTION_DECL_HPP
#define CHARON_DOPING_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Doping_Params.hpp"
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Charon_MMS_AnalyticFunctions.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

namespace charon {

//! obtain a uniform doping
template<typename EvalT, typename Traits>
class Doping_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Doping_Function(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]
  double T0; // temp scaling, [k]

  // input
  PHX::MDField<const ScalarT,Cell,IP> doping_raw;
  PHX::MDField<const ScalarT,Cell,IP> acceptor_raw;
  PHX::MDField<const ScalarT,Cell,IP> donor_raw;
  PHX::MDField<const ScalarT,Cell,IP> latt_temp;
  PHX::MDField<const ScalarT,Cell,IP> dens_e;
  PHX::MDField<const ScalarT,Cell,IP> dens_h;
  PHX::MDField<const ScalarT,Cell,IP> gamma_e;
  PHX::MDField<const ScalarT,Cell,IP> gamma_h;
  PHX::MDField<const ScalarT,Cell,IP> elec_effdos;
  PHX::MDField<const ScalarT,Cell,IP> hole_effdos;

  PHX::MDField<const ScalarT,Cell,BASIS> doping_raw_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> acceptor_raw_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> donor_raw_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> dens_e_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> dens_h_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> gamma_e_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> gamma_h_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos_basis;
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos_basis;

  // output
  PHX::MDField<ScalarT,Cell,IP> doping; // net doping @ IPs
  PHX::MDField<ScalarT,Cell,IP> acceptor;
  PHX::MDField<ScalarT,Cell,IP> donor;

  PHX::MDField<ScalarT,Cell,BASIS> doping_basis; // @ BASIS points
  PHX::MDField<ScalarT,Cell,BASIS> acceptor_basis;
  PHX::MDField<ScalarT,Cell,BASIS> donor_basis;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  //Teuchos::ParameterList dopParamList;
  bool with_IonizAcc;
  bool with_IonizDon;
  Teuchos::ParameterList ionizacc_dopParamList;
  Teuchos::ParameterList ionizdon_dopParamList;

  double gD_acc, en_ionz_acc, NA_crit;
  double gD_don, en_ionz_don, ND_crit;
  bool WithAccEnFromFile, WithDonEnFromFile;
  std::vector<double> accConc;
  std::vector<double> donConc;
  std::map<double,double> accIonizEn;
  std::map<double,double> donIonizEn;

  void initParam(const Teuchos::ParameterList& ionizacc_dopParamList,
                 const Teuchos::ParameterList& ionizdon_dopParamList);

  double evaluateIonizEnFromFile(const std::vector<double>& dop_table,
                                 const std::map<double,double>& IonizEn,
                                 double dop);

  std::vector<uniformDopingParams> udp_vec;
  std::vector<gaussianDopingParams> gdp_vec;
  std::vector<linearDopingParams> ldp_vec;
  std::vector<erfcDopingParams> edp_vec;
  std::vector<mgaussDopingParams> mgdp_vec;

}; // end of class Doping_Function


}

#endif
