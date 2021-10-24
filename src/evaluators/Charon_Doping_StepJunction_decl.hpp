
#ifndef CHARON_DOPING_STEPJUNCTION_DECL_HPP
#define CHARON_DOPING_STEPJUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

namespace charon {

//! obtain net, acceptor, and donor doping for step junction
template<typename EvalT, typename Traits>
class Doping_StepJunction
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Doping_StepJunction(
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

public:
  std::vector<double> evaluateDoping(const double& x, const double& y, const double& z);

private:
  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]

  // output
  PHX::MDField<ScalarT,Cell,IP> doping_raw; // net doping @ IPs
  PHX::MDField<ScalarT,Cell,IP> acceptor_raw;
  PHX::MDField<ScalarT,Cell,IP> donor_raw;

  PHX::MDField<ScalarT,Cell,BASIS> doping_raw_basis; // @ basis points
  PHX::MDField<ScalarT,Cell,BASIS> acceptor_raw_basis;
  PHX::MDField<ScalarT,Cell,BASIS> donor_raw_basis;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  // input from closure model factory
  double acceptorValue; // Acceptor doping value [cm-3]
  double donorValue;    // Donor doping value [cm-3]
  double junctionLoc;   // Junction location [um]

  std::string config;    // Junction configuration: PN or NP
  std::string direction; // Junction direction: X, or Y, or Z

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Doping_StepJunction


}

#endif
