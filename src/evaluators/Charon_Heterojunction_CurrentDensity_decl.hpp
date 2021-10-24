
#ifndef CHARON_HETEROJUNCTION_CURRENTDENSITY_DECL_HPP
#define CHARON_HETEROJUNCTION_CURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Panzer_Evaluator_Macros.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::IP;
using panzer::BASIS;
using panzer::Dim;

namespace charon {

PANZER_EVALUATOR_CLASS(Heterojunction_CurrentDensity)

private:

  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // scaled
  PHX::MDField<const ScalarT,Cell,Point> carr_dens;
  PHX::MDField<const ScalarT,Cell,Point> eff_dos;
  PHX::MDField<const ScalarT,Cell,Point> other_carr_dens;
  PHX::MDField<ScalarT,Cell,IP> hj_flux;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling [K]
  double C0; // concentration scaling [cm^-3]
  double J0; // current density scaling [A/cm^2]

  int num_points;
  int num_ips;
  int num_nodes;
  int detailIndex;

  // for basis
  std::string basis_name;
  std::size_t basis_index;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  double richConst;   // [A/(cm^2.K^2)]
  double bandOffset;  // [eV]
  double fdDensity;   // [cm^-3]
  double multiplier;

  std::string dof_name;
  std::string flux_name;
  std::string other_dof_name;
  std::string discMethod;

  bool isPrimarySideLeft;

  const std::string femsupg = "FEM-SUPG";
  const std::string cvfemsg = "CVFEM-SG";

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

PANZER_EVALUATOR_CLASS_END

}

#endif
