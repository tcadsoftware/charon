
#ifndef CHARON_MOBILITY_LUCENT_DECL_HPP
#define CHARON_MOBILITY_LUCENT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! Lucent mobility model that excludes the transverse-field dependence

template<typename EvalT, typename Traits>
class Mobility_Lucent
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_Lucent(
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

  // initialize mobility parameters
  void initMobilityParams(const std::string& matName, const Teuchos::ParameterList& mobParamList);

  // compute Philips mobility (low field)
  ScalarT evaluatePhilipsMobility(const ScalarT& Na, const ScalarT& Nd,
    const ScalarT& eden, const ScalarT& hden, const ScalarT& latt);

  // compute Lucent mobility for CVFEM-SG and EFFPG-FEM at primary edge center
  ScalarT evalLucentMobForEdgedl(const std::size_t& cell, const int& edge, const ScalarT& edgelfMob,
    const Kokkos::DynRankView<double,PHX::Device>& edgePoints, const ScalarT& edgeLatt);

  // compute Lucent mobility for SUPG-FEM at IP
  ScalarT evalLucentMobForIPdl(const std::size_t& cell, const int& point, const ScalarT& lfMob);


  // output
  PHX::MDField<ScalarT,Cell,Point> mobility; // @ IPs or Edge

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // scaled

  PHX::MDField<const ScalarT,Cell,Point> acceptor;  // scaled
  PHX::MDField<const ScalarT,Cell,Point> donor;     // scaled
  PHX::MDField<const ScalarT,Cell,Point> edensity;  // scaled
  PHX::MDField<const ScalarT,Cell,Point> hdensity;  // scaled

  PHX::MDField<const ScalarT,Cell,Point> intrin_fermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,Point> bandgap;      // band gap w/o BGN in [eV]
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;  // effective band gap w/ BGN in [eV]

  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp;  // gradient of quasi-fermi potential at IP, scaled
  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field; // effective electric field at IP, scaled

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
  double X0;  // length scaling, [cm]
  double E0;  // electric field scaling, [V/cm]
  double T0;  // temperature scaling, [K]

  int num_ips;
  int num_nodes;
  int num_points;
  int num_dims;
  int num_edges;

  std::string basis_name;
  std::size_t basis_index;

  // input from closure model factory
  std::string carrType;
  std::string gParaSet;

  // Philips mobility model parameters in physical units
  double mumax, mumin, nref;
  double gamma, alpha;
  double nref_d, nref_a, cref_d, cref_a;

  // G parameters in physical units
  double ag, bg, cg;
  double alpha_g, beta_g, alpha_prime_g, gamma_g;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  // driving force name for high field
  std::string driveForce;

  // saturation velocity in cm/s
  ScalarT vsat;

  // exponent in the high field model
  double beta_hf;

  // different sign for electrons and holes
  double sign;  // +1 for e and - for h

  // turn on/off high field dependence
  bool hiFieldOn;

  // want mobility at the edge data layout ?
  bool isEdgedl;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_Lucent


}

#endif
