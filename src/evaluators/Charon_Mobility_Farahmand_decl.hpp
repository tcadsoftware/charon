
#ifndef CHARON_MOBILITY_FARAHMAND_DECL_HPP
#define CHARON_MOBILITY_FARAHMAND_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! Farahmand mobility model

template<typename EvalT, typename Traits>
class Mobility_Farahmand
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_Farahmand(
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

  // evaluate the low field mobility
  ScalarT evaluateLowFieldMobility(const ScalarT& Na, const ScalarT& Nd, const ScalarT& latt);

  // compute GaAs mobility for CVFEM-SG and EFFPG-FEM at primary edge center
  ScalarT evaluateMobilityForEdgedl(const std::size_t& cell, const int& edge, const ScalarT& edgelfMob,
    const Kokkos::DynRankView<double,PHX::Device>& edgePoints, const ScalarT& edgeLatt);

  // compute GaAs mobility for SUPG-FEM at IP
  ScalarT evaluateMobilityForIPdl(const std::size_t& cell, const int& point, const ScalarT& lfMob);

  // evaluate high field mobility
  ScalarT evaluateHighFieldMobility(const ScalarT& lfMob, const ScalarT& hiField);

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

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  // carrier type
  std::string carrType;

  // driving force name for high field
  std::string driveForce;

  double sign;   // + for n and - for p

  // Farahmand mobility model parameters
  double mu1, mu2, alpha, beta, delta, gamma, eps, ncrit;

  // High field mobility parameters
  double vsat, ec, n1, n2, an;

  // turn on/off high field dependence
  bool hiFieldOn;

  // want mobility at the edge data layout ?
  bool isEdgedl;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_Farahmand


}

#endif
