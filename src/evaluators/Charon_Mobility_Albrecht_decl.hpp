
#ifndef CHARON_MOBILITY_ALBRECHT_DECL_HPP
#define CHARON_MOBILITY_ALBRECHT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! Albrecht mobility model

template<typename EvalT, typename Traits>
class Mobility_Albrecht
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_Albrecht(
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

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
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

  // low field mobility in [cm^2/(V.s)]
  double lowMob;

  // Albrecht mobility model parameters
  double exa, exb, exc;
  double exN0, exT0, exT1;

  // want mobility at the edge data layout ?
  bool isEdgedl;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_Albrecht


}

#endif
