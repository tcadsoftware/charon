
#ifndef CHARON_MOBILITY_SHIRAHATA_DECL_HPP
#define CHARON_MOBILITY_SHIRAHATA_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! Shirahata mobility model
template<typename EvalT, typename Traits>
class Mobility_Shirahata
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_Shirahata(
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

  // output
  PHX::MDField<ScalarT,Cell,Point> mobility;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp; //lattice temperature [scaled]

  // effective electric field at IP, scaled
  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field; 
  PHX::MDField<const ScalarT,Cell,Point> intrin_fermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,Point> bandgap;      // band gap w/o BGN in [eV]
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;  // effective band gap w/ BGN in [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
  double T0;  // temperature scaling, [K]
  double E0;  // electric field scaling, [K]
  double X0;  // distance scaling, [cm]

  int num_points;
  int num_dims;
  int num_edges;
  bool isEdgedl;
  bool isNodalDL;

  std::string carrType;
  //For nodes
  int basis_index;
  std::string basis_name;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  // Shirahata mobility model parameters
  double muo, E1, E2;
  double P1, P2, theta;
  double criticalDistance;
  double xOIStart,yOIStart,zOIStart,xOIEnd,yOIEnd,zOIEnd;
  double potentialSign;
  std::vector<double> oxideNorm;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_Shirahata


}

#endif
