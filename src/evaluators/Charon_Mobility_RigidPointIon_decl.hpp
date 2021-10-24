
#ifndef CHARON_MOBILITY_RIGIDPOINTION_DECL_HPP
#define CHARON_MOBILITY_RIGIDPOINTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! RigidPointIon mobility model
template<typename EvalT, typename Traits>
class Mobility_RigidPointIon
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_RigidPointIon(
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

  // compute the rigid point ion mobility
  ScalarT computeIonMobility(const ScalarT& kbT, const ScalarT& ionDens);

  // compute the rigid point ion velocity (one component at a time)
  ScalarT computeIonVelocity(const ScalarT& kbT, const ScalarT& ionMob, const ScalarT& ionF);

  // output
  PHX::MDField<ScalarT,Cell,Point> mobility;
  PHX::MDField<ScalarT,Cell,Point> edge_velocity;
  PHX::MDField<ScalarT,Cell,Point,Dim> ip_velocity;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // lattice temperature [scaled]
  PHX::MDField<const ScalarT,Cell,Point> carr_dens; // ion density [scaled]
  PHX::MDField<const ScalarT,Cell,Point> potential; // electric potential [scaled]

  PHX::MDField<const ScalarT,Cell,Point,Dim> ip_ionfield;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double E0;  // electric field scaling, [V/cm]
  double T0;  // temperature scaling, [K]
  double C0;  // conc. scaling, [cm^-3]

  int num_points;
  int num_dims;
  int num_edges;

  // RigidPointIon mobility model parameters
  double escFreq, hopDist, actE, sign;
  double maxIonDens;
  double velMultiplier;

  bool bSetMaxDens;
  bool isEdgedl;

  std::string basis_name;
  std::size_t basis_index;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

}; // end of class Mobility_RigidPointIon


}

#endif
