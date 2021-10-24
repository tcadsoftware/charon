
#ifndef CHARON_IONIZATION_PARTICLE_STRIKE_DECL_HPP
#define CHARON_IONIZATION_PARTICLE_STRIKE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_interp.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain ionization rate due to particle strike
template<typename EvalT, typename Traits>
class Ionization_ParticleStrike
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Ionization_ParticleStrike(
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

  // output
  PHX::MDField<ScalarT,Cell,Point> ionization_particle_strike_rate;

  // input
  PHX::MDField<const ScalarT,Cell,Point> intrin_conc;
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;
  double startX,startY,startZ,endX,endY,endZ;
  double strikeRadius, generationRate, totalCharge;
  double startTime,endTime;
  std::string temporalWaveform;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // [s]
  double C0; // [cm^-3]

  std::vector<int> amILitUp(std::vector<double> nodalX, std::vector<double> nodalY, std::vector<double> nodalZ);
  double distance(double x1, double y1, double z1, double x2, double y2, double z2);
  double PointInCylinderTest(
     double pt1_x, double pt1_y, double pt1_z, double pt2_x, double pt2_y, double pt2_z, 
     double lengthsq, double radius_sq, double testpt_x, double testpt_y, double testpt_z);

  double getTimeFactor(double currentTime);

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  int num_points, num_dim;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Ionization_ParticleStrike


}


#endif
