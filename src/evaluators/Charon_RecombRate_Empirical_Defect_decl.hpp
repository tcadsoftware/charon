
#ifndef CHARON_RECOMBRATE_EMPIRICAL_DEFECT_DECL_HPP
#define CHARON_RECOMBRATE_EMPIRICAL_DEFECT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_empiricalConvolution.hpp"
#include "Charon_PulseDamage_Spec.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

#include "Charon_Scaling_Parameters.hpp"

namespace charon {

//! Calculates the recombination due to neutron damage using the empirical model
template <typename EvalT, typename Traits, typename PointType>
class RecombRate_Empirical_Defect :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{

public:
  RecombRate_Empirical_Defect(panzer::PureBasis const& basis,
                              panzer::IntegrationRule const& ir,
                              Teuchos::ParameterList const& p);

  void evaluateFields(typename Traits::EvalData d);

private:

  // used for member specialization of getCoordField() methods
  template <typename T> struct Type {};

  typedef typename EvalT::ScalarT ScalarT;

  //! Get the coordinates of the IP (panzer::Point) points
  PHX::MDField<const ScalarT,panzer::Cell,PointType,panzer::Dim> getCoordField(Type<panzer::Point>,
                                                                               panzer::PureBasis const& basis,
                                                                               panzer::IntegrationRule const& ir);

  //! Get the coordinates of the basis (panzer::Basis) points
  PHX::MDField<const ScalarT,panzer::Cell,PointType,panzer::Dim> getCoordField(Type<panzer::BASIS>,
                                                                               panzer::PureBasis const& basis,
                                                                               panzer::IntegrationRule const& ir);

  PHX::MDField<const ScalarT,panzer::Cell,PointType,panzer::Dim> coordinates;

  Teuchos::RCP<charon::empiricalConvolution> NfpMu;

  Teuchos::RCP<charon::PulseDamage_Spec> damage_spec;

  Teuchos::RCP<charon::EmpiricalDamage_Data> damage_data;

  // output
  PHX::MDField<ScalarT,panzer::Cell,PointType> empirical_defect_rate;

  // input
  PHX::MDField<const ScalarT,panzer::Cell,PointType> intrin_conc;
  PHX::MDField<const ScalarT,panzer::Cell,PointType> edensity;
  PHX::MDField<const ScalarT,panzer::Cell,PointType> hdensity;

  //scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0;
  double C0;
  double X0;

  int num_points;

  double ebxlo,ebxhi;
  double ebylo,ebyhi;
  double ebzlo,ebzhi;

  double thermalVelocity;
  double crossSection;


  std::string pulseType;
  std::string muDataFile;

  double pulseStart;
  double pulseEnd;
  double pulseMagnitude;

  bool ebOverrideBool;
  double ebVoltageOverride;

  bool cbOverrideBool;
  double cbVoltageOverride;

  int pulsePoints;

  std::string pulseDataFile;
  std::string fileDataPulses;

  bool pulseIsRate;

};

}

#endif
