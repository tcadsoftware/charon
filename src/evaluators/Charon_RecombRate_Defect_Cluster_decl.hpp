
#ifndef CHARON_RECOMBRATE_DEFECT_CLUSTER_DECL_HPP
#define CHARON_RECOMBRATE_DEFECT_CLUSTER_DECL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_interp.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_Vector.hpp"

// Phalanx
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

// Panzer
#include "Panzer_Dimension.hpp"

namespace charon {

using panzer::Cell;
using panzer::Point;

//! obtain SRH recombination rate
template<typename EvalT, typename Traits>
class RecombRate_Defect_Cluster
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RecombRate_Defect_Cluster(
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
  Teuchos::RCP<charon::clusterInterpolator> interpolator;

  // output
  PHX::MDField<ScalarT,Cell,Point> defect_cluster_rate;

  // input
  PHX::MDField<const ScalarT,Cell,Point> intrin_conc;
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0;
  double C0;

  //This is part of the cluster KLUDGE LCM
  PHX::MDField<const ScalarT,Cell,Point> acceptor;
  PHX::MDField<const ScalarT,Cell,Point> donor;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  int num_points;
  bool isIPset;

  // for interpolator
  std::string methodName;
  std::string inputType;
  double cascadeDensity;
  int numberOfFiles;
  float shepardPneg;
  double influenceRadius;

  //For in situ clustering

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class RecombRate_Defect_Cluster


}


#endif
