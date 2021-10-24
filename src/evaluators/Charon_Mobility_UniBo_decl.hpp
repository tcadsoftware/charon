
#ifndef CHARON_MOBILITY_UNIBO_DECL_HPP
#define CHARON_MOBILITY_UNIBO_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! University of Bologna mobility model
template<typename EvalT, typename Traits>
class Mobility_UniBo
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_UniBo(
      const Teuchos::ParameterList& p);

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

  PHX::MDField<const ScalarT,Cell,Point> acceptor;
  PHX::MDField<const ScalarT,Cell,Point> donor;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
  double T0;  // temperature scaling, [K]

  int num_points;
  int num_edges;
  bool isEdgedl;

  std::string carrType;

  // UniBo mobility model parameters
  double mumax, mu0d, mu0a, mu1d, mu1a;
  double mu0d_exp, mu0a_exp, mu1d_exp, mu1a_exp;
  double cr1, cr2, cs1, cs2;
  double cr1_exp, cr2_exp, cs1_exp, cs2_exp;
  double c, gamma, alpha1, alpha2;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_UniBo


}

#endif
