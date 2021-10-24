
#ifndef CHARON_MOBILITY_MASETTI_DECL_HPP
#define CHARON_MOBILITY_MASETTI_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! Masetti mobility model
template<typename EvalT, typename Traits>
class Mobility_Masetti
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_Masetti(
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

  // Masetti mobility model parameters
  double mumax, mumin1, mumin2, mu1;
  double pc, cr, cs;
  double gamma, alpha, beta;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_Masetti


}

#endif
