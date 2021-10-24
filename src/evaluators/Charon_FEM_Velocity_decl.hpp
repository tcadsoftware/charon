
#ifndef CHARON_FEM_VELOCITY_DECL_HPP
#define CHARON_FEM_VELOCITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! obtain electron and hole velocity at IPs
template<typename EvalT, typename Traits>
class FEM_Velocity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    FEM_Velocity(
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
  PHX::MDField<ScalarT,Cell,Point,Dim> velocity; // @ IPs
  PHX::MDField<ScalarT,Cell,Point> edge_velocity; // @ Edges

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> efield;
  PHX::MDField<const ScalarT,Cell,Point> mobility;

  PHX::MDField<const ScalarT,Cell,Point> potential;

  // IP
  std::size_t num_points;
  std::size_t num_dims;
  int num_edges;

  // carrier type
  std::string carrType;

  // sign (- for electron and + for hole)
  double sign;

  bool isEdgedl;

  std::string basis_name;
  std::size_t basis_index;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

}; // end of class FEM_Velocity


}

#endif
