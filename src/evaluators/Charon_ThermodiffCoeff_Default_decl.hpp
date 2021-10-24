
#ifndef CHARON_THERMODIFFCOEFF_DEFAULT_DECL_HPP
#define CHARON_THERMODIFFCOEFF_DEFAULT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"


using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain thermodiffusion coefficient for ions
template<typename EvalT, typename Traits>
class ThermodiffCoeff_Default
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    ThermodiffCoeff_Default(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> thermodiff_coeff; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,Point> mobility;    // scaled
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;   // scaled
  PHX::MDField<const ScalarT,Cell,Point> soret_coeff; // scaled

  int num_points;
  int num_edges;

  bool isEdgedl;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class ThermodiffCoeff_Default


}

#endif
