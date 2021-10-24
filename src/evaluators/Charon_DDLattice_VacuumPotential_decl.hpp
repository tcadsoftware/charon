
#ifndef CHARON_DDLATTICE_VACUUMPOTENTIAL_DECL_HPP
#define CHARON_DDLATTICE_VACUUMPOTENTIAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class DDLattice_VacuumPotential
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLattice_VacuumPotential(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> vac_pot;  // vacuum potential

  // input
  PHX::MDField<const ScalarT,Cell,Point> potential;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> eff_affinity;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // Temperature Scaling in [K]

  std::size_t num_points;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DDLattice_VacuumPotential


}

#endif
