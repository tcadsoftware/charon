
#ifndef CHARON_DDLATTICE_HEATGENERATION_DECL_HPP
#define CHARON_DDLATTICE_HEATGENERATION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! obtain the heat generation
template<typename EvalT, typename Traits>
class DDLattice_HeatGeneration
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLattice_HeatGeneration(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> heat_gen;  // [scaled]

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_curr_dens; // electron current density [scaled]
  PHX::MDField<const ScalarT,Cell,Point,Dim> hole_curr_dens; // hole current density [scaled]
  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_field; // electron electric field [scaled]
  PHX::MDField<const ScalarT,Cell,Point,Dim> hole_field; // hole electric field [scaled]

  PHX::MDField<const ScalarT,Cell,Point> latt_temp;    // lattice temperature [scaled]
  PHX::MDField<const ScalarT,Cell,Point> total_recomb; // total net recombination [scaled]

  PHX::MDField<const ScalarT,Cell,Point> eff_band_gap; // effective band gap [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // Temperature Scaling [K]
  double J0; // current density scaling [A/cm^2]
  double E0; // electric field scaling [V/cm]
  double H0; // heat generation scaling [W/cm^3]

  PHX::MDField<const ScalarT,Cell,Point,Dim> ion_curr_dens; // ion current density [scaled]
  PHX::MDField<const ScalarT,Cell,Point,Dim> ion_field;     // ion electric field [scaled]

  int num_points;
  int num_dims;

  bool solveIon;
  bool solveHole;
  bool solveElectron;
  bool haveSource;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DDLattice_HeatGeneration


}

#endif
