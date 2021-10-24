
#ifndef CHARON_AVALANCHE_CROWELLSZE_DECL_HPP
#define CHARON_AVALANCHE_CROWELLSZE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::BASIS;
using panzer::Dim;

namespace charon {

//! obtain Crowell-Sze avalanche (impact ionization) generation rate
template<typename EvalT, typename Traits>
class Avalanche_CrowellSze
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Avalanche_CrowellSze(
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
  PHX::MDField<ScalarT,Cell,Point> avalanche_rate;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp_h;

  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_pot;

  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field_h;

  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_drForce;
  PHX::MDField<const ScalarT,Cell,Point,Dim> hole_drForce;

  PHX::MDField<const ScalarT,Cell,Point,Dim> curr_dens_e;
  PHX::MDField<const ScalarT,Cell,Point,Dim> curr_dens_h;

  PHX::MDField<const ScalarT,Cell,Point> dens_e;
  PHX::MDField<const ScalarT,Cell,Point> dens_h;

  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // lattice temperature

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double J0;      // current density scaling
  double R0;      // recomb./gen. scaling
  double E0;      // electric field scaling
  double C0;      // concentration scaling
  double T0;      // temperature scaling

  int num_points, num_dims;
  int num_nodes;
  std::size_t basis_index;
  std::string basis_name;
  bool isSGCVFEM;

  // Crowell-Sze model parameters
  double lambda300_e, Ei_e;
  double lambda300_h, Ei_h;
  double E_opt_ph;

  // minimum field below which avalanche generation = 0
  double minField;

  // driving force name
  std::string driveForce;

  // driving force damping parameters [cm^-3]
  double eDrForceRefDens, hDrForceRefDens;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  /**
   * @brief Initialize Avalanche Crowell-Sze parameters
   */
  void initAvaParams(const std::string& matName, const Teuchos::ParameterList& avaParamList);

}; // end of class Avalanche_CrowellSze


}

#endif
