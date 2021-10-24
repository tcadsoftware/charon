
#ifndef CHARON_FEM_ELECTRICFIELD_DECL_HPP
#define CHARON_FEM_ELECTRICFIELD_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::BASIS;
using panzer::Dim;

// Evaluate the effective electric field Fn\p,eff at IPs, and
// the gradient of quasi fermi potential grad(qfp) = -\+ grad(n\p)/(n\p) - Fn\p,eff

namespace charon {

template<typename EvalT, typename Traits>
class FEM_ElectricField
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    FEM_ElectricField(
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
  PHX::MDField<ScalarT,Cell,Point,Dim> efield;
  PHX::MDField<ScalarT,Cell,Point,Dim> grad_qfp;

  // input
  PHX::MDField<const ScalarT,Cell,Point> intrinfermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,Point> potential;   // potential in [V0]
  PHX::MDField<const ScalarT,Cell,Point> bandgap;     // band gap w/o BGN
  PHX::MDField<const ScalarT,Cell,Point> effbandgap;  // effective band gap w/ BGN

  PHX::MDField<const ScalarT,Cell,Point> density;           // carrier density
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_density;  // gradient of carrier density

  PHX::MDField<const ScalarT,Cell,Point> elec_degfactor;
  PHX::MDField<const ScalarT,Cell,Point> hole_degfactor;

  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0;
  double T0;

  // IP
  std::size_t num_points;
  std::size_t num_dims;

  // BASIS
  std::string basis_name;
  std::size_t basis_index;
  int num_basis;

  // carrier type
  std::string carrType;
  double sign;

  Kokkos::DynRankView<ScalarT,PHX::Device> negEffPot;  // intermediate field

  Kokkos::DynRankView<ScalarT,PHX::Device> negPot;  // for testing
  Kokkos::DynRankView<ScalarT,PHX::Device> gradNegPot;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class FEM_ElectricField


}

#endif
