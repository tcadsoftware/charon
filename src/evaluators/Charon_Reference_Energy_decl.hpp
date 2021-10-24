
#ifndef CHARON_REFERENCE_ENERGY_DECL_HPP
#define CHARON_REFERENCE_ENERGY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

/**
 * @brief Compute the Reference Energy in [eV], which is the intrinsic Fermi
 * energy level of the Reference Material from the vacuum level, i.e.,
 * Eref(T0) = Chi(T0) + 0.5*Eg(T0) + 0.5*kb*T0*log(Nc(T0)/Nv(T0)),
 * where T0 = Temperature Scaling).
 * Eref is needed for hetero-geneous devices and valid for homogeneous devices.
 *
 * Added by Suzey on 10/25/2017:
 * By using T0 (Temperature Scaling), not TL (Lattice Temperature), Eref(T0)
 * becomes valid for both isothermal and non-isothermal DD simulations.
 * By using Eref(T0), the potential solved from the Poisson equation is
 * the vacuum potential shifted by Eref(T0) in both isothermal
 * and non-isothermal DD simulations.
 */

//! calculate reference energy Eref needed for heterogenous devices
template<typename EvalT, typename Traits>
class Reference_Energy
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Reference_Energy(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> ref_energy;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // temperature scaling in [K]

  // material parameters
  double Eg300, alpha, beta, Chi300; // Chi300 = electron affinity at 300 K
  double Nc300, Nv300, Nc_F, Nv_F;
  double constBg, constChi;          // const band gap

  bool isBgConst, isChiConst;

  int numPoints;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Reference_Energy


}

#endif
