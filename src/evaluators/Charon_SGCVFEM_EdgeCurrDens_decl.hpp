
#ifndef CHARON_SGCVFEM_EDGECURRDENS_DECL_HPP
#define CHARON_SGCVFEM_EDGECURRDENS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Edge;

namespace charon {

// compute the edge current density (scaled and scalar) in the CVFEM-SG formulation
// evaluated at the center of primary edge

template<typename EvalT, typename Traits>
class SGCVFEM_EdgeCurrDens
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCVFEM_EdgeCurrDens(
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
  PHX::MDField<ScalarT,Cell,Edge> edge_currdens; // edge curr. dens. (scalar)

  // input
  PHX::MDField<const ScalarT,Cell,Edge> diff_coeff;    // diff. coeff.
  PHX::MDField<const ScalarT,Cell,Edge> mobility;      // mobility

  PHX::MDField<const ScalarT,Cell,BASIS> density;      // carrier density

  PHX::MDField<const ScalarT,Cell,BASIS> intrin_fermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> bandgap;      // band gap w/o BGN in [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> eff_bandgap;  // effective band gap w/ BGN in [eV]

  PHX::MDField<const ScalarT,Cell,BASIS> elec_degfac;  // electron degeneracy factor
  PHX::MDField<const ScalarT,Cell,BASIS> hole_degfac;  // hole degeneracy factor
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // lattice temperature

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // voltage scaling in [V]
  double T0; // temperature scaling in [K]

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  std::string carrType;

  // for Edge
  int num_dims;
  int num_edges;

  // different sign for electrons and holes
  double sign;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCVFEM_EdgeCurrDens


}

#endif
