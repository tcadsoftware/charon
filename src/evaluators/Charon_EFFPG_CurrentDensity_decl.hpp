
#ifndef CHARON_EFFPG_CURRENTDENSITY_DECL_HPP
#define CHARON_EFFPG_CURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

#include "Intrepid2_Basis.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;
using panzer::BASIS;
using panzer::Edge;

namespace charon {

//! calculate the EFF-PG version of current density (Scharfetter-Gummel-like flux)
template<typename EvalT, typename Traits>
class EFFPG_CurrentDensity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    EFFPG_CurrentDensity(
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
  PHX::MDField<ScalarT,Cell,IP,Dim> current_density; // @ IPs

  // input
  PHX::MDField<const ScalarT,Cell,Edge> mobility;     // scaled
  PHX::MDField<const ScalarT,Cell,Edge> diff_coeff;
  PHX::MDField<const ScalarT,Cell,BASIS> carrier_density;

  PHX::MDField<const ScalarT,Cell,BASIS> intrin_fermi; // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> bandgap;
  PHX::MDField<const ScalarT,Cell,BASIS> effbandgap;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]

  // IP
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ips;
  int num_dims;

  // BASIS
  std::string basis_name;
  std::size_t basis_index;
  int num_nodes;
  int num_edges;

  // HCURL BASIS
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  // carrier type
  std::string carrType;

  // diffusion coefficient sign
  double diffSign;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class EFFPG_CurrentDensity


}

#endif
