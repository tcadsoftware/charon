
#ifndef CVFEM_INTEGRATOR_SUBCVFLUXDOTNORM_DECL_HPP
#define CVFEM_INTEGRATOR_SUBCVFLUXDOTNORM_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Intrepid2_Basis.hpp"


using panzer::Cell;
using panzer::Dim;
using panzer::BASIS;
using panzer::Edge;


namespace charon {

// Integrate over a subcontrol volume surface \int_cv_surf F \cdot n dS
template<typename EvalT, typename Traits>
class Integrator_SubCVFluxDotNorm
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Integrator_SubCVFluxDotNorm(
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

  // output
  PHX::MDField<ScalarT,Cell,BASIS> residual;

  // input
  PHX::MDField<const ScalarT,Cell,Edge,Dim> subcv_edge_flux;

  std::string basis_name;
  std::size_t basis_index;

  std::size_t int_rule_index;
  std::size_t int_rule_degree;

  int num_nodes;
  int num_edges;
  int num_dims;
  int num_ips;

  double multiplier;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;


}; // end of class Integrator_SubCVFluxDotNorm


}


#endif
