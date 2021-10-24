
#ifndef CHARON_SGCVFEM_SUBCVCENTROIDVECTOR_DECL_HPP
#define CHARON_SGCVFEM_SUBCVCENTROIDVECTOR_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Edge;
using panzer::Dim;

namespace charon {

// using edge basis functions to interpolate the edge current density
// to the centroid of a subcontrol volume

template<typename EvalT, typename Traits>
class SGCVFEM_SubCVCentroidVector
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCVFEM_SubCVCentroidVector(
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

  // output - at subcontrol volume centroids
  PHX::MDField<ScalarT,Cell,IP,Dim> subcv_centvec;

  // input
  PHX::MDField<const ScalarT,Cell,Edge> edge_currdens;

  // for basis points
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  std::string carrType;

  // dimensions
  int num_dims;
  int num_edges;
  int num_ips;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCVFEM_SubCVCentroidVector


}

#endif
