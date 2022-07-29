
#ifndef CHARON_SGCVFEM_CENTROID_DRIVEFORCE_DECL_HPP
#define CHARON_SGCVFEM_CENTROID_DRIVEFORCE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::IP;
//using panzer::BASIS;
using panzer::Edge;
using panzer::Dim;

namespace charon {

// Evaluate the effective electric field Fn\p,eff at IPs (the
// centroids of a subcontrol volumes, and the gradient of quasi
// fermi potential grad(qfp) = -\+ grad(n\p)/(n\p) - Fn\p,eff


template<typename EvalT, typename Traits>
class SGCVFEM_CentroidDriveForce
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    SGCVFEM_CentroidDriveForce(
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
  PHX::MDField<ScalarT,Cell,IP,Dim> field;

  // input
  PHX::MDField<const ScalarT,Cell,Point> potential;   // potential in [V0]
  PHX::MDField<const ScalarT,Cell,Point> intrinfermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,Point> bandgap;     // band gap w/o BGN
  PHX::MDField<const ScalarT,Cell,Point> effbandgap;  // effective band gap w/ BGN
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  PHX::MDField<const ScalarT,Cell,Point> density; // carrier density


  PHX::MDField<const ScalarT,Cell,Point> elec_degfactor;
  PHX::MDField<const ScalarT,Cell,Point> hole_degfactor;

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0;
  double T0;

  // for basis points
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  // reference edge length
  double refEdgeLen;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  std::string carrType;
  std::string drForceType;
  int sign;

  // dimensions
  int num_dims;
  int num_edges;
  int num_ips;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class SGCVFEM_CentroidDriveForce


}

#endif
