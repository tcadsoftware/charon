
#ifndef CHARON_IC_GAUSS_DECL_HPP
#define CHARON_IC_GAUSS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

// Set initial condition to a Gaussian profile in space
namespace charon {

template<typename EvalT, typename Traits>
class IC_Gauss
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IC_Gauss(
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
  PHX::MDField<ScalarT,Cell,BASIS> carrier_density; // scaled

  void initialize(const Teuchos::ParameterList& plist);

  void testcoord(const std::string& axis, const Teuchos::ParameterList& plist,
                 double& width, double& gaussMin, double& gaussMax,
                 double& min, double& max, bool& checkAxis);

  // evaluate the Gaussian initial condition at a given (x,y,z)
  double evaluateGaussIC(const double& x, const double& y, const double& z);

  // evaluate a single-axis Gauss (x, or y, or z)
  double evalSingleGauss(const double& coord, const double& width,
                         const double& gaussMin, const double& gaussMax,
                         const double& min, const double& max, const bool& checkAxis);

  // for basis points
  int num_basis;
  int num_dim;
  std::string basis_name;
  std::size_t basis_index;

  std::string dof_name;

  double maxValue, minValue;
  double xWidth, xGaussMin, xGaussMax, xMin, xMax;
  double yWidth, yGaussMin, yGaussMax, yMin, yMax;
  double zWidth, zGaussMin, zGaussMax, zMin, zMax;

  bool xAxis, yAxis, zAxis;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class IC_Gauss


}

#endif
