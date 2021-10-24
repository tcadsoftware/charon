
#ifndef CHARON_IC_FUNCTION_DECL_HPP
#define CHARON_IC_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Doping_Params.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

// Set initial condition to a function profile in space
namespace charon {


  class uniformICParams
  {
    public:
        void parseUniform (const Teuchos::ParameterList& plist);
        double value;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
  };

  class gaussianICParams
  {
    void testcoord (const std::string& axis, const Teuchos::ParameterList& plist,
                    double& width, double& gaussMin, double& gaussMax,
                    double& min, double& max, bool& checkAxis);

    public:
      gaussianICParams() {;}

      ~gaussianICParams() {;}

      void parseGaussian (const Teuchos::ParameterList& plist, const int num_dim);

      double maxValue;
      double minValue;

      // x direction
      double xWidth;
      double xGaussMin;
      double xGaussMax;
      double xMin;
      double xMax;
      bool xCheckAxis;

      // y direction
      double yWidth;
      double yGaussMin;
      double yGaussMax;
      double yMin;
      double yMax;
      bool yCheckAxis;

      // z direction
      double zWidth;
      double zGaussMin;
      double zGaussMax;
      double zMin;
      double zMax;
      bool zCheckAxis;
  };



template<typename EvalT, typename Traits>
class IC_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IC_Function(
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

  // evaluate the total initial condition (IC) at a given (x,y,z)
  double evaluateIC(const double& x, const double& y, const double& z);

  // evaluate uniform IC at a given (x,y,z)
  double evalUniformIC(const double& x, const double& y, const double& z, const uniformICParams& udp);

  // evaluate gaussian IC at a given (x,y,z)
  double evalGaussianIC(const double& x, const double& y, const double& z, const gaussianICParams& gdp);

  // evaluate a single-axis gaussian (x, or y, or z)
  double evalSingleGaussian(const double& coord, const double& minValue,
                            const double& maxValue, const double& width,
                            const double& gaussMin, const double& gaussMax,
                            const double& min, const double& max, const bool& checkAxis);


  // output
  PHX::MDField<ScalarT,Cell,BASIS> carrier_density; // scaled

  // for basis points
  int num_basis;
  int num_dim;
  std::string basis_name;
  std::size_t basis_index;

  std::string dof_name;

  // Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  std::vector<uniformICParams> udp_vec;
  std::vector<gaussianICParams> gdp_vec;

}; // end of class IC_Function


}

#endif
