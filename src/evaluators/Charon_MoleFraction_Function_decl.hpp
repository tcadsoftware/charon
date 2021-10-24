
#ifndef CHARON_MOLEFRACTION_FUNCTION_DECL_HPP
#define CHARON_MOLEFRACTION_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include <vector>

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Point;

namespace charon {

  class uniformMoleFracParams
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

  class linearMoleFracParams
  {
    void testcoord (const std::string axis, const Teuchos::ParameterList& plist,
                    double& min, double& max, bool& checkAxis);
    public:
      linearMoleFracParams(): y_min(0.0), y_max(0.0), z_min(0.0), z_max(0.0)
      {;};

      ~linearMoleFracParams(){;}

      void parseLinear (const Teuchos::ParameterList& plist, const int num_dim);

      double startVal;
      double endVal;

      // x direction
      double x_min;
      double x_max;
      bool x_checkAxis;

      // y direction
      double y_min;
      double y_max;
      bool y_checkAxis;

      // z direction
      double z_min;
      double z_max;
      bool z_checkAxis;
  };

template<typename EvalT, typename Traits>
class MoleFraction_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    MoleFraction_Function(
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

  // evaluate total x mole fraction (due to various functions) at given (x,y,z)
  double evaluateMoleFraction(const double& x, const double& y, const double& z);

  // evaluate uniform/constant x mole fraction at given (x,y,z)
  double evalUniformMoleFrac(const double& x, const double& y,
    const double& z, const uniformMoleFracParams& umfp);

  // evaluate linear x mole fraction at given (x,y,z)
  double evalLinearMoleFrac(const double& x, const double& y,
    const double& z, const linearMoleFracParams& lmfp);

  // evaluate a single-axis linear (x, or y, or z)
  double evalSingleLinear(const std::string axis, bool& found, const double& coord,
      const double& min, const double& max, const bool& checkAxis);

  // output
  PHX::MDField<ScalarT,Cell,IP> molefrac; // x mole fraction
  PHX::MDField<ScalarT,Cell,BASIS> molefrac_basis;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;
  int num_basis;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  Teuchos::ParameterList moleFracParamList;

  std::vector<uniformMoleFracParams> umfp_vec;
  std::vector<linearMoleFracParams> lmfp_vec;


}; // end of class MoleFraction_Function


}

#endif
