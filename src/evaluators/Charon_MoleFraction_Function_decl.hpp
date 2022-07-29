
#ifndef CHARON_MOLEFRACTION_FUNCTION_DECL_HPP
#define CHARON_MOLEFRACTION_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include <vector>
#include "Charon_DopingRaw_Function.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Point;

namespace charon {

  class uniformMoleFracParams
  {
    public:
    void parseUniform(const Teuchos::ParameterList& plist, const std::string& matArity);
    double value;
    double value1;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    std::string mat_arity;
  };

  class linearMoleFracParams
  {
    void testcoord (const std::string axis, const Teuchos::ParameterList& plist,
                    double& min, double& max, bool& checkAxis);
    public:
      linearMoleFracParams(): y_min(0.0), y_max(0.0), z_min(0.0), z_max(0.0) {};

      ~linearMoleFracParams() {}

      void parseLinear(const Teuchos::ParameterList& plist, const int num_dim, 
		       const std::string& matArity);

      double startVal;
      double endVal;
      double startVal1;
      double endVal1;

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

      std::string mat_arity;
  };


  class gaussMoleFracParams
  {
    void testcoord(const std::string axis, const Teuchos::ParameterList& plist,
                    std::string& dir, double& min, double& max,
                    double& loc, double& width, bool& checkAxis);
    public:
      gaussMoleFracParams(): y_loc(0.0), y_width(0.0), y_min(0.0), y_max(0.0),
			     z_loc(0.0), z_width(0.0), z_min(0.0), z_max(0.0) {};

      ~gaussMoleFracParams(){};

      void parseGaussian(const Teuchos::ParameterList& plist, 
			 const int num_dim, const std::string& matArity);

      double maxVal;
      double minVal;
      double maxVal1;
      double minVal1;
      // x direction
      std::string x_dir;
      double x_loc;
      double x_width;
      double x_min;
      double x_max;
      bool x_checkAxis;
      // y direction
      std::string y_dir;
      double y_loc;
      double y_width;
      double y_min;
      double y_max;
      bool y_checkAxis;
      // z direction
      std::string z_dir;
      double z_loc;
      double z_width;
      double z_min;
      double z_max;
      bool z_checkAxis;

      std::string mat_arity;
  };


  class erfcMoleFracParams
  {
    void testcoord(const std::string axis, const Teuchos::ParameterList& plist,
                    std::string& dir, double& min, double& max,
                    double& loc, double& width, bool& checkAxis);
    public:
      erfcMoleFracParams(): y_loc(0.0), y_width(0.0), y_min(0.0), y_max(0.0),
                              z_loc(0.0), z_width(0.0), z_min(0.0), z_max(0.0) {};

      ~erfcMoleFracParams(){};

    void parseErfc(const Teuchos::ParameterList& plist, 
		   const int num_dim, const std::string& matArity);

      double maxVal;
      double minVal;
      double maxVal1;
      double minVal1;
      // x direction
      std::string x_dir;
      double x_loc;
      double x_width;
      double x_min;
      double x_max;
      bool x_checkAxis;
      // y direction
      std::string y_dir;
      double y_loc;
      double y_width;
      double y_min;
      double y_max;
      bool y_checkAxis;
      // z direction
      std::string z_dir;
      double z_loc;
      double z_width;
      double z_min;
      double z_max;
      bool z_checkAxis;

      std::string mat_arity;
  };  

 
  class mgaussMoleFracParams
  {
    void testcoord(const std::string axis, const Teuchos::ParameterList& plist,
                    std::string& dir, double& min, double& max,
                    bool& checkErfc, double& width, bool& checkAxis);
    public:
      mgaussMoleFracParams(): y_width(0.0), y_min(0.0), y_max(0.0),
                              z_width(0.0), z_min(0.0), z_max(0.0) {};

      ~mgaussMoleFracParams(){};

      void parseMGauss(const Teuchos::ParameterList& plist, 
		       const int num_dim, const std::string& matArity);

      double maxVal;
      double minVal;
      double maxVal1;
      double minVal1;
      // x direction
      std::string x_dir;
      double x_width;
      double x_min;
      double x_max;
      bool x_checkErfc;
      bool x_checkAxis;
      // y direction
      std::string y_dir;
      double y_width;
      double y_min;
      double y_max;
      bool y_checkErfc;
      bool y_checkAxis;
      // z direction
      std::string z_dir;
      double z_width;
      double z_min;
      double z_max;
      bool z_checkErfc;
      bool z_checkAxis;

      std::string mat_arity;
  };  
  

  class haloMoleFracParams
  {
    void testcoord(const std::string axis, const Teuchos::ParameterList& plist);

    public:
      haloMoleFracParams(): x_center(0.0), y_center(0.0), z_center(0.0),r1(0.0), 
			    r2(0.0), rotation(0.0) {};

      ~haloMoleFracParams(){};

      void parseHalo(const Teuchos::ParameterList& plist, 
		       const int num_dim, const std::string& matArity);

      
      std::string distributionType;
      double val;
      double minVal;
      double val1;
      double minVal1;
      double width;
      // x direction
      double x_center;
      bool x_checkAxis;
      // y direction 
      double y_center;
      bool y_checkAxis;
      // z direction
      double z_center;
      bool z_checkAxis;
      // ellipsoid axes
      double r1,r2;
      double rotation;

      std::string mat_arity;
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
  double evaluate_xMoleFraction(const double& x, const double& y, const double& z);

  // evaluate total y mole fraction (due to various functions) at given (x,y,z)
  double evaluate_yMoleFraction(const double& x, const double& y, const double& z);

  // evaluate uniform/constant x mole fraction at given (x,y,z)
  double evalUniform_xMoleFrac(const double& x, const double& y,
    const double& z, const uniformMoleFracParams& umfp);

  // evaluate linear x mole fraction at given (x,y,z)
  double evalLinear_xMoleFrac(const double& x, const double& y,
    const double& z, const linearMoleFracParams& lmfp);

  // evaluate uniform/constant y mole fraction at given (x,y,z)
  double evalUniform_yMoleFrac(const double& x, const double& y,
    const double& z, const uniformMoleFracParams& umfp);

  // evaluate linear y mole fraction at given (x,y,z)
  double evalLinear_yMoleFrac(const double& x, const double& y,
    const double& z, const linearMoleFracParams& lmfp);

  // evaluate a single-axis linear (x, or y, or z)
  double evalSingleLinear(const std::string axis, bool& found, const double& coord,
      const double& min, const double& max, const bool& checkAxis);

  // output
  PHX::MDField<ScalarT,Cell,IP> molefrac; // x mole fraction
  PHX::MDField<ScalarT,Cell,BASIS> molefrac_basis;
  PHX::MDField<ScalarT,Cell,IP> xMoleFrac; // x mole fraction
  PHX::MDField<ScalarT,Cell,BASIS> xMoleFrac_basis;
  PHX::MDField<ScalarT,Cell,IP> yMoleFrac; // y mole fraction
  PHX::MDField<ScalarT,Cell,BASIS> yMoleFrac_basis;

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
  bool withMoleFrac;
  std::string mat_arity;

  std::vector<uniformMoleFracParams> umfp_vec;
  std::vector<linearMoleFracParams> lmfp_vec;
  std::vector<gaussianDopingParams> gmfp_equiv_vec;
  std::vector<gaussianDopingParams> gmfp_equiv_vec1;
  std::vector<erfcDopingParams> emfp_equiv_vec;
  std::vector<erfcDopingParams> emfp_equiv_vec1;
  std::vector<mgaussDopingParams> mmfp_equiv_vec;
  std::vector<mgaussDopingParams> mmfp_equiv_vec1;
  std::vector<haloDopingParams> hmfp_equiv_vec;
  std::vector<haloDopingParams> hmfp_equiv_vec1;

  Teuchos::RCP<ProfileEvals> prof_eval;
}; // end of class MoleFraction_Function


}

#endif
