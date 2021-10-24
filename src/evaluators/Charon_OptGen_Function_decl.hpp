
#ifndef CHARON_OPTGEN_FUNCTION_DECL_HPP
#define CHARON_OPTGEN_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_empiricalConvolution.hpp"
#include "Charon_Scaling_Parameters.hpp"


using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

namespace charon {

  class gaussianOptGenParams
  {
    void testcoord (const std::string axis, const Teuchos::ParameterList& plist,
                    std::string& dir, double& min, double& max,
                    double& loc, double& width, bool& checkAxis);
    public:
      gaussianOptGenParams(): y_loc(0.0), y_width(0.0), y_min(0.0), y_max(0.0),
                              z_loc(0.0), z_width(0.0), z_min(0.0), z_max(0.0)
      {;};

      ~gaussianOptGenParams(){;}

      void parseGaussian (const Teuchos::ParameterList& plist, const int num_dim);

      double maxVal;
      double minVal;

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

      // time direction
      std::string t_dir;
      double t_loc;
      double t_width;
      double t_min;
      double t_max;
      bool t_checkAxis;
  };


template<typename EvalT, typename Traits>
class OptGen_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    OptGen_Function(
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

  // evaluate total opt. gen. (due to various opt. gen. functions) at given (x,y,z,t)
  double evaluateOptGen(const double& x, const double& y, const double& z, const double& t);

  // evaluate gaussian opt. gen. at given (x,y,z,t)
  double evalGaussianOptGen(const double& x, const double& y, const double& z,
                            const double& t, const gaussianOptGenParams& gp);

  // evaluate a single-axis gaussian (x, or y, or z, or t)
  double evalSingleGaussian(const std::string& axis, bool& found, const double& coord,
         const double& minVal, const double& maxVal, const double& min, const double& max,
         const double& loc, const double& width, const bool& checkAxis, const std::string& dir);

  // export Optgen from an external file and map generation values to given (x,y,z,t)
  double evalFileOptGen(const double& x, const double& y, const double& z, const double& t, const Teuchos::ParameterList& plist);

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]
  double R0; // recomb./gen. scaling,[#/(cm^3.s)]

  // output
  PHX::MDField<ScalarT,Cell,IP> optgen; // opt. gen. @ IPs
  PHX::MDField<ScalarT,Cell,BASIS> optgen_basis; // @ BASIS points

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  Teuchos::ParameterList genParamList;

  std::vector<gaussianOptGenParams> gp_vec;



  // This structures are needed for evalFileOptGen
  void readOptGenFile(const Teuchos::ParameterList& plist);

  struct optgen_struct
  {
    double t,g;
    inline bool operator < (const optgen_struct &a) const
    { return ( (t<a.t)); }
    inline bool operator == (const optgen_struct &a) const
    { return ( (t==a.t)); }
    inline double distance(const optgen_struct &a) const
    { return ( fabs(t-a.t) ); }
  };

  std::vector<optgen_struct> GV;

  Teuchos::RCP<charonSpline> ogSpline;

  std::vector<double> s_x, s_y;




}; // end of class OptGen_Function


}

#endif
