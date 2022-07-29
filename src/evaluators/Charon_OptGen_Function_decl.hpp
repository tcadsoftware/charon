
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
    gaussianOptGenParams(): maxVal(1.0), minVal(1.0), y_loc(0.0), y_width(0.0), y_min(0.0),
                            y_max(0.0), z_loc(0.0), z_width(0.0), z_min(0.0), z_max(0.0)
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
  void readOptGenFile(const Teuchos::ParameterList& plist);
  double evalFileOptGen(const double& x, const double& y, const double& z, const double& t, const Teuchos::ParameterList& plist);

  // read in time vs. optgen (t, g) data
  void readTimeFile1D(const std::string& timeFile); 

  // read in 1D spatial profie for optgen, (x, factor) data
  void readSpaceFile1D(const Teuchos::ParameterList& plist);

  // read in 2D spatial profie for optgen, (x, y, factor) data
  void readSpaceFile2D(const Teuchos::ParameterList& plist);

  // evaluate opt. gen. at given (x,y,z,t) using (t,g) data from Time File
  double evaluateTimeFile1D(int counter, const double& x, const double& y, const double& z, 
                            const double& t, const Teuchos::ParameterList& plist);

  // evaluate opt. gen. spatial factor at given (x,y,z) using (x,g) data from Space File1D
  double evaluateSpaceFile1D(int counter, const double& x, const double& y, const double& z, 
                             const Teuchos::ParameterList& plist);

  // evaluate opt. gen. spatial factor at given (x,y,z) using (x,y,g) data from Space File2D
  double evaluateSpaceFile2D(int counter, const double& x, const double& y, const double& z, 
                             const Teuchos::ParameterList& plist);


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


  // This struct is to store optgen 1D time profile (t, g) with g = opt. gen. in cm^(-3).s^(-1)
  struct optgen_time_1D
  {
    double t, g;

    inline bool operator < (const optgen_time_1D &a) const
    { return ( (t < a.t)); }

    inline bool operator == (const optgen_time_1D &a) const
    { return ( (t == a.t)); }

    inline double distance(const optgen_time_1D &a) const
    { return ( fabs(t-a.t) ); }
  };
  std::vector<std::vector<optgen_time_1D> > OGTime1D;


  // This struct is to store optgen 1D spatial profile (x, g) with g = a scaling factor
  struct optgen_space_1D
  {
    double x, g;

    inline bool operator < (const optgen_space_1D &a) const
    { return ( (x < a.x)); }

    inline bool operator == (const optgen_space_1D &a) const
    { return ( (x == a.x)); }

    inline double distance(const optgen_space_1D &a) const
    { return ( fabs(x-a.x) ); }
  };
  std::vector<std::vector<optgen_space_1D> > OGSpace1D;


  // This struct is to store optgen 2D spatial profile, (x, y, g) with g = a scaling factor
  struct optgen_space_2D
  {
    double x, y, g;

    inline bool operator < (const optgen_space_2D &a) const
    { return ( (x < a.x) || ((x == a.x) && (y < a.y))); }

    inline bool operator == (const optgen_space_2D &a) const
    { return ( (x == a.x) && (y == a.y));  }
    
    inline double distance(const optgen_space_2D &a) const
    { return ( sqrt((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)) );  }

  };
  std::vector<std::vector<optgen_space_2D> > OGSpace2D;

  // store external coorindate bounds (automatically determined when reading external files)
  std::vector<double> MinX2D, MaxX2D, MinY2D, MaxY2D; 

  // boolean indicating whether to compute spatial optgen profile at each
  // evaluation or store for the entire workset
  bool store_wkst_optgen = false;


}; // end of class OptGen_Function


}

#endif
