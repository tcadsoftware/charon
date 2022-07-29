
#ifndef CHARON_DOPINGRAW_FUNCTION_DECL_HPP
#define CHARON_DOPINGRAW_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Doping_Params.hpp"
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Charon_MMS_AnalyticFunctions.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include <ostream>

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

namespace charon {


// Generic Profile Calculator
class ProfileEvals
{
private:
  int num_dim;

public:
  ProfileEvals(int dim);

  // evaluate a single-axis linear (x, or y, or z)
  double evalSingleLinear(const std::string axis, bool& found, 
			  const double& coord,
			  const double& min, const double& max, 
			  const bool& checkAxis);
  // evaluate linear profile at given (x,y,z)
  std::vector<double> evalLinearProfile(
	           const double& x, const double& y,
		   const double& z, const linearDopingParams& ldp);
  // evaluate a single-axis gaussian (x, or y, or z)
  double evalSingleGaussian(const std::string axis, bool& found, 
			    const double& coord, const double& minProfVal,
			    const double& maxProfVal, const double& min, 
			    const double& max, const double& loc,
			    const double& width, const bool& checkAxis, 
			    const std::string& dir);
  // evaluate gaussian profile at given (x,y,z)
  std::vector<double> evalGaussianProfile(
	const double& x, const double& y, const double& z, 
	const gaussianDopingParams& gdp);
  // evaluate a single-axis erfc (x, or y, or z)
  double evalSingleErfc(const std::string axis, bool& found, 
			const double& coord, const double& minProfVal,
			const double& maxProfVal, const double& min, 
			const double& max, const double& loc,
			const double& width, const bool& checkAxis, 
			const std::string& dir);
  // evaluate erfc profile at given (x,y,z)
  std::vector<double> evalErfcProfile(
        const double& x, const double& y, const double& z, 
	const erfcDopingParams& edp);
  // evaluate a single-axis MGauss (x, or y, or z)
  double evalSingleMGauss(const std::string axis, bool& found, 
			  const double& coord, const double& minProfVal, 
			  const double& maxProfVal, const double& min, 
			  const double& max, const bool& checkErfc,
			  const double& width, const bool& checkAxis);
  // evaluate gaussian (=MGauss) profile at given (x,y,z)
  std::vector<double> evalMGaussProfile(
	const double& x, const double& y, const double& z, 
        const mgaussDopingParams& mgdp);
  // evaluate halo profile
  std::vector<double> evalHaloProfile(
        const double& x, const double& y, const double& z, 
        const haloDopingParams& hdp);
};



//! obtain a uniform doping
template<typename EvalT, typename Traits>
class DopingRaw_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DopingRaw_Function(
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

  void preEvaluate(typename Traits::PreEvalData d);

  // evaluate total doping (due to various doping functions) at given (x,y,z)
  std::vector<double> evaluateDoping(const double& x, const double& y, const double& z);

  // evaluate uniform doping at given (x,y,z)
  std::vector<double> evalUniformDoping(const double& x, const double& y,
    const double& z, const uniformDopingParams& udp, int udCounter, const Teuchos::ParameterList& plist);

  // evaluate gaussian doping at given (x,y,z)
  std::vector<double> evalGaussianDoping(const double& x, const double& y,
    const double& z, const gaussianDopingParams& gdp);

  // evaluate a single-axis gaussian (x, or y, or z)
  double evalSingleGaussian(const std::string axis, bool& found, const double& coord, const double& minDopVal,
   const double& maxDopVal, const double& min, const double& max, const double& loc,
   const double& width, const bool& checkAxis, const std::string& dir);

  // evaluate linear doping at given (x,y,z)
  std::vector<double> evalLinearDoping(const double& x, const double& y,
    const double& z, const linearDopingParams& ldp);

  // evaluate a single-axis linear (x, or y, or z)
  double evalSingleLinear(const std::string axis, bool& found, const double& coord,
      const double& min, const double& max, const bool& checkAxis);

  // evaluate erfc doping at given (x,y,z)
  std::vector<double> evalErfcDoping(const double& x, const double& y,
    const double& z, const erfcDopingParams& edp);

  // evaluate a single-axis erfc (x, or y, or z)
  double evalSingleErfc(const std::string axis, bool& found, const double& coord, const double& minDopVal,
   const double& maxDopVal, const double& min, const double& max, const double& loc,
   const double& width, const bool& checkAxis, const std::string& dir);

  // evaluate gaussian (=MGauss) doping at given (x,y,z)
  std::vector<double> evalMGaussDoping(const double& x, const double& y, const double& z, const mgaussDopingParams& mgdp);

  // evaluate a single-axis MGauss (x, or y, or z)
  double evalSingleMGauss(const std::string axis, bool& found, const double& coord,
   const double& minDopVal, const double& maxDopVal,
   const double& min, const double& max, const bool& checkErfc,
   const double& width, const bool& checkAxis);

  //evaluate halo doping
  std::vector<double> evalHaloDoping(const double& x, const double& y,
    const double& z, const haloDopingParams& hdp);

  // export 2D doping from an external file and map doping values to given (x,y,z) (Note: doping along z-axis is assumed to be homogeneous)
  std::vector<double> evalFile2DDoping(int counter, const double& x, const double& y,
                                       const double& z, const Teuchos::ParameterList& plist);

  // export 1D doping from an external file and map doping values to given (x,y,z)
  std::vector<double> evalFile1DDoping(int counter, const double& x, const double& y,
    const double& z, const Teuchos::ParameterList& plist);

  // export 2D doping from an external file and map doping values to given (x,y,z)
  std::vector<double> evalFile3DDoping(int counter, const double& x, const double& y,
                                       const double& z, const Teuchos::ParameterList& plist);

  // This is needed for evalFile2DDoping
  void readDopingFile2D(const Teuchos::ParameterList& plist);

  // This is needed for evalFile3DDoping
  void readDopingFile3D(const Teuchos::ParameterList& plist);

  // This is needed for evalFile1DDoping
  void readDopingFile1D(const Teuchos::ParameterList& plist);

  // Read in Gauss decay parameters 
  void readGaussDecayParams(int type, int counter, const Teuchos::ParameterList& plist);   

  // Compute Gauss decay factor to be added to any doping read from file
  double evalGaussDecayFactor(int type, int counter, const double& x, const double& y,
    const double& z);

  //This is needed for manufactured solution doping
  Teuchos::RCP<charon::MMS_DD_RDH_1_AnalyticFunction> mms1;
  Teuchos::RCP<charon::MMS_DD_RDH_2_AnalyticFunction> mms2;
  Teuchos::RCP<charon::MMS_NLP_GLH_1_AnalyticFunction> mms_nlp_func;

  class doping_struct
  {
  public:
    doping_struct():
      x(0.0),
      y(0.0),
      z(0.0),
      d(0.0)
    {
    }
    doping_struct(const doping_struct& dS):
      x(dS.x),
      y(dS.y),
      z(dS.z),
      d(dS.d)
    {
    }

    double x,y,z,d;

    inline bool operator < (const doping_struct &a) const
    { return ( (x<a.x) ||
               ((x==a.x) && (y<a.y)) ||
               ((x==a.x) && (y==a.y) && (z<a.z))); }

    inline bool operator == (const doping_struct &a) const
    { return ( (x==a.x) && (y==a.y) && (z==a.z) ); }

    inline double distance(const doping_struct &a) const
    { return ( sqrt((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y) +(z-a.z)*(z-a.z)) ); }

  };

  std::vector< std::vector<doping_struct> >DV;
  std::vector<double> MinX, MaxX, MinY, MaxY, MinZ, MaxZ; // external doping bounds (automatically determined when read external doping values

  // This struct is needed for evalFile1DDoping
  struct doping_struct_1D
  {
    double x,d;
    inline bool operator < (const doping_struct_1D &a) const
    { return ( (x<a.x)); }
    inline bool operator == (const doping_struct_1D &a) const
    { return ( (x==a.x)); }
    inline double distance(const doping_struct_1D &a) const
    { return ( fabs(x-a.x) ); }
  };

  std::vector<std::vector<doping_struct_1D> > DV1;

  // vectors to store gauss decay parameters 
  // The first vector size = 2 (0 for File1D, 1 for Uniform)
  // The second vector size = number of File1D or Uniform doping profiles
  // The third vector size = number of gauss decay axises (6 total) for each File1D or Uniform doping
  std::vector<std::vector<std::vector<std::string> > > decayDir;  
  std::vector<std::vector<std::vector<double> > > decayPos; 
  std::vector<std::vector<std::vector<double> > >decayWidth; 
  
  bool sweepingIsOn;

  // input
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]

  // output
  PHX::MDField<ScalarT,Cell,IP> doping_raw; // net doping @ IPs
  PHX::MDField<ScalarT,Cell,IP> acceptor_raw;
  PHX::MDField<ScalarT,Cell,IP> donor_raw;

  PHX::MDField<ScalarT,Cell,BASIS> doping_raw_basis; // @ BASIS points
  PHX::MDField<ScalarT,Cell,BASIS> acceptor_raw_basis;
  PHX::MDField<ScalarT,Cell,BASIS> donor_raw_basis;

  // Fields for each workset to store values
  std::vector<PHX::MDField<ScalarT,Cell,IP> > acceptor_raw_wkst;
  std::vector<PHX::MDField<ScalarT,Cell,IP> > donor_raw_wkst;
  std::vector<PHX::MDField<ScalarT,Cell,BASIS> > acceptor_raw_basis_wkst;
  std::vector<PHX::MDField<ScalarT,Cell,BASIS> > donor_raw_basis_wkst;

  // workset id
  int worksetId;

  // boolean indicating whether to compute doping at each
  // evaluation or store for the entire workset
  bool store_wkst_doping;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  Teuchos::ParameterList dopParamList;

  std::vector<uniformDopingParams> udp_vec;
  std::vector<gaussianDopingParams> gdp_vec;
  std::vector<linearDopingParams> ldp_vec;
  std::vector<erfcDopingParams> edp_vec;
  std::vector<mgaussDopingParams> mgdp_vec;
  std::vector<haloDopingParams> hdp_vec;

  //Homotopy stuff
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;

  Teuchos::RCP<ProfileEvals> prof_eval;

}; // end of class DopingRaw_Function


}

#endif
