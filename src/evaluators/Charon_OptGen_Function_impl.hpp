
#ifndef CHARON_OPTGEN_FUNCTION_IMPL_HPP
#define CHARON_OPTGEN_FUNCTION_IMPL_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  parseGaussian()
//
///////////////////////////////////////////////////////////////////////////////
void gaussianOptGenParams::parseGaussian (const Teuchos::ParameterList& plist,
                                          const int num_dim)
{
  using std::string;
  using Teuchos::ParameterList;

  maxVal = plist.get<double>("Max Value");
  minVal = plist.get<double>("Min Value");

  if ((maxVal < 0.) || (minVal < 0.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gaussian Optical Generation Max and Min Values must be greater than 0.");

  testcoord("X", plist, x_dir, x_min, x_max, x_loc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
    testcoord("Z", plist, y_dir, z_min, z_max, z_loc, z_width, z_checkAxis);
  }
  testcoord("Time", plist, t_dir, t_min, t_max, t_loc, t_width, t_checkAxis);
}


///////////////////////////////////////////////////////////////////////////////
//
//  testcoord()
//
///////////////////////////////////////////////////////////////////////////////
void gaussianOptGenParams::testcoord (const std::string axis,
  const Teuchos::ParameterList& plist, std::string& dir, double& min,
  double& max, double& loc, double& width, bool& checkAxis)
{
  using std::string;
  using Teuchos::ParameterList;

  // axis+" Peak Location" and axis+" Width" must be specified together
  if (plist.isParameter(axis+" Peak Location") && !plist.isParameter(axis+" Width"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! " << axis << " Peak Location must be specified together with "
      << axis << " Width !" << std::endl);

  if (!plist.isParameter(axis+" Peak Location") && plist.isParameter(axis+" Width"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! " << axis << " Width must be specified together with "
      << axis << " Peak Location !" << std::endl);

  checkAxis = false;
  if (plist.isParameter(axis+" Peak Location") && plist.isParameter(axis+" Width"))
  {
    loc = plist.get<double>(axis+" Peak Location");
    width = plist.get<double>(axis+" Width");

    checkAxis = true;
    dir = "Both";  // both axis directions by default

    // check if X/Y/Z Direction is specified (Both, Positive, or Negative)
    if (plist.isParameter(axis+" Direction")) dir = plist.get<string>(axis+" Direction");

    // check if X/Y/Zmin and/or X/Y/Zmax are specified
    max = 1e100;
    min = -1e100;
    if ( plist.isParameter(axis+" Min") )  min = plist.get<double>(axis+" Min");
    if ( plist.isParameter(axis+" Max") )  max = plist.get<double>(axis+" Max");
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
OptGen_Function<EvalT, Traits>::
OptGen_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();

  // opt. gen. parameterlist
  genParamList = p.sublist("Optical Generation ParameterList");

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  R0 = scaleParams->scale_params.R0;

  // evaluated fields
  optgen = MDField<ScalarT,Cell,IP>(n.field.opt_gen,scalar);
  optgen_basis = MDField<ScalarT,Cell,BASIS>(n.field.opt_gen,basis_scalar);
  this->addEvaluatedField(optgen);
  this->addEvaluatedField(optgen_basis);

  std::string name = "OptGen_Function";
  this->setName(name);

  for (ParameterList::ConstIterator model_it = genParamList.begin();
       model_it != genParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      if (funcType == "File") readOptGenFile(funcParamList);

      if ( (funcType == "Gauss") || (funcType == "Gaussian") )
      {
        gaussianOptGenParams gp_;
        gp_.parseGaussian(funcParamList,num_dim);
        gp_vec.push_back(gp_);
      }
    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop
}


///////////////////////////////////////////////////////////////////////////////
//
//  readOptGenFile()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void OptGen_Function<EvalT, Traits>::readOptGenFile(const Teuchos::ParameterList& plist)
{
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::string;

   //read OptGen profile from an external file and assign it to FileDoping array
   const string FileName = plist.get<string>("File Name");
   optgen_struct r;

   ifstream in(FileName.c_str());

   if (!in)
     TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
     << "Error ! Cannot read OptGen file '" << FileName << "'" << std::endl);

   while( in >> r.t >> r.g ) GV.push_back(r);

   // sort the OptGen points by t
   sort(GV.begin(),GV.end());
   // eliminate duplicate (t,g) entries
   GV.resize( unique (GV.begin(), GV.end()) - GV.begin() );

   //initiate Spline interpolation parameters here

   ogSpline = Teuchos::rcp(new charon::charonSpline());
   s_x.resize(GV.size());
   s_y.resize(GV.size());
   for (size_t i=1; i<=GV.size()-1; i++)
   {
        s_x[i] = GV[i].t;
        s_y[i] = GV[i].g;
   }
   ogSpline->createSpline(s_x,s_y);

}







///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
OptGen_Function<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
OptGen_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  size_type num_basis = optgen_basis.dimension(1);

  // get the present time in [s]. Note that input times are also in [s].
  double t = workset.time * t0;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // optical generation at IPs
    for (int ip = 0; ip < num_ip; ++ip)
    {
      double x = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,0);
      double y = 0.0, z = 0.0;
      if (num_dim == 2)
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
      if (num_dim == 3)
      {
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
        z = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,2);
      }

      // evaluate the opt. gen.
      double genValue = evaluateOptGen(x,y,z,t);
      optgen(cell,ip) = genValue / R0;
    }

    // optical generation at basis points
    for (size_type basis = 0; basis < num_basis; ++basis)
    {
      double x = (workset.bases[basis_index])->basis_coordinates(cell,basis,0);
      double y = 0.0, z = 0.0;
      if (num_dim == 2)
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
      if (num_dim == 3)
      {
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
        z = (workset.bases[basis_index])->basis_coordinates(cell,basis,2);
      }

      double genValue = evaluateOptGen(x,y,z,t);
      optgen_basis(cell,basis) = genValue / R0;
    }

  } // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateOptGen()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double OptGen_Function<EvalT, Traits>::evaluateOptGen(const double& x,
    const double& y, const double& z, const double& t)
{
  using std::string;
  using Teuchos::ParameterList;

  double genVal = 0.0;
  double tmpVal = 0.0;

  for (std::size_t i = 0; i < gp_vec.size(); ++i)
  {
    tmpVal = evalGaussianOptGen(x,y,z,t,gp_vec[i]);
    genVal += tmpVal;
  }

// This section is needed to process OptGen from non-analytical (external) files
  for (ParameterList::ConstIterator model_it = genParamList.begin();
       model_it != genParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");


      if ( (funcType == "Gauss") || (funcType == "Gaussian") )
        continue;

      else if (funcType == "File")
        tmpVal = evalFileOptGen(x,y,z,t,funcParamList);

      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid OptGen function type!"
          << "Must be Gauss or File." << std::endl );

      genVal += tmpVal;

    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop




  return genVal;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalGaussianOptGen()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double OptGen_Function<EvalT, Traits>::evalGaussianOptGen(const double& x, const double& y,
  const double& z, const double& t, const gaussianOptGenParams& gp)
{
  using std::string;
  using Teuchos::ParameterList;

  const double maxVal = gp.maxVal;
  const double minVal = gp.minVal;

  // x direction
  const string x_dir = gp.x_dir;
  const double x_loc = gp.x_loc;
  const double x_width = gp.x_width;
  const double x_min = gp.x_min;
  const double x_max = gp.x_max;
  const bool x_checkAxis = gp.x_checkAxis;

  // y direction
  const string y_dir = gp.y_dir;
  const double y_loc = gp.y_loc;
  const double y_width = gp.y_width;
  const double y_min = gp.y_min;
  const double y_max = gp.y_max;
  const bool y_checkAxis = gp.y_checkAxis;

  // z direction
  const string z_dir = gp.z_dir;
  const double z_loc = gp.z_loc;
  const double z_width = gp.z_width;
  const double z_min = gp.z_min;
  const double z_max = gp.z_max;
  const bool z_checkAxis = gp.z_checkAxis;

  // time direction
  const string t_dir = gp.t_dir;
  const double t_loc = gp.t_loc;
  const double t_width = gp.t_width;
  const double t_min = gp.t_min;
  const double t_max = gp.t_max;
  const bool t_checkAxis = gp.t_checkAxis;

  bool found = false;
  double xGaussVal = 1.0, yGaussVal = 1.0, zGaussVal = 1.0, tGaussVal = 1.0;

  xGaussVal = evalSingleGaussian("X", found, x, minVal, maxVal, x_min, x_max, x_loc, x_width, x_checkAxis, x_dir);
  if (num_dim == 2)
    yGaussVal = evalSingleGaussian("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
  if (num_dim == 3)
  {
    yGaussVal = evalSingleGaussian("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
    zGaussVal = evalSingleGaussian("Z", found, z, minVal, maxVal, z_min, z_max, z_loc, z_width, z_checkAxis, z_dir);
  }
  tGaussVal = evalSingleGaussian("Time", found, t, minVal, maxVal, t_min, t_max, t_loc, t_width, t_checkAxis, t_dir);

  // throw exception if NO Gaussian profile is specified
  if (!found)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Gaussian is specified "
      << "for Optical Generation / Function Type of Gauss/Gaussian! At least one Gaussian along "
      << "x, y, z, or t must be specified! ");

  // compute the opt. gen. rate in [#/(cm^3.s)]
  double genValue = maxVal*xGaussVal*yGaussVal*zGaussVal*tGaussVal;
  return genValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleGaussian()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double OptGen_Function<EvalT, Traits>::evalSingleGaussian
  (const std::string& axis, bool& found, const double& coord,
   const double& minVal, const double& maxVal, const double& min, const double& max,
   const double& loc, const double& width, const bool& checkAxis, const std::string& dir)
{
  using std::string;
  using Teuchos::ParameterList;

  // If Gaussian is NOT specified for a certain axis (say X), returns 1.0
  double gaussVal = 1.0;

  // Gaussian is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if a Gaussian is set along an axis, then found = true

    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
    {
      if (dir == "Both")
        gaussVal = exp(-log(maxVal/minVal) * pow((coord-loc)/width, 2.0) );
      else if (dir == "Positive")
      {
        if (coord >= loc)
          gaussVal = exp(-log(maxVal/minVal) * pow((coord-loc)/width, 2.0));
        else
          gaussVal = 1.0;
      }
      else if (dir == "Negative")
      {
        if (coord <= loc)
          gaussVal = exp(-log(maxVal/minVal) * pow((coord-loc)/width, 2.0));
        else
          gaussVal = 1.0;
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Direction must be either Both, Positive, or Negative !");
    }

    // If Gaussian is specified for a certain axis but coord is outside [min, max], returns 0.0
    else
      gaussVal = 0.0; // 0 outside the [min, max) range

  }

  return gaussVal;
}




//**********************************************************************
// Exports OptGen profile from an external file and maps values
// to each (t) point if (x,y,z) are within the user-defined range
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  evalFileOptGen()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double OptGen_Function<EvalT, Traits>::evalFileOptGen(const double& x, const double& y, const double& z, const double& t, const Teuchos::ParameterList& plist)

{
  using Teuchos::ParameterList;
  using std::string;
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::vector;


  optgen_struct r;

  double OptGenVal=0.0;
  double t0, t1, v0, v1;
  string method;

  double xmin = -1e100, ymin = -1e100, zmin = -1e100;
  double xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if (plist.isParameter("X Min"))  xmin = plist.get<double>("X Min");
  if (plist.isParameter("X Max"))  xmax = plist.get<double>("X Max");
  if (plist.isParameter("Y Min"))  ymin = plist.get<double>("Y Min");
  if (plist.isParameter("Y Max"))  ymax = plist.get<double>("Y Max");
  if (plist.isParameter("Z Min"))  zmin = plist.get<double>("Z Min");
  if (plist.isParameter("Z Max"))  zmax = plist.get<double>("Zmax");

//  cout << "t=" << t << "GV.size()=" << GV.size() << " values:" << GV[0].t << "," << GV[0].g << ";" << GV[4].t << "," << GV[4].g << endl;


  if (x>=xmin && x<=xmax && y>=ymin && y<=ymax && z>=zmin && z<=zmax)
  {
    r.t=t;

    size_t closest_point = 0;
    double closest_distance = GV[0].distance(r);

    for (size_t i=1; i<=GV.size()-1; i++)
    {
     double distance = GV[i].distance(r);
     if (distance < closest_distance)
     {
      closest_distance = distance;
      closest_point = i;
     }
    }
//  if (t>2e-5) cout << "t=" << t << ", x=" << x <<", y=" << y << ", GV.t=" << GV[closest_point].t << ", GV.g=" << GV[closest_point].g << endl;


  // generation in Charon is POSITIVE
   OptGenVal = fabs(GV[closest_point].g);

   if (plist.isParameter("Interpolation"))
   {
       method = plist.get<string>("Interpolation");
       if ( method == "Linear")
         {
           if (closest_point>1 && closest_point<GV.size())
             {
               if (t>GV[closest_point].t)
                 {
                  t0=GV[closest_point].t;
                  t1=GV[closest_point+1].t;
                  v0=fabs(GV[closest_point].g);
                  v1=fabs(GV[closest_point+1].g);
                  OptGenVal = v0+(v1-v0)*(t-t0)/(t1-t0);
                 }
               else if (t < GV[closest_point].t)
                 {
                  t0=GV[closest_point-1].t;
                  t1=GV[closest_point].t;
                  v0=fabs(GV[closest_point-1].g);
                  v1=fabs(GV[closest_point].g);
                  OptGenVal = v0+(v1-v0)*(t-t0)/(t1-t0);
                 }

//               if (t>2e-5) cout << "t=" << t << " t0=" << t0 << " t1=" << t1 << " v0=" << v0 << " v1=" << v1 << " V=" << OptGenVal << endl;

             }
         }
       else if (method == "Spline")
         {

           OptGenVal = fabs(ogSpline->evaluateSpline(t));

         }

       else
           TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid Interpolation method in OptGen File!"
            << "If present it must be Linear or Spline." << std::endl );
   }





  }  //if x>=xmin...

  return OptGenVal;
}





} // namespace charon

#endif
