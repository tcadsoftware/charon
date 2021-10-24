
#ifndef CHARON_IC_FUNCTION_IMPL_HPP
#define CHARON_IC_FUNCTION_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"

/*
This evaluator implements a combination of analytic spatial profiles for initial
condition. It currently supports Uniform and Gaussian/Gauss profiles.

The specification of this model is similar to the Doping_Function, except it uses
X Min, X Max, Y Min, Y Max, Z Min, and Z Max, instead of Xmin, Xmax, Ymin, Ymax,
Zmin, and Zmax that are used in Doping_Function.

To be tested.
*/


void charon::uniformICParams::parseUniform (const Teuchos::ParameterList& plist)
{
  value = plist.get<double>("IC Value");

  // default valuues
  xmin = -1e100, ymin = -1e100, zmin = -1e100;
  xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if (plist.isParameter("X Min"))  xmin = plist.get<double>("X Min");
  if (plist.isParameter("X Max"))  xmax = plist.get<double>("X Max");
  if (plist.isParameter("Y Min"))  ymin = plist.get<double>("Y Min");
  if (plist.isParameter("Y Max"))  ymax = plist.get<double>("Y Max");
  if (plist.isParameter("Z Min"))  zmin = plist.get<double>("Z Min");
  if (plist.isParameter("Z Max"))  zmax = plist.get<double>("Z Max");
}

void charon::gaussianICParams::parseGaussian (const Teuchos::ParameterList& plist,
                                              const int num_dim)
{
  // Max Value and Min Value must be given
  maxValue = plist.get<double>("Max Value");
  minValue = plist.get<double>("Min Value");

  // initialize all coordinates to zero first
  xWidth = 0.; xGaussMin = 0.; xGaussMax = 0.; xMin = 0.; xMax = 0.;
  yWidth = 0.; yGaussMin = 0.; yGaussMax = 0.; yMin = 0.; yMax = 0.;
  zWidth = 0.; zGaussMin = 0.; zGaussMax = 0.; zMin = 0.; zMax = 0.;

  // initialize booleans to false first
  xCheckAxis = false; yCheckAxis = false; zCheckAxis = false;

  // obtain actual values from input xml
  testcoord("X", plist, xWidth, xGaussMin, xGaussMax, xMin, xMax, xCheckAxis);
  if (num_dim == 2)
    testcoord("Y", plist, yWidth, yGaussMin, yGaussMax, yMin, yMax, yCheckAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, yWidth, yGaussMin, yGaussMax, yMin, yMax, yCheckAxis);
    testcoord("Z", plist, zWidth, zGaussMin, zGaussMax, zMin, zMax, zCheckAxis);
  }

  if ((!xCheckAxis) && (!yCheckAxis) && (!zCheckAxis))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Gauss is specified! "
     << "At least one Gauss along x, y, or z must be specified!");

  return;
}

void charon::gaussianICParams::testcoord (const std::string& axis, const Teuchos::ParameterList& plist,
                                          double& width, double& gaussMin, double& gaussMax,
                                          double& min, double& max, bool& checkAxis)
{
  if (plist.isParameter(axis+" Width"))
  {
    width = plist.get<double>(axis+" Width");
    gaussMin = plist.get<double>(axis+" Gauss Min");
    gaussMax = plist.get<double>(axis+" Gauss Max");
    checkAxis = true;

    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+" Min")) min = plist.get<double>(axis+" Min");
    if (plist.isParameter(axis+" Max")) max = plist.get<double>(axis+" Max");

    if (min >= max)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Min must be smaller than " << axis << " Max !");

    if (gaussMin > gaussMax)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Gauss Min must be smaller or equal to " << axis << " Gauss Max !");

    if (gaussMin < min)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Gauss Min must be greater or equal to " << axis << " Min !");

    if (gaussMax > max)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Gauss Max must be smaller or equal to " << axis << " Max !");
  }
  else
  {
    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+" Min")) min = plist.get<double>(axis+" Min");
    if (plist.isParameter(axis+" Max")) max = plist.get<double>(axis+" Max");
  }

}


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IC_Function<EvalT, Traits>::
IC_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  // RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  // p.validateParameters(*valid_params);

  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  num_basis = data_layout->dimension(1);
  basis_name = basis->name();
  num_dim = basis->dimension();

  dof_name = p.get<string>("DOF Name");

  // Obtain the parameterlist
  const ParameterList& paramList = p.sublist("Function ParameterList");

  // Evaluated field
  carrier_density = MDField<ScalarT,Cell,BASIS>(dof_name, data_layout);
  this->addEvaluatedField(carrier_density);

  std::string name = "IC_Function";
  this->setName(name);

  for (ParameterList::ConstIterator model_it = paramList.begin();
       model_it !=paramList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      if (funcType == "Uniform")
      {
        uniformICParams udp_;
        udp_.parseUniform(funcParamList);
        udp_vec.push_back(udp_);
      }
      if ( (funcType == "Gauss") || (funcType == "Gaussian") )
      {
        gaussianICParams gdp_;
        gdp_.parseGaussian(funcParamList,num_dim);
        gdp_vec.push_back(gdp_);
      }

    }  // end of if (key.compare(0, 8, "Function") == 0)
  }  // end of for loop

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Function<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
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

      double icValue = evaluateIC(x,y,z);
      carrier_density(cell,basis) = icValue;  // scaled
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateIC()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Function<EvalT, Traits>::evaluateIC(const double& x, const double& y, const double& z)
{
  double tempValue = 0.0;
  double icValue = 0.0;

  for (std::size_t i = 0; i < udp_vec.size(); ++i)
  {
    tempValue = evalUniformIC(x,y,z,udp_vec[i]);
    icValue += tempValue;
  }
  for (std::size_t i = 0; i < gdp_vec.size(); ++i)
  {
    tempValue = evalGaussianIC(x,y,z,gdp_vec[i]);
    icValue += tempValue;
  }
  return icValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalUniformIC()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Function<EvalT, Traits>::evalUniformIC
  (const double& x, const double& y, const double& z, const uniformICParams& udp)
{
  const double value = udp.value;
  const double xmin = udp.xmin;
  const double ymin = udp.ymin;
  const double zmin = udp.zmin;
  const double xmax = udp.xmax;
  const double ymax = udp.ymax;
  const double zmax = udp.zmax;

  double icValue = 0.0;
  if ( (x >= xmin) && (x <= xmax) && (y >= ymin) && (y <= ymax) && (z >= zmin) && (z <= zmax) )
    icValue = value;

  return icValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalGaussianIC()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Function<EvalT, Traits>::evalGaussianIC
  (const double& x, const double& y, const double& z, const gaussianICParams& gdp)
{
  const double minValue = gdp.minValue;
  const double maxValue = gdp.maxValue;

  // x direction
  const double xWidth = gdp.xWidth;
  const double xGaussMin = gdp.xGaussMin;
  const double xGaussMax = gdp.xGaussMax;
  const double xMin = gdp.xMin;
  const double xMax = gdp.xMax;
  const bool xCheckAxis  = gdp.xCheckAxis;

  // y direction
  const double yWidth = gdp.yWidth;
  const double yGaussMin = gdp.yGaussMin;
  const double yGaussMax = gdp.yGaussMax;
  const double yMin = gdp.yMin;
  const double yMax = gdp.yMax;
  const bool yCheckAxis  = gdp.yCheckAxis;

  // z direction
  const double zWidth = gdp.zWidth;
  const double zGaussMin = gdp.zGaussMin;
  const double zGaussMax = gdp.zGaussMax;
  const double zMin = gdp.zMin;
  const double zMax = gdp.zMax;
  const bool zCheckAxis  = gdp.zCheckAxis;

  double xGaussVal = 1.0, yGaussVal = 1.0, zGaussVal = 1.0;

  xGaussVal = evalSingleGaussian(x, minValue, maxValue, xWidth, xGaussMin, xGaussMax, xMin, xMax, xCheckAxis);
  if (num_dim == 2)
    yGaussVal = evalSingleGaussian(y, minValue, maxValue, yWidth, yGaussMin, yGaussMax, yMin, yMax, yCheckAxis);
  if (num_dim == 3)
  {
    yGaussVal = evalSingleGaussian(y, minValue, maxValue, yWidth, yGaussMin, yGaussMax, yMin, yMax, yCheckAxis);
    zGaussVal = evalSingleGaussian(z, minValue, maxValue, zWidth, zGaussMin, zGaussMax, zMin, zMax, zCheckAxis);
  }

  double icValue = maxValue*xGaussVal*yGaussVal*zGaussVal;

  return icValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleGaussian()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Function<EvalT, Traits>::evalSingleGaussian
  (const double& coord, const double& minValue, const double& maxValue,
   const double& width, const double& gaussMin, const double& gaussMax,
   const double& min, const double& max, const bool& checkAxis)
{
  // If Gauss is NOT specified for a certain axis (say X), returns 1.0
  double gaussVal = 1.0;

  // set gaussVal to 0 when coord is outside [min max] range
  if ( (coord < min) || (coord > max) )  gaussVal = 0.0;

  // Gauss is specified along a given axis
  if (checkAxis)
  {
    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
    {
      // Gauss decay for coord < gaussMin
      if (coord < gaussMin)
        gaussVal = std::exp(-std::log(maxValue/minValue) * std::pow((coord-gaussMin)/width, 2.0));

      // Gauss decay
      else if (coord > gaussMax)
        gaussVal = std::exp(-std::log(maxValue/minValue) * std::pow((coord-gaussMax)/width, 2.0));

      // uniform within [gaussMin gaussMax] range
      else
        gaussVal = 1.0;
    }

    // If Gauss is specified for a certain axis but coord is outside [min, max], returns 0.0
    else
      gaussVal = 0.0; // 0 outside the [min, max] range
  }

  return gaussVal;
}

}

#endif

