
#ifndef CHARON_IC_GAUSS_IMPL_HPP
#define CHARON_IC_GAUSS_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"

/*
This evaluator implements a spatial Gaussian profile for initial condition.
N = Nmax * a(x) * b(y) * c(z)
a(x) = exp[-ln(Nmax/Nmin) * ((x-xgmin) / xwidth)^2] for x < xgmin
a(x) = 1                                            for xgmin <= x <= xgmax
a(x) = exp[-ln(Nmax/Nmin) * ((x-xgmax) / xwidth)^2] for x > xgmax
a(x) = 0                                            for x > xmax and x < xmin
b(y) and c(z) follow similar expressions.

Nmax = Max Value, Nmin = Min Value, xwidth = X Width, xgmin = X Gauss Min,
xgmax = X Gauss Max, xmin = X Min, xmax = X Max.

When a Gauss is not given for an axis, e.g., z axis, c(z) = 1 for zmin <= z <= zmax,
and c(z) = 0 for outside [zmin zmax] range.

An example of using the model is given below, where Max Value and Min Value are
scaled values, e.g., in units of C0 (concentration scaling), while all coordinate
related values are in unit of [um].
            <ParameterList name="ION_DENSITY">
                <Parameter name="Value" type="string" value="Gauss"/>
                <Parameter name="Max Value" type="double" value="1.0"/>
                <Parameter name="Min Value" type="double" value="0.5"/>
                <Parameter name="X Width" type="double" value="0.004"/>
                <Parameter name="X Gauss Min" type="double" value="-0.0001"/>
                <Parameter name="X Gauss Max" type="double" value="0.0"/>
                <Parameter name="X Min" type="double" value="-0.01"/>
                <Parameter name="X Max" type="double" value="0.01"/>
                <Parameter name="Y Width" type="double" value="0.004"/>
                <Parameter name="Y Gauss Min" type="double" value="-0.0001"/>
                <Parameter name="Y Gauss Max" type="double" value="0.0"/>
                <Parameter name="Y Min" type="double" value="-0.01"/>
                <Parameter name="Y Max" type="double" value="0.01"/>
            </ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IC_Gauss<EvalT, Traits>::
IC_Gauss(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  num_basis = data_layout->dimension(1);
  basis_name = basis->name();
  num_dim = basis->dimension();

  dof_name = p.get<string>("DOF Name");

  // Obtain the parameterlist
  const ParameterList& paramList = p.sublist("Gauss ParameterList");

  // Initialize the parameters
  initialize(paramList);

  // Evaluated field
  carrier_density = MDField<ScalarT,Cell,BASIS>(dof_name, data_layout);
  this->addEvaluatedField(carrier_density);

  std::string name = "IC_Gauss";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Gauss<EvalT, Traits>::
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
IC_Gauss<EvalT, Traits>::
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

      double icValue = evaluateGaussIC(x,y,z);
      carrier_density(cell,basis) = icValue;  // scaled
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void IC_Gauss<EvalT, Traits>::initialize(const Teuchos::ParameterList& plist)
{
  // Set up parameters for the Gaussian profile
  // Max Value and Min Value must be given
  maxValue = plist.get<double>("Max Value");
  minValue = plist.get<double>("Min Value");

  // initialize all coordinates to zero first
  xWidth = 0.; xGaussMin = 0.; xGaussMax = 0.; xMin = 0.; xMax = 0.;
  yWidth = 0.; yGaussMin = 0.; yGaussMax = 0.; yMin = 0.; yMax = 0.;
  zWidth = 0.; zGaussMin = 0.; zGaussMax = 0.; zMin = 0.; zMax = 0.;

  // initialize booleans to false first
  xAxis = false; yAxis = false; zAxis = false;

  // obtain actual values from input xml
  testcoord("X", plist, xWidth, xGaussMin, xGaussMax, xMin, xMax, xAxis);
  if (num_dim == 2)
    testcoord("Y", plist, yWidth, yGaussMin, yGaussMax, yMin, yMax, yAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, yWidth, yGaussMin, yGaussMax, yMin, yMax, yAxis);
    testcoord("Z", plist, zWidth, zGaussMin, zGaussMax, zMin, zMax, zAxis);
  }

  if ((!xAxis) && (!yAxis) && (!zAxis))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Gauss is specified! "
     << "At least one Gauss along x, y, or z must be specified!");

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  testcoord()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void IC_Gauss<EvalT, Traits>::testcoord(const std::string& axis, const Teuchos::ParameterList& plist,
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
          << "Error ! " << axis << " Gauss Min must be smaller than or equal to " << axis << " Gauss Max !");

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

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateGaussIC()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Gauss<EvalT, Traits>::evaluateGaussIC(const double& x, const double& y, const double& z)
{
  double xGaussVal = 1.0, yGaussVal = 1.0, zGaussVal = 1.0;

  xGaussVal = evalSingleGauss(x, xWidth, xGaussMin, xGaussMax, xMin, xMax, xAxis);
  if (num_dim == 2)
    yGaussVal = evalSingleGauss(y, yWidth, yGaussMin, yGaussMax, yMin, yMax, yAxis);
  if (num_dim == 3)
  {
    yGaussVal = evalSingleGauss(y, yWidth, yGaussMin, yGaussMax, yMin, yMax, yAxis);
    zGaussVal = evalSingleGauss(z, zWidth, zGaussMin, zGaussMax, zMin, zMax, zAxis);
  }

  double icValue = maxValue*xGaussVal*yGaussVal*zGaussVal;

  return icValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleGauss()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double IC_Gauss<EvalT, Traits>::evalSingleGauss (const double& coord, const double& width,
                          const double& gaussMin, const double& gaussMax,
                          const double& min, const double& max, const bool& checkAxis)
{
  // If Gauss is NOT specified for a certain axis (say X), returns 1.0
  double gaussVal = 1.0;

  // Set gaussVal to 0 when coord is outside [min max] range
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


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
IC_Gauss<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("DOF Name", "?");

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->sublist("Gauss ParameterList", false, "");
  p->sublist("Gauss ParameterList").set<std::string>("Value", "Gauss", "Gaussian or Gauss profile for initial condition");
  p->sublist("Gauss ParameterList").set<double>("Max Value", 0., "in unit of C0");
  p->sublist("Gauss ParameterList").set<double>("Min Value", 0., "in unit of C0");
  p->sublist("Gauss ParameterList").set<double>("X Width", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("X Gauss Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("X Gauss Max", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("X Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("X Max", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Y Width", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Y Gauss Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Y Gauss Max", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Y Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Y Max", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Z Width", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Z Gauss Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Z Gauss Max", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Z Min", 0., "[um]");
  p->sublist("Gauss ParameterList").set<double>("Z Max", 0., "[um]");

  return p;
}

}

#endif

