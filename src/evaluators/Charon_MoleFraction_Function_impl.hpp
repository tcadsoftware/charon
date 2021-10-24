
#ifndef CHARON_MOLEFRACTION_FUNCTION_IMPL_HPP
#define CHARON_MOLEFRACTION_FUNCTION_IMPL_HPP

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

/*
Mole Fraction specification is similar to that of Doping, except that Mole Fraction
currently supports only Linear and Uniform profiles.
An example of specifying Mole Fraction is given below:

<ParameterList name="Mole Fraction">
    <Parameter name="Value" type="string" value="Function"/>
    <ParameterList name="Function 1">
        <Parameter name="Function Type" type="string" value="Linear"/>
        <Parameter name="Mole Start Value" type="double" value="0.3"/>
        <Parameter name="Mole End Value" type="double" value="0.0"/>
        <Parameter name="Xmin" type="double" value="0.3"/>
        <Parameter name="Xmax" type="double" value="0.7"/>
    </ParameterList>
    <ParameterList name="Function 2">
        <Parameter name="Function Type" type="string" value="Uniform"/>
        <Parameter name="Mole Value" type="double" value="0.3"/>
        <Parameter name="Xmin" type="double" value="-4.0"/>
        <Parameter name="Xmax" type="double" value="-3.5"/>
        <Parameter name="Ymin" type="double" value="-0.2"/>
        <Parameter name="Ymax" type="double" value="0.0"/>
    </ParameterList>
</ParameterList>
*/

namespace charon {

void uniformMoleFracParams::parseUniform (const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;

  value = plist.get<double>("Mole Value");
  xmin = -1e100, ymin = -1e100, zmin = -1e100;
  xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if ((value < 0.) || (value > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole Value must be within [0 1] !");

  if (plist.isParameter("Xmin"))  xmin = plist.get<double>("Xmin");
  if (plist.isParameter("Xmax"))  xmax = plist.get<double>("Xmax");
  if (plist.isParameter("Ymin"))  ymin = plist.get<double>("Ymin");
  if (plist.isParameter("Ymax"))  ymax = plist.get<double>("Ymax");
  if (plist.isParameter("Zmin"))  zmin = plist.get<double>("Zmin");
  if (plist.isParameter("Zmax"))  zmax = plist.get<double>("Zmax");
}

void linearMoleFracParams::parseLinear (const Teuchos::ParameterList& plist, const int num_dim)
{

  using std::string;
  using Teuchos::ParameterList;

  startVal = plist.get<double>("Mole Start Value");
  endVal = plist.get<double>("Mole End Value");

  if ((startVal < 0.) || (startVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole Start Value must be within [0 1] !");

  if ((endVal < 0.) || (endVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole End Value must be within [0 1] !");

  testcoord("X", plist, x_min, x_max, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
    testcoord("Z", plist, z_min, z_max, z_checkAxis);
  }
}

void linearMoleFracParams::testcoord (const std::string axis, const Teuchos::ParameterList& plist,
     double& min, double& max, bool& checkAxis)
{
  checkAxis = false;
  if (plist.isParameter(axis+"min") || plist.isParameter(axis+"max"))
  {
    checkAxis = true;
    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
MoleFraction_Function<EvalT, Traits>::
MoleFraction_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

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
  num_basis = basis_scalar->dimension(1);

  // mole fraction parameterlist
  moleFracParamList = p.sublist("Mole Fraction ParameterList");

  // evaluated fields
  molefrac = MDField<ScalarT,Cell,IP>(n.field.mole_frac,scalar);
  molefrac_basis = MDField<ScalarT,Cell,BASIS>(n.field.mole_frac,basis_scalar);

  this->addEvaluatedField(molefrac);
  this->addEvaluatedField(molefrac_basis);

  std::string name = "MoleFraction_Function";
  this->setName(name);

  for (ParameterList::ConstIterator model_it = moleFracParamList.begin();
       model_it != moleFracParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      if (funcType == "Uniform")
      {
        uniformMoleFracParams umfp_;
        umfp_.parseUniform(funcParamList);
        umfp_vec.push_back(umfp_);
      }
      if (funcType == "Linear")
      {
        linearMoleFracParams lmfp_;
        lmfp_.parseLinear(funcParamList,num_dim);
        lmfp_vec.push_back(lmfp_);
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
MoleFraction_Function<EvalT, Traits>::
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
MoleFraction_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // x mole fraction at IPs
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

      // evaluate the x mole fraction
      double value = evaluateMoleFraction(x,y,z);
      molefrac(cell,ip) = value;
    }

    // x mole fraction at basis points
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

      double value = evaluateMoleFraction(x,y,z);
      molefrac_basis(cell,basis) = value;
    }

  } // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateMoleFraction()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evaluateMoleFraction
  (const double& x, const double& y, const double& z)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  double tempVal = 0.0;

  for (std::size_t i = 0; i < umfp_vec.size(); ++i)
  {
    tempVal = evalUniformMoleFrac(x,y,z,umfp_vec[i]);
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < lmfp_vec.size(); ++i)
  {
    tempVal = evalLinearMoleFrac(x,y,z,lmfp_vec[i]);
    mfValue += tempVal;
  }

  return mfValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalUniformMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalUniformMoleFrac
  (const double& x, const double& y, const double& z, const uniformMoleFracParams& umfp)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  const double val = umfp.value;
  const double xmin = umfp.xmin;
  const double ymin = umfp.ymin;
  const double zmin = umfp.zmin;
  const double xmax = umfp.xmax;
  const double ymax = umfp.ymax;
  const double zmax = umfp.zmax;

  if ( (x >= xmin) && (x <= xmax) && (y >= ymin) && (y <= ymax) && (z >= zmin) && (z <= zmax) )
    mfValue = val;

  // return 0 if (x,y,z) is outside the box region

  return mfValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalLinearMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalLinearMoleFrac
  (const double& x, const double& y, const double& z, const linearMoleFracParams& lmfp)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  const double startVal = lmfp.startVal;
  const double endVal = lmfp.endVal;

  // x direction
  const double x_min = lmfp.x_min;
  const double x_max = lmfp.x_max;
  const bool x_checkAxis = lmfp.x_checkAxis;

  // y direction
  const double y_min = lmfp.y_min;
  const double y_max = lmfp.y_max;
  const bool y_checkAxis = lmfp.y_checkAxis;

  // z direction
  const double z_min = lmfp.z_min;
  const double z_max = lmfp.z_max;
  const bool z_checkAxis = lmfp.z_checkAxis;

  bool found = false;
  double xLinearVal = 1.0, yLinearVal = 1.0, zLinearVal = 1.0;

  xLinearVal = evalSingleLinear("X", found, x, x_min, x_max, x_checkAxis);
  if (num_dim == 2)
    yLinearVal = evalSingleLinear("Y", found, y, y_min, y_max, y_checkAxis);
  if (num_dim == 3)
  {
    yLinearVal = evalSingleLinear("Y", found, y, y_min, y_max, y_checkAxis);
    zLinearVal = evalSingleLinear("Z", found, z, z_min, z_max, z_checkAxis);
  }

  // throw exception if NO Linear profile is specified
  if (!found)
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No linear function is specified "
     << "for Function Type = Linear! At least one linear function along "
     << "x, y, or z must be specified! ");

  if ( (xLinearVal >= 0.0) && (yLinearVal >= 0.0) && (zLinearVal >= 0.0) )
    mfValue = startVal + (endVal-startVal)*xLinearVal*yLinearVal*zLinearVal;

  // return 0 if (x,y,z) is outside the box region

  return mfValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleLinear()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalSingleLinear
  (const std::string /* axis */, bool& found, const double& coord, const double& min,
   const double& max, const bool& checkAxis)
{
  // If Linear is NOT specified for a certain axis (say X), returns 1.0
  double LinearVal = 1.0;

  // Linear is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if Linear is set along an axis, then found = true

    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
     LinearVal = (coord-min)/(max-min);
    else
     LinearVal = -1.0; // -1 is a flag indicating outside the [min,max] range
  }

  return LinearVal;
}


} // namespace charon

#endif
