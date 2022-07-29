
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
#include "Charon_Material_Properties.hpp"

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

For binary and ternary materials:

<ParameterList name="Mole Fraction">
    <Parameter name="Value" type="string" value="Function"/>
    <ParameterList name="Function 1">
        <Parameter name="Function Type" type="string" value="Linear"/>
        <!--
        <Parameter name="Mole Start Value" type="double" value="0.3"/>
        <Parameter name="Mole End Value" type="double" value="0.0"/>
        -->
        <Parameter name="xMoleFrac Start Value" type="double" value="0.3"/>
        <Parameter name="xMoleFrac End Value" type="double" value="0.0"/>
        <Parameter name="Xmin" type="double" value="0.3"/>
        <Parameter name="Xmax" type="double" value="0.7"/>
    </ParameterList>
    <ParameterList name="Function 2">
        <Parameter name="Function Type" type="string" value="Uniform"/>
        <!-- <Parameter name="Mole Value" type="double" value="0.3"/> -->
        <Parameter name="xMoleFrac Value" type="double" value="0.3"/>
        <Parameter name="Xmin" type="double" value="-4.0"/>
        <Parameter name="Xmax" type="double" value="-3.5"/>
        <Parameter name="Ymin" type="double" value="-0.2"/>
        <Parameter name="Ymax" type="double" value="0.0"/>
    </ParameterList>
</ParameterList>

For quaternary materials:

<ParameterList name="Mole Fraction">
    <Parameter name="Value" type="string" value="Function"/>
    <ParameterList name="Function 1">
        <Parameter name="Function Type" type="string" value="Linear"/>
        <!--
        <Parameter name="Mole Start Value" type="double" value="0.3"/>
        <Parameter name="Mole End Value" type="double" value="0.0"/>
        -->
        <Parameter name="xMoleFrac Start Value" type="double" value="0.3"/>
        <Parameter name="xMoleFrac End Value" type="double" value="0.0"/>
        <Parameter name="yMoleFrac Start Value" type="double" value="0.2"/>
        <Parameter name="yMoleFrac End Value" type="double" value="0.0"/>
        <Parameter name="Xmin" type="double" value="0.3"/>
        <Parameter name="Xmax" type="double" value="0.7"/>
    </ParameterList>
    <ParameterList name="Function 2">
        <Parameter name="Function Type" type="string" value="Uniform"/>
        <!-- <Parameter name="Mole Value" type="double" value="0.3"/> -->
        <Parameter name="xMoleFrac Value" type="double" value="0.3"/>
        <Parameter name="yMoleFrac Value" type="double" value="0.2"/>
        <Parameter name="Xmin" type="double" value="-4.0"/>
        <Parameter name="Xmax" type="double" value="-3.5"/>
        <Parameter name="Ymin" type="double" value="-0.2"/>
        <Parameter name="Ymax" type="double" value="0.0"/>
    </ParameterList>
</ParameterList>
*/

namespace charon {

void uniformMoleFracParams::parseUniform(const Teuchos::ParameterList& plist, 
                                         const std::string& matArity) {
  using Teuchos::ParameterList;
  using std::string;

  mat_arity = matArity;

  // value = plist.get<double>("Mole Value");
  if (plist.isParameter("Mole Value"))
    value = plist.get<double>("Mole Value");
  else if (plist.isParameter("xMoleFrac Value"))
    value = plist.get<double>("xMoleFrac Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Value' or 'xMoleFrac Value' must be specified!");
  if (mat_arity == "Quaternary") {
    value1 = plist.get<double>("yMoleFrac Value");
  }
  xmin = -1e100, ymin = -1e100, zmin = -1e100;
  xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if ((value < 0.) || (value > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole Value must be within [0 1] !");

  if (mat_arity == "Quaternary") {
    if ((value1 < 0.) || (value1 > 1.))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	<< "Error ! yMoleFrac Value must be within [0 1] !");  
  }

  if (plist.isParameter("Xmin"))  xmin = plist.get<double>("Xmin");
  if (plist.isParameter("Xmax"))  xmax = plist.get<double>("Xmax");
  if (plist.isParameter("Ymin"))  ymin = plist.get<double>("Ymin");
  if (plist.isParameter("Ymax"))  ymax = plist.get<double>("Ymax");
  if (plist.isParameter("Zmin"))  zmin = plist.get<double>("Zmin");
  if (plist.isParameter("Zmax"))  zmax = plist.get<double>("Zmax");
}


void linearMoleFracParams::parseLinear(const Teuchos::ParameterList& plist, 
                                       const int num_dim, const std::string& matArity) {
  using std::string;
  using Teuchos::ParameterList;

  mat_arity = matArity;

  /*
  startVal = plist.get<double>("Mole Start Value");
  endVal = plist.get<double>("Mole End Value");
  */
  if (plist.isParameter("Mole Start Value"))
    startVal = plist.get<double>("Mole Start Value");
  else if (plist.isParameter("xMoleFrac Start Value"))
    startVal = plist.get<double>("xMoleFrac Start Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Start Value' or 'xMoleFrac Start Value' must be specified!");
  if (plist.isParameter("Mole End Value"))
    endVal = plist.get<double>("Mole End Value");
  else if (plist.isParameter("xMoleFrac End Value"))
    endVal = plist.get<double>("xMoleFrac End Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole End Value' or 'xMoleFrac End Value' must be specified!");
  if (mat_arity == "Quaternary") {
    if (!plist.isParameter("yMoleFrac Start Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac Start Value' must be specified!");
    if (!plist.isParameter("yMoleFrac End Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac End Value' must be specified!");
    startVal1 = plist.get<double>("yMoleFrac Start Value");
    endVal1 = plist.get<double>("yMoleFrac End Value");
  }

  if ((startVal < 0.) || (startVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole Start Value must be within [0 1] !");

  if ((endVal < 0.) || (endVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Mole End Value must be within [0 1] !");

  if (mat_arity == "Quaternary") {
    if ((startVal1 < 0.) || (startVal1 > 1.))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error ! yMoleFrac Start Value must be within [0 1] !");
    if ((endVal1 < 0.) || (endVal1 > 1.))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error ! yMoleFrac End Value must be within [0 1] !");
  }

  testcoord("X", plist, x_min, x_max, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
    testcoord("Z", plist, z_min, z_max, z_checkAxis);
  }
}

void linearMoleFracParams::testcoord(const std::string axis, const Teuchos::ParameterList& plist,
     double& min, double& max, bool& checkAxis) {
  checkAxis = false;
  if (plist.isParameter(axis+"min") || plist.isParameter(axis+"max"))
  {
    checkAxis = true;
    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }
}


void gaussMoleFracParams::parseGaussian(const Teuchos::ParameterList& plist, 
                                       const int num_dim, const std::string& matArity) {
  using std::string;
  using Teuchos::ParameterList;

  mat_arity = matArity;

  if (plist.isParameter("Mole Max Value"))
    maxVal = plist.get<double>("Mole Max Value");
  else if (plist.isParameter("xMoleFrac Max Value"))
    maxVal = plist.get<double>("xMoleFrac Max Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Max Value' or 'xMoleFrac Max Value' must be specified!");
  if (plist.isParameter("Mole Min Value"))
    minVal = plist.get<double>("Mole Min Value");
  else if (plist.isParameter("xMoleFrac Min Value"))
    minVal = plist.get<double>("xMoleFrac Min Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Min Value' or 'xMoleFrac Min Value' must be specified!");
  if ((minVal < 0.) || (maxVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gaussian xMoleFrac Values must be within [0 1] !");
  if (mat_arity == "Quaternary") {
    if (!plist.isParameter("yMoleFrac Max Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac Max Value' must be specified!");
    if (!plist.isParameter("yMoleFrac Min Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac Min Value' must be specified!");
    maxVal1 = plist.get<double>("yMoleFrac Max Value");
    minVal1 = plist.get<double>("yMoleFrac Min Value");
    if ((minVal1 < 0.) || (maxVal1 > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gaussian yMoleFrac Values must be within [0 1] !");
  }

  testcoord("X", plist, x_dir, x_min, x_max, x_loc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_loc, z_width, z_checkAxis);
  }
}


void gaussMoleFracParams::testcoord(const std::string axis, 
				     const Teuchos::ParameterList& plist,
				     std::string& dir, double& min,
				     double& max, double& loc, double& width,
				     bool& checkAxis) {
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
  if (plist.isParameter(axis+" Peak Location") && plist.isParameter(axis+" Width")) {
    loc = plist.get<double>(axis+" Peak Location");
    width = plist.get<double>(axis+" Width");
    checkAxis = true;
    dir = "Both";  // both axis directions by default
    // check if X/Y/Z Direction is specified (Both, Positive, or Negative)
    if (plist.isParameter(axis+" Direction"))
      dir = plist.get<string>(axis+" Direction");
    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  } else {
    loc = 0.0;
    width = 0.0;
    dir = "Both";
    min = -1e100;
    max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }
}



void erfcMoleFracParams::parseErfc(const Teuchos::ParameterList& plist, 
				   const int num_dim, const std::string& matArity) {
  using std::string;
  using Teuchos::ParameterList;

  mat_arity = matArity;

  if (plist.isParameter("Mole Max Value"))
    maxVal = plist.get<double>("Mole Max Value");
  else if (plist.isParameter("xMoleFrac Max Value"))
    maxVal = plist.get<double>("xMoleFrac Max Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Max Value' or 'xMoleFrac Max Value' must be specified!");
  if (plist.isParameter("Mole Min Value"))
    minVal = plist.get<double>("Mole Min Value");
  else if (plist.isParameter("xMoleFrac Min Value"))
    minVal = plist.get<double>("xMoleFrac Min Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Min Value' or 'xMoleFrac Min Value' must be specified!");
  if ((minVal < 0.) || (maxVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Erfc xMoleFrac Values must be within [0 1] !");
  if (mat_arity == "Quaternary") {
    if (!plist.isParameter("yMoleFrac Max Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac Max Value' must be specified!");
    if (!plist.isParameter("yMoleFrac Min Value"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	 << "Error ! 'yMoleFrac Min Value' must be specified!");
    maxVal1 = plist.get<double>("yMoleFrac Max Value");
    minVal1 = plist.get<double>("yMoleFrac Min Value");
    if ((minVal1 < 0.) || (maxVal1 > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Erfc yMoleFrac Values must be within [0 1] !");
  }

  testcoord("X", plist, x_dir, x_min, x_max, x_loc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_loc, z_width, z_checkAxis);
  }

}


void erfcMoleFracParams::testcoord(const std::string axis, 
				   const Teuchos::ParameterList& plist,
				   std::string& dir, double& min,
				   double& max, double& loc, double& width,
				   bool& checkAxis) {
  using std::string;
  using Teuchos::ParameterList;

  // axis+" Bend Location" and axis+" Width" must be specified together
  if (plist.isParameter(axis+" Bend Location") && !plist.isParameter(axis+" Width"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! " << axis << " Bend Location must be specified together with "
      << axis << " Width !" << std::endl);

  if (!plist.isParameter(axis+" Bend Location") && plist.isParameter(axis+" Width"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! " << axis << " Width must be specified together with "
      << axis << " Bend Location !" << std::endl);

  checkAxis = false;
  // Erfc is specified along a given axis
  if (plist.isParameter(axis+" Bend Location") && plist.isParameter(axis+" Width"))
  {
    loc = plist.get<double>(axis+" Bend Location");
    width = plist.get<double>(axis+" Width");

    checkAxis = true;
    dir = "Positive";  // positive axis directions by default

    // check if X/Y/Z Direction is specified (Positive or Negative)
    if (plist.isParameter(axis+" Direction"))
      dir = plist.get<string>(axis+" Direction");

    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }
}


void mgaussMoleFracParams::parseMGauss(const Teuchos::ParameterList& plist, 
				       const int num_dim, const std::string& matArity) {
  using std::string;
  using Teuchos::ParameterList;

  mat_arity = matArity;

  if (plist.isParameter("Mole Max Value"))
    maxVal = plist.get<double>("Mole Max Value");
  else if (plist.isParameter("xMoleFrac Max Value"))
    maxVal = plist.get<double>("xMoleFrac Max Value");
  else 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! 'Mole Max Value' must be specified!");
  minVal = 0.0;
  if (plist.isParameter("Mole Min Value"))
    minVal = plist.get<double>("Mole Min Value");
  else if (plist.isParameter("xMoleFrac Min Value"))
    minVal = plist.get<double>("xMoleFrac Min Value");
  if ((minVal < 0.) || (maxVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! MGauss xMoleFrac Values must be within [0 1] !");
  if (mat_arity == "Quaternary") {
    if (plist.isParameter("yMoleFrac Max Value"))
      maxVal1 = plist.get<double>("yMoleFrac Max Value");
    else 
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error ! 'yMoleFrac Max Value' must be specified!");
    minVal1 = 0.0;
    if (plist.isParameter("yMoleFrac Min Value"))
      minVal1 = plist.get<double>("Mole Min Value");
    if ((minVal1 < 0.) || (maxVal1 > 1.))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! MGauss yMoleFrac Values must be within [0 1] !");
  }

  testcoord("X", plist, x_dir, x_min, x_max, x_checkErfc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_checkErfc, z_width, z_checkAxis);
  }
}


void mgaussMoleFracParams::testcoord(
			       const std::string axis,
			       const Teuchos::ParameterList& plist,
			       std::string& dir, double& min,
			       double& max, bool& checkErfc, double& width,
			       bool& checkAxis) {
  using std::string;
  using Teuchos::ParameterList;

  // MGauss is specified along a given axis
  checkAxis = false;
  if (plist.isParameter(axis+" Width"))
  {
    width = plist.get<double>(axis+" Width");

    checkAxis = true;
    dir = "Both";  // both axis directions by default

    // check if X/Y/Z Direction is specified (Both, Positive, or Negative)
    if (plist.isParameter(axis+" Direction"))
      dir = plist.get<string>(axis+" Direction");

    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");

    checkErfc = plist.isParameter(axis+" ERFC_ON") ? plist.get<bool>(axis+" ERFC_ON") : false;
  }
}


void haloMoleFracParams::parseHalo(const Teuchos::ParameterList& plist, 
				   const int num_dim, const std::string& matArity) {
  using std::string;
  using Teuchos::ParameterList;

  mat_arity = matArity;

  val = 0.0;
  minVal = 0.0;
  if (plist.isParameter("Mole Value"))
    val = plist.get<double>("Mole Value");
  else if (plist.isParameter("xMoleFrac Value"))
    val = plist.get<double>("xMoleFrac Value");
  if (plist.isParameter("Mole Max Value"))
    val = plist.get<double>("Mole Max Value");
  else if (plist.isParameter("xMoleFrac Max Value"))
    val = plist.get<double>("xMoleFrac Max Value");
  if (plist.isParameter("Mole Min Value"))
    minVal = plist.get<double>("Mole Min Value");
  else if (plist.isParameter("xMoleFrac Min Value"))
    minVal = plist.get<double>("xMoleFrac Min Value");
  if ((val < 0.) || (val > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Halo xMoleFrac Values must be within [0 1] !");
  if ((minVal < 0.) || (minVal > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Halo xMoleFrac Values must be within [0 1] !");
  if (!plist.isParameter("Mole Distribution")) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	<< "Error! You must specify a distribution in the halo: Uniform or Gaussian.");
  } else {
    distributionType = plist.get<string>("Mole Distribution");
  }
  if(mat_arity == "Quaternary") {
    val1 = 0.0;
    minVal1 = 0.0;
    if (plist.isParameter("yMoleFrac Value"))
      val1 = plist.get<double>("yMoleFrac Value");
    if (plist.isParameter("yMoleFrac Max Value"))
      val1 = plist.get<double>("yMoleFrac Max Value");
    if (plist.isParameter("yMoleFrac Min Value"))
    minVal1 = plist.get<double>("yMoleFrac Min Value");
    if ((val1 < 0.) || (val1 > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Halo yMoleFrac Values must be within [0 1] !");
  if ((minVal1 < 0.) || (minVal1 > 1.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Halo yMoleFrac Values must be within [0 1] !");
  }

  if (!plist.isParameter("R1"))
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	      << "Error! The primary radius R1 must be specified for halo mole fraction.");
  if (!plist.isParameter("R2"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
	      << "Error! The primary radius R2 must be specified for halo mole fraction.");
  r1 = plist.get<double>("R1");
  r2 = plist.get<double>("R2");
  if (!plist.isParameter("X Center"))
    x_center = 0.0;
  else
    x_center = plist.get<double>("X Center");
  if (!plist.isParameter("Y Center"))
    y_center = 0.0;
  else
    y_center = plist.get<double>("Y Center");
  if (!plist.isParameter("Z Center"))
    z_center = 0.0;
  else
    z_center = plist.get<double>("Z Center");
  if (!plist.isParameter("Rotation"))
    rotation = 0.0;
  else
    rotation = plist.get<double>("Rotation");
  if (!plist.isParameter("Width"))
    width = 0.0;
  else
    width = plist.get<double>("Width");
}


void haloMoleFracParams::testcoord(
			  const std::string axis,
			  const Teuchos::ParameterList& plist) {
  using std::string;
  using Teuchos::ParameterList;
 
  //Blank for now.  May come back and figure this out later.  :LCM
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

  // mole fraction related
  moleFracParamList = p.sublist("Mole Fraction ParameterList");
  string matName = p.get<string>("Material Name");
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();
  mat_arity = matProperty.getArityType(matName);
  // flag for composition-aware material
  withMoleFrac = matProperty.hasMoleFracDependence(matName);
  
  // evaluated fields
  molefrac = MDField<ScalarT,Cell,IP>(n.field.mole_frac,scalar);
  molefrac_basis = MDField<ScalarT,Cell,BASIS>(n.field.mole_frac,basis_scalar);
  this->addEvaluatedField(molefrac);
  this->addEvaluatedField(molefrac_basis);

  if (withMoleFrac) {
    xMoleFrac = MDField<ScalarT,Cell,IP>(n.field.xMoleFrac,scalar);
    xMoleFrac_basis = MDField<ScalarT,Cell,BASIS>(n.field.xMoleFrac,basis_scalar);
    this->addEvaluatedField(xMoleFrac);
    this->addEvaluatedField(xMoleFrac_basis);
    if (mat_arity == "Quaternary") {
      yMoleFrac = MDField<ScalarT,Cell,IP>(n.field.yMoleFrac,scalar);
      yMoleFrac_basis = MDField<ScalarT,Cell,BASIS>(n.field.yMoleFrac,basis_scalar);
      this->addEvaluatedField(yMoleFrac);
      this->addEvaluatedField(yMoleFrac_basis);
    }
  }

  std::string name = "MoleFraction_Function";
  this->setName(name);

  for (ParameterList::ConstIterator model_it = moleFracParamList.begin();
       model_it != moleFracParamList.end(); ++model_it) {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0) {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      if (funcType == "Uniform")
      {
        uniformMoleFracParams umfp_;
        umfp_.parseUniform(funcParamList,mat_arity);
        umfp_vec.push_back(umfp_);
      }
      else if (funcType == "Linear")
      {
        linearMoleFracParams lmfp_;
        lmfp_.parseLinear(funcParamList,num_dim,mat_arity);
        lmfp_vec.push_back(lmfp_);
      }
      else if (funcType == "Gaussian" or funcType == "Gauss")
      {
	gaussMoleFracParams gmfp_;
        gmfp_.parseGaussian(funcParamList,num_dim,mat_arity);
	// create an equivalent doping params to 
	// be used with the profile calculator
	gaussianDopingParams gdp;
	gdp.dopType = "Acceptor";
	gdp.maxVal = gmfp_.maxVal;
	gdp.minVal = gmfp_.minVal;
	// x direction
	gdp.x_dir = gmfp_.x_dir;
	gdp.x_loc = gmfp_.x_loc;
	gdp.x_width = gmfp_.x_width;
	gdp.x_min = gmfp_.x_min;
	gdp.x_max = gmfp_.x_max;
	gdp.x_checkAxis = gmfp_.x_checkAxis;
	// y direction
	gdp.y_dir = gmfp_.y_dir;
	gdp.y_loc = gmfp_.y_loc;
	gdp.y_width = gmfp_.y_width;
	gdp.y_min = gmfp_.y_min;
	gdp.y_max = gmfp_.y_max;
	gdp.y_checkAxis = gmfp_.y_checkAxis;
	// z direction
	gdp.z_dir = gmfp_.z_dir;
	gdp.z_loc = gmfp_.z_loc;
	gdp.z_width = gmfp_.z_width;
	gdp.z_min = gmfp_.z_min;
	gdp.z_max = gmfp_.z_max;
	gdp.z_checkAxis = gmfp_.z_checkAxis;
	gmfp_equiv_vec.push_back(gdp);
	if (mat_arity == "Quaternary") {
	  gaussianDopingParams gdp1;
	  gdp1.dopType = "Acceptor";
	  gdp1.maxVal = gmfp_.maxVal1;
	  gdp1.minVal = gmfp_.minVal1;
	  // x direction
	  gdp1.x_dir = gmfp_.x_dir;
	  gdp1.x_loc = gmfp_.x_loc;
	  gdp1.x_width = gmfp_.x_width;
	  gdp1.x_min = gmfp_.x_min;
	  gdp1.x_max = gmfp_.x_max;
	  gdp1.x_checkAxis = gmfp_.x_checkAxis;
	  // y direction
	  gdp1.y_dir = gmfp_.y_dir;
	  gdp1.y_loc = gmfp_.y_loc;
	  gdp1.y_width = gmfp_.y_width;
	  gdp1.y_min = gmfp_.y_min;
	  gdp1.y_max = gmfp_.y_max;
	  gdp1.y_checkAxis = gmfp_.y_checkAxis;
	  // z direction
	  gdp1.z_dir = gmfp_.z_dir;
	  gdp1.z_loc = gmfp_.z_loc;
	  gdp1.z_width = gmfp_.z_width;
	  gdp1.z_min = gmfp_.z_min;
	  gdp1.z_max = gmfp_.z_max;
	  gdp1.z_checkAxis = gmfp_.z_checkAxis;
	  gmfp_equiv_vec1.push_back(gdp1);
	}
      }
      else if (funcType == "Erfc")
      {
	erfcMoleFracParams emfp_;
        emfp_.parseErfc(funcParamList,num_dim,mat_arity);
	// create an equivalent doping params to 
	// be used with the profile calculator
        erfcDopingParams edp;
	edp.dopType = "Acceptor";
	edp.maxVal = emfp_.maxVal;
	edp.minVal = emfp_.minVal;
	// x direction
	edp.x_dir = emfp_.x_dir;
	edp.x_loc = emfp_.x_loc;
	edp.x_width = emfp_.x_width;
	edp.x_min = emfp_.x_min;
	edp.x_max = emfp_.x_max;
	edp.x_checkAxis = emfp_.x_checkAxis;
	// y direction
	edp.y_dir = emfp_.y_dir;
	edp.y_loc = emfp_.y_loc;
	edp.y_width = emfp_.y_width;
	edp.y_min = emfp_.y_min;
	edp.y_max = emfp_.y_max;
	edp.y_checkAxis = emfp_.y_checkAxis;
	// z direction
	edp.z_dir = emfp_.z_dir;
	edp.z_loc = emfp_.z_loc;
	edp.z_width = emfp_.z_width;
	edp.z_min = emfp_.z_min;
	edp.z_max = emfp_.z_max;
	edp.z_checkAxis = emfp_.z_checkAxis;
	emfp_equiv_vec.push_back(edp);
	if (mat_arity == "Quaternary") {
	  erfcDopingParams edp1;
	  edp1.dopType = "Acceptor";
	  edp1.maxVal = emfp_.maxVal1;
	  edp1.minVal = emfp_.minVal1;
	  // x direction
	  edp1.x_dir = emfp_.x_dir;
	  edp1.x_loc = emfp_.x_loc;
	  edp1.x_width = emfp_.x_width;
	  edp1.x_min = emfp_.x_min;
	  edp1.x_max = emfp_.x_max;
	  edp1.x_checkAxis = emfp_.x_checkAxis;
	  // y direction
	  edp1.y_dir = emfp_.y_dir;
	  edp1.y_loc = emfp_.y_loc;
	  edp1.y_width = emfp_.y_width;
	  edp1.y_min = emfp_.y_min;
	  edp1.y_max = emfp_.y_max;
	  edp1.y_checkAxis = emfp_.y_checkAxis;
	  // z direction
	  edp1.z_dir = emfp_.z_dir;
	  edp1.z_loc = emfp_.z_loc;
	  edp1.z_width = emfp_.z_width;
	  edp1.z_min = emfp_.z_min;
	  edp1.z_max = emfp_.z_max;
	  edp1.z_checkAxis = emfp_.z_checkAxis;
	  emfp_equiv_vec1.push_back(edp1);
	}
      }
      else if (funcType == "MGauss")
      {
	mgaussMoleFracParams mmfp_;
        mmfp_.parseMGauss(funcParamList,num_dim,mat_arity);
	// create an equivalent doping params to 
	// be used with the profile calculator
        mgaussDopingParams mdp;
	mdp.dopType = "Acceptor";
	mdp.maxVal = mmfp_.maxVal;
	mdp.minVal = mmfp_.minVal;
	// x direction
	mdp.x_dir = mmfp_.x_dir;
	mdp.x_width = mmfp_.x_width;
	mdp.x_min = mmfp_.x_min;
	mdp.x_max = mmfp_.x_max;
	mdp.x_checkErfc = mmfp_.x_checkErfc;
	mdp.x_checkAxis = mmfp_.x_checkAxis;
	// y direction
	mdp.y_dir = mmfp_.y_dir;
	mdp.y_width = mmfp_.y_width;
	mdp.y_min = mmfp_.y_min;
	mdp.y_max = mmfp_.y_max;
	mdp.y_checkErfc = mmfp_.y_checkErfc;
	mdp.y_checkAxis = mmfp_.y_checkAxis;
	// z direction
	mdp.z_dir = mmfp_.z_dir;
	mdp.z_width = mmfp_.z_width;
	mdp.z_min = mmfp_.z_min;
	mdp.z_max = mmfp_.z_max;
	mdp.z_checkErfc = mmfp_.z_checkErfc;
	mdp.z_checkAxis = mmfp_.z_checkAxis;
	mmfp_equiv_vec.push_back(mdp);
	if (mat_arity == "Quaternary") {
	  mgaussDopingParams mdp1;
	  mdp1.dopType = "Acceptor";
	  mdp1.maxVal = mmfp_.maxVal1;
	  mdp1.minVal = mmfp_.minVal1;
	  // x direction
	  mdp1.x_dir = mmfp_.x_dir;
	  mdp1.x_width = mmfp_.x_width;
	  mdp1.x_min = mmfp_.x_min;
	  mdp1.x_max = mmfp_.x_max;
	  mdp1.x_checkErfc = mmfp_.x_checkErfc;
	  mdp1.x_checkAxis = mmfp_.x_checkAxis;
	  // y direction
	  mdp1.y_dir = mmfp_.y_dir;
	  mdp1.y_width = mmfp_.y_width;
	  mdp1.y_min = mmfp_.y_min;
	  mdp1.y_max = mmfp_.y_max;
	  mdp1.y_checkErfc = mmfp_.y_checkErfc;
	  mdp1.y_checkAxis = mmfp_.y_checkAxis;
	  // z direction
	  mdp1.z_dir = mmfp_.z_dir;
	  mdp1.z_width = mmfp_.z_width;
	  mdp1.z_min = mmfp_.z_min;
	  mdp1.z_max = mmfp_.z_max;
	  mdp1.z_checkErfc = mmfp_.z_checkErfc;
	  mdp1.z_checkAxis = mmfp_.z_checkAxis;
	  mmfp_equiv_vec1.push_back(mdp1);
	}
      }
      else if (funcType == "Halo")
      {
	haloMoleFracParams hmfp_;
        hmfp_.parseHalo(funcParamList,num_dim,mat_arity);
	// create an equivalent doping params to 
	// be used with the profile calculator
        haloDopingParams hdp;
	hdp.dopType = "Acceptor";
	hdp.distributionType = hmfp_.distributionType;
	hdp.dopingVal = hmfp_.val;
	hdp.minDopingVal = hmfp_.minVal;
	hdp.width = hmfp_.width;
	hdp.x_center = hmfp_.x_center;
	hdp.x_checkAxis = hmfp_.x_checkAxis;
	hdp.y_center = hmfp_.y_center;
	hdp.y_checkAxis = hmfp_.y_checkAxis;
	hdp.z_center = hmfp_.z_center;
	hdp.z_checkAxis = hmfp_.z_checkAxis;
	hdp.r1 = hmfp_.r1;
	hdp.r2 = hmfp_.r2;
	hdp.rotation = hmfp_.rotation;
	hmfp_equiv_vec.push_back(hdp);
	if (mat_arity == "Quaternary") {
	  haloDopingParams hdp1;
	  hdp1.dopType = "Acceptor";
	  hdp1.distributionType = hmfp_.distributionType;
	  hdp1.dopingVal = hmfp_.val1;
	  hdp1.minDopingVal = hmfp_.minVal1;
	  hdp1.width = hmfp_.width;
	  hdp1.x_center = hmfp_.x_center;
	  hdp1.x_checkAxis = hmfp_.x_checkAxis;
	  hdp1.y_center = hmfp_.y_center;
	  hdp1.y_checkAxis = hmfp_.y_checkAxis;
	  hdp1.z_center = hmfp_.z_center;
	  hdp1.z_checkAxis = hmfp_.z_checkAxis;
	  hdp1.r1 = hmfp_.r1;
	  hdp1.r2 = hmfp_.r2;
	  hdp1.rotation = hmfp_.rotation;
	  hmfp_equiv_vec1.push_back(hdp1);
	}
      }
    } // end of if (key.compare(0, 8, "Function") == 0)
  }
  
  // enable profile calculator if needed
  if (gmfp_equiv_vec.size() > 0 || emfp_equiv_vec.size() > 0 ||
      mmfp_equiv_vec.size() > 0 || hmfp_equiv_vec.size() > 0)
    prof_eval = Teuchos::rcp(new ProfileEvals(num_dim));
  
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
    // x,y mole fractions at IPs
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

      // evaluate the x,y mole fraction
      double value = evaluate_xMoleFraction(x,y,z);
      molefrac(cell,ip) = value;
      if (withMoleFrac) {
	xMoleFrac(cell,ip) = value;
	if (mat_arity == "Quaternary") 
	  yMoleFrac(cell,ip) = evaluate_yMoleFraction(x,y,z);
      }
    }

    // x,y mole fractions at basis points
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

      double value = evaluate_xMoleFraction(x,y,z);
      molefrac_basis(cell,basis) = value;
      if (withMoleFrac) {
	xMoleFrac_basis(cell,basis) = value;
	if (mat_arity == "Quaternary") 
	  yMoleFrac_basis(cell,basis) = evaluate_yMoleFraction(x,y,z);
      }
    }

  } // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluate_xMoleFraction()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evaluate_xMoleFraction
  (const double& x, const double& y, const double& z)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  double tempVal = 0.0;

  for (std::size_t i = 0; i < umfp_vec.size(); ++i)
  {
    tempVal = evalUniform_xMoleFrac(x,y,z,umfp_vec[i]);
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < lmfp_vec.size(); ++i)
  {
    tempVal = evalLinear_xMoleFrac(x,y,z,lmfp_vec[i]);
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < gmfp_equiv_vec.size(); ++i)
  {
    tempVal = prof_eval->evalGaussianProfile(x,y,z,gmfp_equiv_vec[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < emfp_equiv_vec.size(); ++i)
  {
    tempVal = prof_eval->evalErfcProfile(x,y,z,emfp_equiv_vec[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < mmfp_equiv_vec.size(); ++i)
  {
    tempVal = prof_eval->evalMGaussProfile(x,y,z,mmfp_equiv_vec[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < hmfp_equiv_vec.size(); ++i)
  {
    tempVal = prof_eval->evalHaloProfile(x,y,z,hmfp_equiv_vec[i])[0];
    mfValue += tempVal;
  }

  return mfValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluate_yMoleFraction()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evaluate_yMoleFraction
  (const double& x, const double& y, const double& z)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  double tempVal = 0.0;

  for (std::size_t i = 0; i < umfp_vec.size(); ++i)
  {
    tempVal = evalUniform_yMoleFrac(x,y,z,umfp_vec[i]);
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < lmfp_vec.size(); ++i)
  {
    tempVal = evalLinear_yMoleFrac(x,y,z,lmfp_vec[i]);
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < gmfp_equiv_vec1.size(); ++i)
  {
    tempVal = prof_eval->evalGaussianProfile(x,y,z,gmfp_equiv_vec1[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < emfp_equiv_vec1.size(); ++i)
  {
    tempVal = prof_eval->evalErfcProfile(x,y,z,emfp_equiv_vec1[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < mmfp_equiv_vec1.size(); ++i)
  {
    tempVal = prof_eval->evalMGaussProfile(x,y,z,mmfp_equiv_vec1[i])[0];
    mfValue += tempVal;
  }
  for (std::size_t i = 0; i < hmfp_equiv_vec1.size(); ++i)
  {
    tempVal = prof_eval->evalHaloProfile(x,y,z,hmfp_equiv_vec1[i])[0];
    mfValue += tempVal;
  }

  return mfValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalUniform_xMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalUniform_xMoleFrac
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
//  evalUniform_yMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalUniform_yMoleFrac
  (const double& x, const double& y, const double& z, const uniformMoleFracParams& umfp)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  const double val = umfp.value1;
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
//  evalLinear_xMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalLinear_xMoleFrac
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
//  evalLinear_yMoleFrac()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double MoleFraction_Function<EvalT, Traits>::evalLinear_yMoleFrac
  (const double& x, const double& y, const double& z, const linearMoleFracParams& lmfp)
{
  using std::string;
  using Teuchos::ParameterList;

  double mfValue = 0.0;
  const double startVal = lmfp.startVal1;
  const double endVal = lmfp.endVal1;

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
