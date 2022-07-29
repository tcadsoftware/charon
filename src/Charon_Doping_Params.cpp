
#include "Charon_Doping_Params.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"

#include <iostream>

void charon::uniformDopingParams::parseUniform (const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;
  if (plist.isParameter("SweepMe"))
    if (plist.get<bool>("SweepMe"))
      {
	sweepMe = plist.get<bool>("SweepMe");
      }


  std::string funcType = plist.get<string>("Function Type");
  dopType = plist.get<string>("Doping Type");
  if (not sweepMe)
    dopVal = plist.get<double>("Doping Value");
    
  xmin = -1e100, ymin = -1e100, zmin = -1e100;
  xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if (plist.isParameter("Xmin"))  xmin = plist.get<double>("Xmin");
  if (plist.isParameter("Xmax"))  xmax = plist.get<double>("Xmax");
  if (plist.isParameter("Ymin"))  ymin = plist.get<double>("Ymin");
  if (plist.isParameter("Ymax"))  ymax = plist.get<double>("Ymax");
  if (plist.isParameter("Zmin"))  zmin = plist.get<double>("Zmin");
  if (plist.isParameter("Zmax"))  zmax = plist.get<double>("Zmax");
  if (plist.isParameter("Initial Doping Value"))  initialDopVal = plist.get<double>("Initial Doping Value");
  if (plist.isParameter("Final Doping Value"))    finalDopVal   = plist.get<double>("Final Doping Value");
}

void charon::gaussianDopingParams::parseGaussian (const Teuchos::ParameterList& plist,
                                                  const int num_dim)
{
  using std::string;
  using Teuchos::ParameterList;

  dopType = plist.get<string>("Doping Type");
  maxVal = plist.get<double>("Doping Max Value");
  minVal = plist.get<double>("Doping Min Value");

  if ((maxVal < 0.) || (minVal < 0.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gaussian doping Max and Min Values must be greater than 0.");


  testcoord("X", plist, x_dir, x_min, x_max, x_loc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_loc, z_width, z_checkAxis);
  }
}

void charon::gaussianDopingParams::testcoord (const std::string axis,
                                              const Teuchos::ParameterList& plist,
                                              std::string& dir, double& min,
                                              double& max, double& loc, double& width,
                                              bool& checkAxis)
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
    if (plist.isParameter(axis+" Direction"))
      dir = plist.get<string>(axis+" Direction");

    min = -1e100, max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }

  else
  {
    loc = 0.0;
    width = 0.0;
    dir = "Both";
    min = -1e100;
    max = 1e100;
    if (plist.isParameter(axis+"min")) min = plist.get<double>(axis+"min");
    if (plist.isParameter(axis+"max")) max = plist.get<double>(axis+"max");
  }

}

void charon::linearDopingParams::parseLinear (const Teuchos::ParameterList& plist, const int num_dim)
{

  using std::string;
  using Teuchos::ParameterList;

  dopType = plist.get<string>("Doping Type");
  startVal = plist.get<double>("Doping Start Value");
  endVal = plist.get<double>("Doping End Value");

  if ((startVal < 0.) || (endVal < 0.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Linear doping Start and End Values must be greater than 0.");

  testcoord("X", plist, x_min, x_max, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_min, y_max, y_checkAxis);
    testcoord("Z", plist, z_min, z_max, z_checkAxis);
  }

}

void charon::linearDopingParams::testcoord (const std::string axis, const Teuchos::ParameterList& plist,
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

void charon::erfcDopingParams::parseErfc (const Teuchos::ParameterList& plist,
                                          const int num_dim)
{
  using std::string;
  using Teuchos::ParameterList;


  dopType = plist.get<string>("Doping Type");
  maxVal = plist.get<double>("Doping Max Value");
  minVal = plist.get<double>("Doping Min Value");
  if ((maxVal < 0.) || (minVal < 0.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gauassian doping Max and Min Values must be greater than 0.");

  testcoord("X", plist, x_dir, x_min, x_max, x_loc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_loc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_loc, z_width, z_checkAxis);
  }
}

void charon::erfcDopingParams::testcoord (const std::string axis,
                                          const Teuchos::ParameterList& plist,
                                          std::string& dir, double& min,
                                          double& max, double& loc, double& width,
                                          bool& checkAxis)
{
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

void charon::mgaussDopingParams::parseMGauss (const Teuchos::ParameterList& plist,
                                                    const int num_dim)
{
  using std::string;
  using Teuchos::ParameterList;

  dopType = plist.get<string>("Doping Type");
  maxVal = plist.get<double>("Doping Max Value");

  minVal = 0.0;
  if (plist.isParameter("Doping Min Value"))
    minVal = plist.get<double>("Doping Min Value");

  if ((maxVal < 0.) )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Gaussian doping Max value must be greater than 0.");

  testcoord("X", plist, x_dir, x_min, x_max, x_checkErfc, x_width, x_checkAxis);
  if (num_dim == 2)
    testcoord("Y", plist, y_dir, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    testcoord("Y", plist, y_dir, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
    testcoord("Z", plist, z_dir, z_min, z_max, z_checkErfc, z_width, z_checkAxis);
  }
}

void charon::mgaussDopingParams::testcoord (const std::string axis,
                                            const Teuchos::ParameterList& plist,
                                            std::string& dir, double& min,
                                            double& max, bool& checkErfc, double& width,
                                            bool& checkAxis)
{
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


void charon::haloDopingParams::parseHalo (const Teuchos::ParameterList& plist,
					  const int num_dim)
{
  using std::string;
  using Teuchos::ParameterList;

  dopType = plist.get<string>("Doping Type");
  dopingVal = 0.0;
  minDopingVal = 0.0;
  if (plist.isParameter("Doping Value"))
    dopingVal = plist.get<double>("Doping Value");

  if (plist.isParameter("Doping Max Value"))
    dopingVal = plist.get<double>("Doping Max Value");

 if (plist.isParameter("Doping Min Value"))
    minDopingVal = plist.get<double>("Doping Min Value");

  if (!plist.isParameter("R1"))
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
					 << "Error!  The primary radius R1 must be specified for halo doping.");

  if (!plist.isParameter("R2"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
			       << "Error!  The primary radius R2 must be specified for halo doping.");

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

  if (!plist.isParameter("Doping Distribution"))
    {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
			       << "Error!  You must specify a distribution in the halo.  Uniform or Gaussian.");
    }
  else
    distributionType = plist.get<string>("Doping Distribution");

  if (!plist.isParameter("Width"))
    width = 0.0;
  else
    width = plist.get<double>("Width");


}

void charon::haloDopingParams::testcoord (const std::string axis,
					  const Teuchos::ParameterList& plist)
{
  using std::string;
  using Teuchos::ParameterList;
 
  //Blank for now.  May come back and figure this out later.  :LCM


}

