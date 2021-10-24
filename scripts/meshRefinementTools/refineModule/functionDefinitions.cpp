

#include <cmath>
#include "functionDefinitions.hpp"
#include <iostream>

//-----------------------------------------------------------------------
// set bounds of function not sure this will be necessary
//-----------------------------------------------------------------------
void functionDefinitions::setXBounds(std::string Axis, double min, double max)
{
  
  if(Axis == "x")
    {
      xmin = min;
      xmax = max;
      if(xmin != xmax)
	calcX = true;
      xminDefined = true;
      xmaxDefined = true;
    }
  else if(Axis == "y")
    {
      ymin = min;
      ymax = max;
      if(ymin != ymax)
	calcY = true;
      yminDefined = true;
      ymaxDefined = true;
    }
  else if(Axis == "z")
    {
      zmin = min;
      zmax = max;
      if(zmin != zmax)
	calcZ = true;
      zminDefined = true;
      zmaxDefined = true;
    }
  else
    {
      //error out
      std::cout<<" Illegal coordinate axis set in the doping function, "<<functionName<<std::endl;
    }
}

//-----------------------------------------------------------------------
// set the minimum and maximum doping concentrations for this function
//-----------------------------------------------------------------------

void functionDefinitions::setDopingBounds(double min, double max)
{
  dopingMin = min;
  dopingMax = max;
  dopingMinDefined = true;
  dopingMaxDefined = true;
}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//   ERFC METHODS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Spit out the xml equivalent syntax of the erfc doping function
//-----------------------------------------------------------------------

void erfcFunctionDefinition::printXML(std::ofstream& ofile, int functionIndex)
{

  ofile<<" <ParameterList name=\"Function"<<functionIndex<<"\">"<<std::endl;

  ofile<<"   <Parameter name=\"Function Type\" type=\"string\" value=\"Erfc\"/>"<<std::endl;

  ofile<<"   <Parameter name=\"Doping Type\" type=\"string\" value=\""<<dopingType<<"\"/>"<<std::endl;

  if(dopingMinDefined)
    ofile<<"   <Parameter name=\"Doping Min Value\" type=\"double\" value=\""<<dopingMin<<"\"/>"<<std::endl;

  if(dopingMaxDefined)
    ofile<<"   <Parameter name=\"Doping Max Value\" type=\"double\" value=\""<<dopingMax<<"\"/>"<<std::endl;

  //ofile<<"   <Parameter name=\"X Direction\" type=\"string\" value=\"Positive\"/>"<<std::endl;

  if(xFunctionWidthDefined)
    ofile<<"   <Parameter name=\"X Width\" type=\"double\" value=\""<<xFunctionWidth<<"\"/>"<<std::endl;

  if(yFunctionWidthDefined)
    ofile<<"   <Parameter name=\"Y Width\" type=\"double\" value=\""<<yFunctionWidth<<"\"/>"<<std::endl;

  if(zFunctionWidthDefined)
    ofile<<"   <Parameter name=\"Z Width\" type=\"double\" value=\""<<zFunctionWidth<<"\"/>"<<std::endl;


  if(xFunctionMiddleLocationDefined)
    ofile<<"   <Parameter name=\"X Bend Location\" type=\"double\" value=\""<<xFunctionMiddleLocation<<"\"/>"<<std::endl;

  if(yFunctionMiddleLocationDefined)
    ofile<<"   <Parameter name=\"Y Bend Location\" type=\"double\" value=\""<<yFunctionMiddleLocation<<"\"/>"<<std::endl;

  if(zFunctionMiddleLocationDefined)
    ofile<<"   <Parameter name=\"Z Bend Location\" type=\"double\" value=\""<<zFunctionMiddleLocation<<"\"/>"<<std::endl;

  ofile<<" </ParameterList>"<<std::endl;

}


//-----------------------------------------------------------------------
// evaluate the doping concentration 2D
//-----------------------------------------------------------------------

double erfcFunctionDefinition::evaluateFunction(double x, double y)
{

  double xerfc = 1.0;
  double yerfc = 1.0;

  double sign = 1.0;
  if(dopingType =="acceptor")
    sign = -1.0;

  double functionArgument = (x-xFunctionMiddleLocation)/xFunctionWidth;
  if (functionDirection == "negative")
    functionArgument *= -1.0;

  if(calcX)
    xerfc = 0.5*erfc(functionArgument);

  functionArgument = (y-yFunctionMiddleLocation)/yFunctionWidth;
  if (functionDirection == "negative")
    functionArgument *= -1.0;

  if(calcY)
    yerfc = 0.5*erfc(functionArgument);

  double argument1 = dopingMax/dopingMin;
  double argument2 = xerfc*yerfc;
  double doping = sign*dopingMin*pow(argument1,argument2);

  return doping;

}



//-----------------------------------------------------------------------
// evaluate the doping concentration 3D
//-----------------------------------------------------------------------

double erfcFunctionDefinition::evaluateFunction(double x, double y, double z)
{

  double xerfc = 1.0;
  double yerfc = 1.0;
  double zerfc = 1.0;

  double sign = 1.0;
  if(dopingType =="acceptor")
    sign = -1.0;

  double functionArgument = (x-xFunctionMiddleLocation)/xFunctionWidth;
  if (functionDirection == "negative")
    functionArgument *= -1.0;

  if(calcX)
    xerfc = 0.5*erfc(functionArgument);

  functionArgument = (y-yFunctionMiddleLocation)/yFunctionWidth;
  if (functionDirection == "negative")
    functionArgument *= -1.0;

  if(calcY)
    yerfc = 0.5*erfc(functionArgument);

  functionArgument = (z-zFunctionMiddleLocation)/zFunctionWidth;
  if (functionDirection == "negative")
    functionArgument *= -1.0;

  if(calcZ)
    zerfc = 0.5*erfc(functionArgument);

  double argument1 = dopingMax/dopingMin;
  double argument2 = xerfc*yerfc*zerfc;
  double doping = sign*dopingMin*pow(argument1,argument2);

  return doping;

}



//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//   MGAUSS METHODS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Spit out the xml equivalent syntax of the doping function
//-----------------------------------------------------------------------

void mGaussFunctionDefinition::printXML(std::ofstream& ofile, int functionIndex)
{

  ofile<<" <ParameterList name=\"Function"<<functionIndex<<"\">"<<std::endl;
  ofile<<"   <Parameter name=\"Function Type\" type=\"string\" value=\"MGauss\"/>"<<std::endl;

  ofile<<"   <Parameter name=\"Doping Type\" type=\"string\" value=\""<<dopingType<<"\"/>"<<std::endl;

  if(dopingMinDefined)
    ofile<<"   <Parameter name=\"Doping Min Value\" type=\"double\" value=\""<<dopingMin<<"\"/>"<<std::endl;

  if(dopingMaxDefined)
    ofile<<"   <Parameter name=\"Doping Max Value\" type=\"double\" value=\""<<dopingMax<<"\"/>"<<std::endl;

  if(xminDefined)
    ofile<<"   <Parameter name=\"Xmin\" type=\"double\" value=\""<<xmin<<"\"/>"<<std::endl;

  if(xminDefined)
    ofile<<"   <Parameter name=\"Xmax\" type=\"double\" value=\""<<xmax<<"\"/>"<<std::endl;

  if(yminDefined)
    ofile<<"   <Parameter name=\"Ymin\" type=\"double\" value=\""<<ymin<<"\"/>"<<std::endl;

  if(ymaxDefined)
    ofile<<"   <Parameter name=\"Ymax\" type=\"double\" value=\""<<ymax<<"\"/>"<<std::endl;

  if(zminDefined)
    ofile<<"   <Parameter name=\"Zmin\" type=\"double\" value=\""<<zmin<<"\"/>"<<std::endl;

  if(zmaxDefined)
    ofile<<"   <Parameter name=\"Zmax\" type=\"double\" value=\""<<zmax<<"\"/>"<<std::endl;

  //ofile<<"   <Parameter name=\"X Direction\" type=\"string\" value=\"Positive\"/>"<<std::endl;

  if(xFunctionWidthDefined)
    ofile<<"   <Parameter name=\"X Width\" type=\"double\" value=\""<<xFunctionWidth<<"\"/>"<<std::endl;

  if(yFunctionWidthDefined)
    ofile<<"   <Parameter name=\"Y Width\" type=\"double\" value=\""<<yFunctionWidth<<"\"/>"<<std::endl;

  if(zFunctionWidthDefined)
    ofile<<"   <Parameter name=\"Z Width\" type=\"double\" value=\""<<zFunctionWidth<<"\"/>"<<std::endl;

  ofile<<"   <Parameter name=\"X ERFC_ON\" type=\"bool\" value=\"false\"/>"<<std::endl;

  ofile<<" </ParameterList>"<<std::endl;

}


//-----------------------------------------------------------------------
// evaluate the doping concentration 2D
//-----------------------------------------------------------------------

double mGaussFunctionDefinition::evaluateFunction(double x, double y)
{

  double xmGauss = 1.0;
  double ymGauss = 1.0;

  double sign = 1.0;
  if(dopingType =="acceptor")
    sign = -1.0;

  if(x <= xmin)
    {
      if(dopingMin > 0.0)
	{
	  xmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((x-xmin)/xFunctionWidth, 2.0));
	}
      else
	{
	  xmGauss = std::exp(-(x-xmin)*(x-xmin)/xFunctionWidth/xFunctionWidth);
	}
    }
  else if (x >= xmax)
    {
      if(dopingMin > 0.0)
	{
	  xmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((x-xmax)/xFunctionWidth, 2.0));
	}
      else
	{
	  xmGauss = std::exp(-(x-xmax)*(x-xmax)/xFunctionWidth/xFunctionWidth);
	}
    }
  else
    {
      //xmGauss = 0.5*(erfc((x-xmax)/xFunctionWidth)-erfc((x-xmin)/xFunctionWidth));
    }


  if(y <= ymin)
    {
      if(dopingMin > 0.0)
	{
	  ymGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((y-ymin)/yFunctionWidth, 2.0));
	}
      else
	{
	  ymGauss = std::exp(-(y-ymin)*(y-ymin)/yFunctionWidth/yFunctionWidth);
	}
    }
  else if (y >= ymax)
    {
      if(dopingMin > 0.0)
	{
	  ymGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((y-ymax)/yFunctionWidth, 2.0));
	}
      else
	{
	  ymGauss = std::exp(-(y-ymax)*(y-ymax)/yFunctionWidth/yFunctionWidth);
	}
    }
  else
    {
      //ymGauss = 0.5*(erfc((y-ymax)/yFunctionWidth)-erfc((y-ymin)/yFunctionWidth));
    }

  double doping = sign*dopingMax*xmGauss*ymGauss;

  return doping;

}



//-----------------------------------------------------------------------
// evaluate the doping concentration 3D
//-----------------------------------------------------------------------

double mGaussFunctionDefinition::evaluateFunction(double x, double y, double z)
{

  double xmGauss = 1.0;
  double ymGauss = 1.0;
  double zmGauss = 1.0;

  double sign = 1.0;
  if(dopingType =="acceptor")
    sign = -1.0;

  if(x < xmin)
    {
      if(dopingMin > 0.0)
	{
	  xmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((x-xmin)/xFunctionWidth, 2.0));
	}
      else
	{
	  xmGauss = std::exp(-(x-xmin)*(x-xmin)/xFunctionWidth/xFunctionWidth);
	}
    }
  else if (x > xmax)
    {
      if(dopingMin > 0.0)
	{
	  xmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((x-xmax)/xFunctionWidth, 2.0));
	}
      else
	{
	  xmGauss = std::exp(-(x-xmax)*(x-xmax)/xFunctionWidth/xFunctionWidth);
	}
    }
  else
    {
      xmGauss = 0.5*(erfc((x-xmax)/xFunctionWidth)-erfc((x-xmin)/xFunctionWidth));
    }
  

  if(y < ymin)
    {
      if(dopingMin > 0.0)
	{
	  ymGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((y-ymin)/yFunctionWidth, 2.0));
	}
      else
	{
	  ymGauss = std::exp(-(y-ymin)*(y-ymin)/yFunctionWidth/yFunctionWidth);
	}
    }
  else if (y > ymax)
    {
      if(dopingMin > 0.0)
	{
	  ymGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((y-ymax)/yFunctionWidth, 2.0));
	}
      else
	{
	  ymGauss = std::exp(-(y-ymax)*(y-ymax)/yFunctionWidth/yFunctionWidth);
	}
    }
  else
    {
      ymGauss = 0.5*(erfc((y-ymax)/yFunctionWidth)-erfc((y-ymin)/yFunctionWidth));
    }


  if(z < zmin)
    {
      if(dopingMin > 0.0)
	{
	  zmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((z-zmin)/zFunctionWidth, 2.0));
	}
      else
	{
	  zmGauss = std::exp(-(z-zmin)*(z-zmin)/zFunctionWidth/zFunctionWidth);
	}
    }
  else if (z > zmax)
    {
      if(dopingMin > 0.0)
	{
	  zmGauss = std::exp(-std::log(dopingMax/dopingMin) * std::pow((z-zmax)/zFunctionWidth, 2.0));
	}
      else
	{
	  zmGauss = std::exp(-(z-zmax)*(z-zmax)/zFunctionWidth/zFunctionWidth);
	}
    }
  else
    {
      zmGauss = 0.5*(erfc((z-zmax)/zFunctionWidth)-erfc((z-zmin)/zFunctionWidth));
    }

  double doping = sign*dopingMax*xmGauss*ymGauss*zmGauss;

  //  std::cout<<" My doping = "<<doping<<std::endl;

  return doping;


}


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//   UNIFORM METHODS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Spit out the xml equivalent syntax of the erfc doping function
//-----------------------------------------------------------------------

void uniformFunctionDefinition::printXML(std::ofstream& ofile, int functionIndex)
{

  ofile<<" <ParameterList name=\"Function"<<functionIndex<<"\">"<<std::endl;

  ofile<<"   <Parameter name=\"Function Type\" type=\"string\" value=\"Uniform\"/>"<<std::endl;

  ofile<<"   <Parameter name=\"Doping Type\" type=\"string\" value=\""<<dopingType<<"\"/>"<<std::endl;

  if(dopingMaxDefined)
    ofile<<"   <Parameter name=\"Doping Value\" type=\"double\" value=\""<<dopingMax<<"\"/>"<<std::endl;

  if(xminDefined)
    ofile<<"   <Parameter name=\"Xmin\" type=\"double\" value=\""<<xmin<<"\"/>"<<std::endl;

  if(xminDefined)
    ofile<<"   <Parameter name=\"Xmax\" type=\"double\" value=\""<<xmax<<"\"/>"<<std::endl;

  if(yminDefined)
    ofile<<"   <Parameter name=\"Ymin\" type=\"double\" value=\""<<ymin<<"\"/>"<<std::endl;

  if(ymaxDefined)
    ofile<<"   <Parameter name=\"Ymax\" type=\"double\" value=\""<<ymax<<"\"/>"<<std::endl;

  if(zminDefined)
    ofile<<"   <Parameter name=\"Zmin\" type=\"double\" value=\""<<zmin<<"\"/>"<<std::endl;

  if(zmaxDefined)
    ofile<<"   <Parameter name=\"Zmax\" type=\"double\" value=\""<<zmax<<"\"/>"<<std::endl;

  ofile<<" </ParameterList>"<<std::endl;

}


//-----------------------------------------------------------------------
// evaluate the uniform doping concentration 2D
//-----------------------------------------------------------------------

double uniformFunctionDefinition::evaluateFunction(double x, double y)
{

  double sign = 1.0;
  double doping;

  if(dopingType =="acceptor")
    sign = -1.0;

  if (x <= xmin || x >= xmax ||
      y <= ymin || y >= ymax)
    doping = 0.0;
  else
    doping = sign*dopingMin;

  return doping;

}


//-----------------------------------------------------------------------
// evaluate the uniform doping concentration 3D
//-----------------------------------------------------------------------

double uniformFunctionDefinition::evaluateFunction(double x, double y, double z)
{

  double sign = 1.0;
  double doping;

  if(dopingType =="acceptor")
    sign = -1.0;

  if (x <= xmin || x >= xmax ||
      y <= ymin || y >= ymax ||
      z <= zmin || z >= zmax)
    doping = 0.0;
  else
    doping = sign*dopingMin;

  return doping;

}


