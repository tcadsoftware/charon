
#include "scalarFunction.hpp"
#include "lusolve.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

double scalarFunction::evaluateFunction(std::vector<double> x)
{

  double value=0.0;

  for(size_t func=0 ; func<dopingFunctions.size() ; ++func)
    {  //Ugh
      if(x.size() == 2)
	{
	  value += dopingFunctions[func]->evaluateFunction(x[0],x[1]);
	}
      else
	{
	  value += dopingFunctions[func]->evaluateFunction(x[0],x[1],x[2]);
	}
    }

  return value;

}
//-----------------------------------------------------------------------
// write the xml doping function
//-----------------------------------------------------------------------

void scalarFunction::printXMLFunctions()
{

  std::ofstream ofile("xmlDopingFunctions.xml");

  for(size_t func=0 ; func<dopingFunctions.size() ; ++func)
    dopingFunctions[func]->printXML(ofile,func);

  ofile.close();

}


//-----------------------------------------------------------------------
// create the doping function
//-----------------------------------------------------------------------

void scalarFunction::setDopingFunction(std::string name, std::string functionType, std::string dopingType)
{

  if(dopingType != "acceptor" && dopingType != "donor")
    {
      std::cout<<" ERROR!! Ilegal doping type-> "<<name<<"    "<<dopingType<<std::endl;
      return;
    }

  if(functionType == "erfc")
    {
      dopingFunctions.push_back(new erfcFunctionDefinition(name, functionType, dopingType));
    }
  else if(functionType == "mgauss")
    {
      dopingFunctions.push_back(new mGaussFunctionDefinition(name, functionType, dopingType));
    }
  else if(functionType == "uniform")
    {
      dopingFunctions.push_back(new uniformFunctionDefinition(name, functionType, dopingType));
    }
  else
    {
      std::cout<<" ERROR!! Ilegal doping function definition-> "<<name<<"    "<<functionType<<std::endl;
      return;
    }
}

//-----------------------------------------------------------------------
// set the minimum and maximum geometric bounds
//-----------------------------------------------------------------------
void scalarFunction::setXBounds(std::string name, std::string axis, double min, double max)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setXBounds(axis, min, max);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}

//-----------------------------------------------------------------------
// set the minimum and maximum doping concentrations for this function
//-----------------------------------------------------------------------
void scalarFunction::setDopingBounds(std::string name, double min, double max)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setDopingBounds( min,  max);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}



//-----------------------------------------------------------------------
// set the coordinate direction (+/-) of the function
//-----------------------------------------------------------------------

void scalarFunction::setDirection(std::string name, std::string direction)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setDirection(direction);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}



//-----------------------------------------------------------------------
// set the function location--a translation from the origin
//-----------------------------------------------------------------------
void scalarFunction::setLocation(std::string name, double x, double y, double z)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setLocation(x, y, z);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}



//-----------------------------------------------------------------------
// set the function location--a translation from the origin
//-----------------------------------------------------------------------
void scalarFunction::setLocation(std::string name, double x, double y)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setLocation(x, y);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}



//-----------------------------------------------------------------------
// set the function width
//-----------------------------------------------------------------------
void scalarFunction::setWidth(std::string name, double x, double y, double z)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setWidth(x, y, z);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}



//-----------------------------------------------------------------------
// set the function width
//-----------------------------------------------------------------------
void scalarFunction::setWidth(std::string name, double x, double y)
{

  for(size_t func=0 ; func < dopingFunctions.size() ; ++func)
    {
      if(dopingFunctions[func]->getName() == name)
	{
	  dopingFunctions[func]->setWidth(x, y);
	  return;
	}
    }

  std::cout<<" ERROR!!!  I could not find a function with the name of "<<name<<std::endl;

}


