

#ifndef SCALARFUNCTION_H
#define SCALARFUNCTION_H

#include <vector>
#include <string>
#include "functionDefinitions.hpp"
#include <iostream>
#include <fstream>

class scalarFunction
{

  double epsilon;

public:

  scalarFunction() : epsilon(1e-6)
  {}

  double evaluateFunction(std::vector<double> x);

  void setDopingFunction(std::string name, std::string functionType, std::string dopingType);
  void setXBounds(std::string name, std::string axis, double min, double max);
  void setDopingBounds(std::string name, double min, double max);
  void setDirection(std::string name, std::string direction);
  void setLocation(std::string name, double x, double y, double z);
  void setLocation(std::string name, double x, double y);
  void setWidth(std::string name, double x, double y, double z);
  void setWidth(std::string name, double x, double y);

  void listFunctions()
  {
    std::cout<<" Listing the details of "<<dopingFunctions.size()<<" Function(s)"<<std::endl;
    for(size_t func=0 ; func<dopingFunctions.size() ; ++func)
      std::cout<<*dopingFunctions[func]<<std::endl;
    std::cout<<" Functions listed "<<std::endl;
  }


  void printXMLFunctions();

private:

  std::vector<functionDefinitions*> dopingFunctions;

};


#endif
