
#include <string>
#include <iostream>
#include <fstream>

class functionDefinitions
{

public:

  functionDefinitions(std::string name, std::string fType, std::string dType) :
    functionName(name),
    functionType(fType),
    dopingType(dType),
    functionDirection("positive"),
    calcX(false),
    calcY(false),
    calcZ(false),
    xmin(-1e10),
    ymin(-1e10),
    zmin(-1e10),
    xmax(1e10),
    ymax(1e10),
    zmax(1e10),
    dopingMin(0),
    dopingMax(0),
    xFunctionWidth(0),
    yFunctionWidth(0),
    zFunctionWidth(0),
    xFunctionMiddleLocation(0),
    yFunctionMiddleLocation(0),
    zFunctionMiddleLocation(0),
    xminDefined(false),
    yminDefined(false),
    zminDefined(false),
    xmaxDefined(false),
    ymaxDefined(false),
    zmaxDefined(false),
    dopingMinDefined(false),
    dopingMaxDefined(false),
    xFunctionWidthDefined(false),
    yFunctionWidthDefined(false),
    zFunctionWidthDefined(false),
    xFunctionMiddleLocationDefined(false),
    yFunctionMiddleLocationDefined(false),
    zFunctionMiddleLocationDefined(false)
  {  }

  void setXBounds(std::string Axis, double min, double max);
  void setDopingBounds(double min, double max);
  void setDirection(std::string direction){functionDirection = direction;}

  std::string getName(){return functionName;}


  virtual double evaluateFunction(double x, double y) = 0;
  virtual double evaluateFunction(double x, double y, double z) = 0;
  virtual void printXML(std::ofstream& ofile, int functionIndex) = 0;

  void setWidth(double xWidth, double yWidth, double zWidth)
  {
    xFunctionWidth = xWidth;
    yFunctionWidth = yWidth;
    zFunctionWidth = zWidth;
    xFunctionWidthDefined = true;
    yFunctionWidthDefined = true;
    zFunctionWidthDefined = true;
  }

  void setWidth(double xWidth, double yWidth)
  {
    xFunctionWidth = xWidth;
    yFunctionWidth = yWidth;
    xFunctionWidthDefined = true;
    yFunctionWidthDefined = true;
    zFunctionWidthDefined = false;
  }

  void setLocation(double xLocation, double yLocation, double zLocation)
  {
    xFunctionMiddleLocation = xLocation;
    yFunctionMiddleLocation = yLocation;
    zFunctionMiddleLocation = zLocation;
    xFunctionMiddleLocationDefined = true;
    yFunctionMiddleLocationDefined = true;
    zFunctionMiddleLocationDefined = true;
  }

  void setLocation(double xLocation, double yLocation)
  {
    xFunctionMiddleLocation = xLocation;
    yFunctionMiddleLocation = yLocation;
    xFunctionMiddleLocationDefined = true;
    yFunctionMiddleLocationDefined = true;
    zFunctionMiddleLocationDefined = false;
  }

  friend std::ostream& operator<<(
				  std::ostream& os,
				  const functionDefinitions& obj)
  {
    os<<obj.functionName<<std::endl;
    os<<obj.functionType<<std::endl;
    os<<obj.dopingType<<std::endl;
    os<<"X Bounds "<<obj.xmin<<"  to  "<<obj.xmax<<std::endl;
    os<<"Y Bounds "<<obj.ymin<<"  to  "<<obj.ymax<<std::endl;
    os<<"Z Bounds "<<obj.zmin<<"  to  "<<obj.zmax<<std::endl;
    os<<"Doping Bounds "<<obj.dopingMin<<"    "<<obj.dopingMax<<std::endl;
    os<<"X Y Z Doping Widths "<<obj.xFunctionWidth<<"  "<<obj.yFunctionWidth<<"    "<<obj.zFunctionWidth<<std::endl;
    os<<"X Y Z Bend Locations "<<obj.xFunctionMiddleLocation<<"  "<<obj.yFunctionMiddleLocation<<"    "<<obj.zFunctionMiddleLocation<<std::endl;
    os<<"End of function "<<obj.functionName<<std::endl;
    os<<std::endl;
    return os;
  }


protected:

  std::string functionName,functionType,functionAxis,functionDirection,dopingType;
  double dopingMin,dopingMax;
  double xmin,ymin,zmin,xmax,ymax,zmax;
  bool calcX,calcY,calcZ;

  double xFunctionWidth,yFunctionWidth;
  double zFunctionWidth,xFunctionMiddleLocation;
  double yFunctionMiddleLocation,zFunctionMiddleLocation;

  bool xminDefined,yminDefined,zminDefined;
  bool xmaxDefined,ymaxDefined,zmaxDefined;
  bool dopingMinDefined,dopingMaxDefined;
  bool xFunctionWidthDefined,yFunctionWidthDefined,zFunctionWidthDefined;
  bool xFunctionMiddleLocationDefined,yFunctionMiddleLocationDefined,zFunctionMiddleLocationDefined;

} ;

//----------------------------------------------------------------------------------
// erfc function
//----------------------------------------------------------------------------------

class erfcFunctionDefinition : public functionDefinitions
{

public:

  erfcFunctionDefinition(std::string name, std::string fType, std::string dType) : 
    functionDefinitions(name, fType, dType)
  {}


  double evaluateFunction(double x, double y);
  double evaluateFunction(double x, double y, double z);
  void printXML(std::ofstream& ofile, int functionIndex);

private:

} ;

//----------------------------------------------------------------------------------
// MGauss function
//----------------------------------------------------------------------------------

class mGaussFunctionDefinition : public functionDefinitions
{

public:

  mGaussFunctionDefinition(std::string name, std::string fType, std::string dType) : 
    functionDefinitions(name, fType, dType)
  {}


  double evaluateFunction(double x, double y);
  double evaluateFunction(double x, double y, double z);
  void printXML(std::ofstream& ofile, int functionIndex);

private:

} ;


//----------------------------------------------------------------------------------
// Uniform function
//----------------------------------------------------------------------------------

class uniformFunctionDefinition : public functionDefinitions
{

public:

  uniformFunctionDefinition(std::string name, std::string fType, std::string dType) : 
    functionDefinitions(name, fType, dType)
  {}
  

  double evaluateFunction(double x, double y);
  double evaluateFunction(double x, double y, double z);
  void printXML(std::ofstream& ofile, int functionIndex);

private:

} ;

