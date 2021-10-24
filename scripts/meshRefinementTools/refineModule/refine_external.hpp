
#include "cgeometries.hpp"
#include <cmath>
#include <string>
#include <algorithm>
#include "scalarFunction.hpp"
#include "isolineFinder.hpp"

class meshRefine
{

public:

  int meshDimension;
  int sizeMeasure;  //0-average, 1-minimum, 2-maximum, 3-average of max 2 sides
  int numSurfs;
  std::vector<surfaceInfo> surfs;
  //double **cellNodeCoords;
  std::vector< std::vector<double> > cellNodeCoords;
  int numCellNodes;
  int nodeCounter;
  double refinementDistance;
  double refinementDistanceN;
  double refinementDistanceP;
  double refinementDistanceIL;
  double refinementFactor;

  double BBxmin,BBxmax;
  double BBymin,BBymax;
  double BBzmin,BBzmax;

  
  double cellMinSize;
  double cellMinSizeN;
  double cellMinSizeP;
  double cellMinSizeIL;
  bool cellMinNSet;
  bool cellMinPSet;
  bool cellMinILSet;
  bool is3D;
  bool surfacesFileRead;
  bool writeJunctions;

  int maxSurfaceRecursions,guaranteedSurfaceRecursions;
  int tetNum;

  scalarFunction function;

  //Temporary
  double cxc,cyc,czc;
  bool flagTet;

  //Ctor
  meshRefine();

  void helloWorld();

  //void init(surfaceInfo *surfs);
  void init(int flag3D);

  void printSurfs();

  double * createDouble(int size);
  void sizeCellNodes(int numCellNodes_);
  void fillCoordinates(std::vector<double> x_);
  void resetNodeCounter();
  void printCoords();
  void freeCoords();
  void freeSurfs();
  void setRefinementDistance(double dist_);
  double getRefinementDistance();
  bool doIRefine();
  bool doIrefineXplane(double x);
  bool doIrefineCentroid();
  bool doIrefineCentroidSided();
  bool doIrefine3D();
  bool doIrefine2D();
  bool doIrefine2DDepricated();
  double doIrefineCentroidSidedFinishMetric();
  int getNumSurfs();
  double getMaxSide();
  double getMinSide();
  double getAveSide(bool);
  void setCellMinimum(double min);
  void setCellMinimumN(double min);
  void setCellMinimumP(double min);
  void setCellMinimumIL(double min);
  void setRefinementDistanceN(double dist_);
  void setRefinementDistanceP(double dist_);
  void setRefinementDistanceIL(double dist_);
  void setRefinementFactor(double rF);
  void findClosest(int numSubs,std::vector<int> &subSurfList);
  double getVectorRadAngle(double x1, double y1, double z1, 
			   double x2, double y2, double z2,
			   double normx, double normy, double normz);
  
  double distance(double x1,double y1,double z1,double x2,double y2,double z2);
  
  void testX(std::vector<double> x);
  void printTestX(std::vector<double> x);
  bool signedDistanceSameSideCheck(std::vector<double> pointm1, std::vector<double> point0, std::vector<double> point1, 
				   std::vector<double> point2, std::vector<double> point3);
  int checkSide(std::vector<double> norm, std::vector<double> point, std::vector<double> coords);
  int getSurfType(double x, double y, double z, surfaceInfo surf);
  
  void readSurfaces(std::string);
  void read2D(std::string filename);
  void read3D(std::string filename);

  void setDimension(int dim);
  int getDimension();
  void setSizeMeasure(std::string sM);
  double getLengthMetric();

  void createSurfaces();
  void createSurfaces2D();
  void createSurfaces3D();
  void journalSurfaces();

  void setDopingFunction(std::string name, std::string functionType, std::string dopingType);
  void setXBounds(std::string name, std::string axis, double min, double max);
  void setDopingBounds(std::string name, double min, double max);
  void setDirection(std::string name, std::string direction);
  void setDopingWidth(std::string name,  double xWidth, double yWidth, double zWidth); 
  void setDopingWidth(std::string name,  double xWidth, double yWidth); 
  void setDopingLocation(std::string name, double x, double y, double z);
  void setDopingLocation(std::string name, double x, double y);

  void addRefineToLine(double xmin, double ymin, double xmax, double ymax, double ilThick);
  void addRefineToSurface(double x0, double y0, double z0, double x1, double y1, double z1, 
			  double x2, double y2, double z2, double x3, double y3, double z3,
			  double ilThick);
  void setRefinementBoundingBox(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
  {
    BBxmin = xmin;
    BBymin = ymin;
    BBzmin = zmin;
    BBxmax = xmax;
    BBymax = ymax;
    BBzmax = zmax;
  }
  void setTetNum(int tn)
  {
    tetNum = tn;
  }
  int getTetNum(){return tetNum;}


  void listFunctions();

  void setMaxSurfaceRecursions(int level){maxSurfaceRecursions=level;}
  void setGuaranteedSurfaceRecursions(int level){guaranteedSurfaceRecursions=level;}
  void setWriteJunctions(bool setJunctions);
  void printXMLDopingFunctions();

};

