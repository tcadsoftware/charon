
#ifndef GEOMETRIES
#define GEOMETRIES

#include <vector>

class surfaceInfo
{

public:

  int numNodes;
  double xc,yc,zc;
  double normx,normy,normz;
  int surfType;
  std::vector<double> X,Y,Z;
  double size;
  double nThick,pThick;
  double refinementDistance;

};

#endif
