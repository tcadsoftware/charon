
#ifndef GEOMETRIES
#define GEOMETRIES

#include <vector>

class surfaceInfo
{

public:

  int numNodes;
  double xc,yc,zc;
  double normx,normy,normz;
  int dimension;
  std::vector<double> X,Y,Z;

  void setNode(double x, double y, double z)
  {
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);
  }

  void setNode(double x, double y)
  {
    X.push_back(x);
    Y.push_back(y);
  }

  void clearNode()
  {
    X.clear();
    Y.clear();
    Z.clear();
  }

  ~surfaceInfo()
  {
    X.clear();
    Y.clear();
    Z.clear();
  }

};

#endif
