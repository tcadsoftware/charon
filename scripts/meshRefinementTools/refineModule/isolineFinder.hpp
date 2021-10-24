
#ifndef ISOLINEFINDER_H
#define ISOLINEFINDER_H

#include <vector>
#include "cgeometries.hpp"

//This class takes a four-noded quad cell and computes isolines of a given value across it.


class isolineFinder
{

public:

  isolineFinder(){}

  void setCell(std::vector<double> x_, std::vector<double> y_, std::vector<double> fValues_)
  {
    x = x_;
    y = y_;
    f = fValues_;
  }

  std::vector<surfaceInfo> createSurfaces(double isovalue);

private:

  std::vector<double> x,y,f;
  double isoValue;

  std::vector<surfaceInfo> twoNodeSolution();
  std::vector<surfaceInfo> oneNodeSolution(int sign);

};


#endif

