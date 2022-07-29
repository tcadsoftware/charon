
#ifndef DISTANCEFUNCTIONS_HPP
#define DISTANCEFUNCTIONS_HPP

#include <cmath>
#include <vector>
#include "Charon_SurfaceGeometry.hpp"


class distanceFunctions
{

public:

  double twoPointDistance(double x1,double y1,double x2,double y2)
  {
    double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    return distance;
  }

  double twoPointDistance(double x1,double y1,double z1,double x2,double y2,double z2)
  {
    double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    return distance;
  }

  double normalDistanceToLine(double xpoint, double ypoint, double& xnew, double& ynew,
			      std::vector<double> x_loc, std::vector<double> y_loc);

  double distanceToLine(double xpoint, double ypoint, double& xnew, double& ynew,
			std::vector<double> x_loc, std::vector<double> y_loc);

  double normalDistanceToSurface(double x, double y, double z, 
				 double& xnew, double& ynew, double& znew,
				 surfaceInfo surf);

  double distanceToSurface(double x, double y, double z, 
			   double& xnew, double& ynew, double& znew,
			   surfaceInfo surf);



};

#endif
