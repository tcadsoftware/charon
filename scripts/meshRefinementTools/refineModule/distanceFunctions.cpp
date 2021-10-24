

#include "distanceFunctions.hpp"
#include <cmath>
#include <algorithm>
#include "transform.hpp"
#include <iostream>

double distanceFunctions::normalDistanceToLine(double xpoint, double ypoint, double& x_line, double& y_line,
					       std::vector<double> x_loc, const std::vector<double> y_loc)
{

  //Want this to be strictly normal distance...anything off the line segment doesnt count.

  /** This method computes the normal distance from a point to a line segment.  It is strictly a normal
      distance.  If there can't be drawn a line from the point to the line segment that makes a normal angle
      to the line segment, this method returns a distance value of 1e10.

      input: 
            <double> xpoint,ypoint - The popint in question
	    std::vector<double> x_loc,y_loc - size 2 arrays that contain the end points of the line segment

      return: <double> Normal distance from point to line.

  **/


  double A,B,C,value;
  double buffer = 1.e-11;
  double x1,x2,y1,y2;
  double xmin,xmax,ymin,ymax;
  //double x_line,y_line;
  double lam;

  x1 = x_loc[0];
  x2 = x_loc[1];
  y1 = y_loc[0];
  y2 = y_loc[1];


  double middist,xmid,ymid;

  xmid = (x1+x2)/2;
  ymid = (y1+y2)/2;

  middist = sqrt((xpoint-xmid)*(xpoint-xmid) + 
		 (ypoint-ymid)*(ypoint-ymid));

  //buffer = boxavg;

  if(middist > buffer)buffer = 0;

  A = y2-y1;
  B = -(x2-x1);
  C = (x2-x1)*y1 - (y2-y1)*x1;
  lam = -(2*(A*xpoint + B*ypoint + C))/(A*A+B*B);

  x_line = 0.5*lam*A + xpoint;
  y_line = 0.5*lam*B + ypoint;


  xmin = std::min(x1,x2);
  xmax = std::max(x1,x2);
  //Is this a perfectly vertical line?
  if(xmax == xmin)
    {
      ymin = std::min(y1,y2);
      ymax = std::max(y1,y2);
      if (y_line < ymin - buffer)  //The following blocks set x,y line equal to endpoint
	{
	  y_line = ymin;
	  if(y_line == y1)
	    x_line = x1;
	  else
	    x_line = x2;
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);
	  value = 1e10;
	}
      else if (y_line > ymax + buffer)
	{
	  y_line = ymax;
	  if(y_line == y1)
	    x_line = x1;
	  else
	    x_line = x2;
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);
	  value = 1e10;
	}
      else
	{
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);  //This actually checks a true normal distance
	  return value;
	}
    }
  else if (x_line < xmin - buffer)
    {
      x_line = xmin;
      if(x_line == x1)
	y_line = y1;
      else 
	y_line = y2;
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);
      value = 1e10;
    }
  else if (x_line > xmax + buffer)
    {
      x_line = xmax;
      if(x_line == x1)
	y_line = y1;
      else 
	y_line = y2;
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);
      value = 1e10;
    }
  else
    {
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);  //This actually checks a true normal distance
    }
  
  return value;
  
}


double distanceFunctions::distanceToLine(double xpoint, double ypoint, double& x_line, double& y_line,
					 std::vector<double> x_loc, const std::vector<double> y_loc)
{

  //This is a looser version of normal distance...distance off the line segment counts.

  /** This method computes the normal distance from a point to a line segment.  If there can't be drawn 
      a line from the point to the line segment that makes a normal angle to the line segment, this method 
      returns a distance value equal to the distance from the point to the nearest end point of the line.

      input: 
            <double> xpoint,ypoint - The popint in question
	    std::vector<double> x_loc,y_loc - size 2 arrays that contain the end points of the line segment

      return: <double> Distance from point to line.

  **/


  double A,B,C,value;
  double buffer = 1.e-11;
  double x1,x2,y1,y2;
  double xmin,xmax,ymin,ymax;
  //double x_line,y_line;
  double lam;

  x1 = x_loc[0];
  x2 = x_loc[1];
  y1 = y_loc[0];
  y2 = y_loc[1];


  double middist,xmid,ymid;

  xmid = (x1+x2)/2;
  ymid = (y1+y2)/2;

  middist = sqrt((xpoint-xmid)*(xpoint-xmid) + 
		 (ypoint-ymid)*(ypoint-ymid));

  //buffer = boxavg;

  if(middist > buffer)buffer = 0;

  A = y2-y1;
  B = -(x2-x1);
  C = (x2-x1)*y1 - (y2-y1)*x1;
  lam = -(2*(A*xpoint + B*ypoint + C))/(A*A+B*B);

  x_line = 0.5*lam*A + xpoint;
  y_line = 0.5*lam*B + ypoint;


  xmin = std::min(x1,x2);
  xmax = std::max(x1,x2);
  //Is this a perfectly vertical line?
  if(xmax == xmin)
    {
      ymin = std::min(y1,y2);
      ymax = std::max(y1,y2);
      if (y_line < ymin - buffer)  //The following blocks set x,y line equal to endpoint
	{
	  y_line = ymin;
	  if(y_line == y1)
	    x_line = x1;
	  else
	    x_line = x2;
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);
	}
      else if (y_line > ymax + buffer)
	{
	  y_line = ymax;
	  if(y_line == y1)
	    x_line = x1;
	  else
	    x_line = x2;
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);
	}
      else
	{
	  value = twoPointDistance(xpoint,ypoint,x_line,y_line);  //This actually checks a true normal distance
	  return value;
	}
    }
  else if (x_line < xmin - buffer)
    {
      x_line = xmin;
      if(x_line == x1)
	y_line = y1;
      else 
	y_line = y2;
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);
    }
  else if (x_line > xmax + buffer)
    {
      x_line = xmax;
      if(x_line == x1)
	y_line = y1;
      else 
	y_line = y2;
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);
    }
  else
    {
      value = twoPointDistance(xpoint,ypoint,x_line,y_line);  //This actually checks a true normal distance
    }
  
  return value;
  
}

//---------------------------------------------------------------------------------
//3D version of the normal distance from a point to a triangular surface element
//---------------------------------------------------------------------------------



double distanceFunctions::normalDistanceToSurface(double x, double y, double z,
						  double& xnew, double& ynew, double& znew,
						  surfaceInfo surf)
{

  /** This method computes the normal distance from a point to a line segment.  It is strictly a normal
      distance.  If there can't be drawn a line from the point to the line segment that makes a normal angle
      to the line segment, this method returns a distance value of 1e10.

      input: 
            <double> xpoint,ypoint - The popint in question
	    std::vector<double> x_loc,y_loc - size 2 arrays that contain the end points of the line segment

      return: <double> Normal distance from point to line.

  **/


  double A,B,C,D,lambda;
  double x1,x2,x3;
  double y1,y2,y3;
  double z1,z2,z3;
  double numer,denom;

  double xlo=-1.0e30,xhi=1.0e30;
  double ylo=-1.0e30,yhi=1.0e30;
  double zlo=-1.0e30,zhi=1.0e30;

  //find minimum distance as constrained min. problem.

  x1 = surf.X[0];
  y1 = surf.Y[0];
  z1 = surf.Z[0];

  x2 = surf.X[1];
  y2 = surf.Y[1];
  z2 = surf.Z[1];

  x3 = surf.X[2];
  y3 = surf.Y[2];
  z3 = surf.Z[2];

  A = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
  B = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1);
  C = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
  D = -(A*x1 + B*y1 + C*z1);

  numer = -2.0*(A*x + B*y + C*z + D);
  denom = (A*A + B*B + C*C);

  if(denom == 0)
    return 1.e30;

  lambda = numer/denom;

  xnew = 0.5*lambda*A + x;
  ynew = 0.5*lambda*B + y;
  znew = 0.5*lambda*C + z;

  double norm_distance;

  norm_distance = twoPointDistance(x,y,z,xnew,ynew,znew);

  //must transform into planar coord. system and
  //check to see if new point is in triangle

  double xp,yp,xpnew,ypnew,xp0,yp0,xp1,yp1,xp2,yp2;

  lcm_lib::Transform xformer;


  //establish transformation to check if new point lies in surface


  xformer.transform(x1,y1,z1,x2,y2,z2,x3,y3,z3,
		    xp2,yp2,xp0,yp0,xp1,yp1);

  xformer.planar_coords(xpnew,ypnew,xnew,ynew,znew);

  //Calculate distance from point to midpoint of surface

  double middist,xmid,ymid,zmid;

  xmid = (x1+x2+x3)/3;
  ymid = (y1+y2+y3)/3;
  zmid = (z1+z2+z3)/3;

  middist = twoPointDistance(x,y,z,xmid,ymid,zmid);

  //Set up a buffer around triangle

  double buffer,xp0b,yp0b,xp1b,yp1b,xp2b,yp2b;

  //buffer = boxavg;

  if(middist > buffer)buffer = 0;

  //Enlarge triangle by a buffer -- note:  There shouldn't be any dilation
  //in the transformation, so buffer area should be consistent

  double vec1x,vec1y,vecmag;

  if(buffer != 0)
    {
      vec1x = xp0 - xp1 + xp0 - xp2;
      vec1y = yp0 - yp1 + yp0 - yp2;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp0b = xp0 + buffer*vec1x/vecmag;
      yp0b = yp0 + buffer*vec1y/vecmag;

      vec1x = xp1 - xp0 + xp1 - xp2;
      vec1y = yp1 - yp0 + yp1 - yp2;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp1b = xp1 + buffer*vec1x/vecmag;
      yp1b = yp1 + buffer*vec1y/vecmag;

      vec1x = xp2 - xp1 + xp2 - xp0;
      vec1y = yp2 - yp1 + yp2 - yp0;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp2b = xp2 + buffer*vec1x/vecmag;
      yp2b = yp2 + buffer*vec1y/vecmag;
    }
  else
    {
      xp0b = xp0;
      yp0b = yp0;
      xp1b = xp1;
      yp1b = yp1;
      xp2b = xp2;
      yp2b = yp2;
    }

  //Now check to see if point lies in buffered triangle

  bool check_tri;

  check_tri =  lcm_lib::tri_check(xpnew,ypnew,xp0b,yp0b,xp1b,yp1b,xp2b,yp2b);

  if(!check_tri)return 1e30;

  //Check to see if xnew lies outside domain
  //use that distance if point is approximately in cell
  //with surf

  double ptsurfdist,surfside;

  surfside = 0;
  surfside += twoPointDistance(x1,y1,z1,x2,y2,z2);
  surfside += twoPointDistance(x1,y1,z1,x3,y3,z3);
  surfside += twoPointDistance(x3,y3,z3,x2,y2,z2);

  surfside = surfside/3;

  ptsurfdist = twoPointDistance(x,y,z,x1,y1,z1);
  ptsurfdist = std::min(ptsurfdist,twoPointDistance(x,y,z,x2,y2,z2));
  ptsurfdist = std::min(ptsurfdist,twoPointDistance(x,y,z,x3,y3,z3));

  //if(ptsurfdist/surfside < 1 && !initialLSCalculation)
  if(ptsurfdist/surfside < 1)
    {
      if(xnew < xlo || xnew > xhi)check_tri=true;
      if(ynew < ylo || ynew > yhi)check_tri=true;
      if(znew < zlo || znew > zhi)check_tri=true;
    }

  if(check_tri)
    return norm_distance;


  //if no normal distance exists, must find distance to each
  //of the three edges and take minimum.

  double dist1,dist2,dist3;

  std::vector<double> xCoords,yCoords;
  xCoords.resize(2);
  yCoords.resize(2);

  double min_dist,xnt,ynt,znt;


  //establish transform to coord system with side 1

  xformer.transform(x,y,z,x1,y1,z1,x2,y2,z2,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist1 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnew,ynew,znew);

  min_dist = twoPointDistance(x,y,z,xnew,ynew,znew);

  //establish transform to coord system with side 2

  xformer.transform(x,y,z,x1,y1,z1,x3,y3,z3,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist2 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnt,ynt,znt);

  dist2 = twoPointDistance(x,y,z,xnt,ynt,znt);

  if(dist2 < min_dist)
    {
      min_dist = dist2;
      xnew = xnt;
      ynew = ynt;
      znew = znt;
    }


  //establish transform to coord system with side 3

  xformer.transform(x,y,z,x3,y3,z3,x2,y2,z2,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist3 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnt,ynt,znt);

  dist3 = twoPointDistance(x,y,z,xnt,ynt,znt);

  if(dist3 < min_dist)
    {
      min_dist = dist3;
      xnew = xnt;
      ynew = ynt;
      znew = znt;
    }

  //Return

  return min_dist;

}

//---------------------------------------------------------------------------------
//3D version of the distance from a point to a triangular surface element
//---------------------------------------------------------------------------------



double distanceFunctions::distanceToSurface(double x, double y, double z,
						  double& xnew, double& ynew, double& znew,
						  surfaceInfo surf)
{

  /** This method computes the normal distance from a point to a line segment.  It is strictly a normal
      distance.  If there can't be drawn a line from the point to the line segment that makes a normal angle
      to the line segment, this method returns a distance value of 1e10.

      input: 
            <double> xpoint,ypoint - The popint in question
	    std::vector<double> x_loc,y_loc - size 2 arrays that contain the end points of the line segment

      return: <double> Normal distance from point to line.

  **/


  double A,B,C,D,lambda;
  double x1,x2,x3;
  double y1,y2,y3;
  double z1,z2,z3;
  double numer,denom;

  double xlo=-1.0e30,xhi=1.0e30;
  double ylo=-1.0e30,yhi=1.0e30;
  double zlo=-1.0e30,zhi=1.0e30;

  //find minimum distance as constrained min. problem.

  x1 = surf.X[0];
  y1 = surf.Y[0];
  z1 = surf.Z[0];

  x2 = surf.X[1];
  y2 = surf.Y[1];
  z2 = surf.Z[1];

  x3 = surf.X[2];
  y3 = surf.Y[2];
  z3 = surf.Z[2];

  A = (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
  B = (z2-z1)*(x3-x1) - (z3-z1)*(x2-x1);
  C = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
  D = -(A*x1 + B*y1 + C*z1);

  numer = -2.0*(A*x + B*y + C*z + D);
  denom = (A*A + B*B + C*C);

  if(denom == 0)
    return 1.e30;

  lambda = numer/denom;

  xnew = 0.5*lambda*A + x;
  ynew = 0.5*lambda*B + y;
  znew = 0.5*lambda*C + z;

  double norm_distance;

  norm_distance = twoPointDistance(x,y,z,xnew,ynew,znew);

  //must transform into planar coord. system and
  //check to see if new point is in triangle

  double xp,yp,xpnew,ypnew,xp0,yp0,xp1,yp1,xp2,yp2;

  lcm_lib::Transform xformer;


  //establish transformation to check if new point lies in surface


  xformer.transform(x1,y1,z1,x2,y2,z2,x3,y3,z3,
		    xp2,yp2,xp0,yp0,xp1,yp1);

  xformer.planar_coords(xpnew,ypnew,xnew,ynew,znew);

  //Calculate distance from point to midpoint of surface

  double middist,xmid,ymid,zmid;

  xmid = (x1+x2+x3)/3;
  ymid = (y1+y2+y3)/3;
  zmid = (z1+z2+z3)/3;

  middist = twoPointDistance(x,y,z,xmid,ymid,zmid);

  //Set up a buffer around triangle

  double buffer,xp0b,yp0b,xp1b,yp1b,xp2b,yp2b;

  //buffer = boxavg;

  if(middist > buffer)buffer = 0;

  //Enlarge triangle by a buffer -- note:  There shouldn't be any dilation
  //in the transformation, so buffer area should be consistent

  double vec1x,vec1y,vecmag;

  if(buffer != 0)
    {
      vec1x = xp0 - xp1 + xp0 - xp2;
      vec1y = yp0 - yp1 + yp0 - yp2;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp0b = xp0 + buffer*vec1x/vecmag;
      yp0b = yp0 + buffer*vec1y/vecmag;

      vec1x = xp1 - xp0 + xp1 - xp2;
      vec1y = yp1 - yp0 + yp1 - yp2;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp1b = xp1 + buffer*vec1x/vecmag;
      yp1b = yp1 + buffer*vec1y/vecmag;

      vec1x = xp2 - xp1 + xp2 - xp0;
      vec1y = yp2 - yp1 + yp2 - yp0;
      vecmag = sqrt(vec1x*vec1x + vec1y*vec1y);
      xp2b = xp2 + buffer*vec1x/vecmag;
      yp2b = yp2 + buffer*vec1y/vecmag;
    }
  else
    {
      xp0b = xp0;
      yp0b = yp0;
      xp1b = xp1;
      yp1b = yp1;
      xp2b = xp2;
      yp2b = yp2;
    }

  //Now check to see if point lies in buffered triangle

  bool check_tri;

  check_tri =  lcm_lib::tri_check(xpnew,ypnew,xp0b,yp0b,xp1b,yp1b,xp2b,yp2b);

  //Check to see if xnew lies outside domain
  //use that distance if point is approximately in cell
  //with surf

  double ptsurfdist,surfside;

  surfside = 0;
  surfside += twoPointDistance(x1,y1,z1,x2,y2,z2);
  surfside += twoPointDistance(x1,y1,z1,x3,y3,z3);
  surfside += twoPointDistance(x3,y3,z3,x2,y2,z2);

  surfside = surfside/3;

  ptsurfdist = twoPointDistance(x,y,z,x1,y1,z1);
  ptsurfdist = std::min(ptsurfdist,twoPointDistance(x,y,z,x2,y2,z2));
  ptsurfdist = std::min(ptsurfdist,twoPointDistance(x,y,z,x3,y3,z3));

  //if(ptsurfdist/surfside < 1 && !initialLSCalculation)
  /*
  if(ptsurfdist/surfside < 1)
    {
      if(xnew < xlo || xnew > xhi)check_tri=true;
      if(ynew < ylo || ynew > yhi)check_tri=true;
      if(znew < zlo || znew > zhi)check_tri=true;
    }
  */

  if(check_tri)
    return norm_distance;

  //if no normal distance exists, must find distance to each
  //of the three edges and take minimum.

  double dist1,dist2,dist3;

  std::vector<double> xCoords,yCoords;
  xCoords.resize(2);
  yCoords.resize(2);

  double min_dist,xnt,ynt,znt;


  //establish transform to coord system with side 1

  xformer.transform(x,y,z,x1,y1,z1,x2,y2,z2,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist1 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnew,ynew,znew);

  min_dist = twoPointDistance(x,y,z,xnew,ynew,znew);

  //establish transform to coord system with side 2

  xformer.transform(x,y,z,x1,y1,z1,x3,y3,z3,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist2 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnt,ynt,znt);

  dist2 = twoPointDistance(x,y,z,xnt,ynt,znt);

  if(dist2 < min_dist)
    {
      min_dist = dist2;
      xnew = xnt;
      ynew = ynt;
      znew = znt;
    }


  //establish transform to coord system with side 3

  xformer.transform(x,y,z,x3,y3,z3,x2,y2,z2,
		    xp,yp,xp0,yp0,xp1,yp1);

  xCoords[0] = xp0;
  xCoords[1] = xp1;
  yCoords[0] = yp0;
  yCoords[1] = yp1;

  dist3 = normalDistanceToLine(xp,yp,xpnew,ypnew,xCoords,yCoords);

  xformer.deplanar_coords(xpnew,ypnew,xnt,ynt,znt);

  dist3 = twoPointDistance(x,y,z,xnt,ynt,znt);

  if(dist3 < min_dist)
    {
      min_dist = dist3;
      xnew = xnt;
      ynew = ynt;
      znew = znt;
    }

  //Return

  return min_dist;

}

