
#include "Charon_TransformSurface.hpp"
#include <iostream>
#include <sstream>
#include <cmath>

lcm_lib::Transform::Transform()
{
  xform_check = 0;  // zero if direction cosines have not been computed
}

lcm_lib::Transform::~Transform()
{
}

void lcm_lib::Transform::transform(double x,double y,double z,double x1,
	       double y1,double z1,double x2,double y2,
	       double z2,double& xx0,double& yy0,double& xx1,
	       double& yy1,double& xx2,double& yy2)
{

  double r1[3],r2[3],r2p[3],r3[3];

  double point1[3],point2[3];

  //First step is to translate the origin of the coordinate 
  //system to point x1,y1,z1

  //Store the translation values for later

  translation[0] = x1;
  translation[1] = y1;
  translation[2] = z1;

  x -= x1;
  y -= y1;
  z -= z1;

  x2 -= x1;
  y2 -= y1;
  z2 -= z1;

  x1 = 0;
  y1 = 0;
  z1 = 0;

  //compute bases of new coordinate system 


  r1[0] = x2-x1;
  r1[1] = y2-y1;
  r1[2] = z2-z1;

  r2p[0] = x-x1;
  r2p[1] = y-y1;
  r2p[2] = z-z1;

  normalize(r1,3);
  normalize(r2p,3);

  cross_product(r1,r2p,r3,3);

  normalize(r3,3);

  cross_product(r3,r1,r2,3);

  normalize(r2,3);

  //set basis vectors of original system and prepare for transformation

  double i1[3],i2[3],i3[3];

  i1[0] = 1.0;
  i1[1] = 0;
  i1[2] = 0;

  i2[0] = 0;
  i2[1] = 1.0;
  i2[2] = 0;

  i3[0] = 0;
  i3[1] = 0;
  i3[2] = 1.0;

  //compute direction cosines for transformation

  get_dir_cos(r1,r2,r3,i1,i2,i3,3);

  //set point two equal to x y z and transform

  point2[0] = x;
  point2[1] = y;
  point2[2] = z;

  mat_vec(point1,point2,3,3);

  xx0 = point1[0];
  yy0 = point1[1];

  //set point two equal to x1 y1 z1 and transform

  point2[0] = x1;
  point2[1] = y1;
  point2[2] = z1;

  mat_vec(point1,point2,3,3);

  xx1 = point1[0];
  yy1 = point1[1];

  //set point two equal to x2 y2 z2 and transform

  point2[0] = x2;
  point2[1] = y2;
  point2[2] = z2;

  mat_vec(point1,point2,3,3);

  xx2 = point1[0];
  yy2 = point1[1];

  return;

}



void lcm_lib::Transform::normalize(double vec[], int size)
{
  double magnitude;

  magnitude = sqrt(dot_product(vec,vec,size));

  if(magnitude==0.0)return;

  for(int i=0 ; i<size ; ++i)
    vec[i] = vec[i]/magnitude;

  return;
}


double lcm_lib::Transform::dot_product(double vec1[], double vec2[], int size)
{

  double sum = 0.0;

  for(int i=0 ; i<size ; ++i)
    sum += vec1[i]*vec2[i];

  return sum;

}



void lcm_lib::Transform::cross_product(double vec1[], double vec2[],
				       double vec3[], int size)
{

  if(size != 3)return;

  vec3[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  vec3[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  vec3[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return;

}


void  lcm_lib::Transform::get_dir_cos( double r1[], double r2[], 
				       double r3[], double i1[], 
				       double i2[], double i3[],
				       int size)
{

  if(size != 3)return;


  dir_cos[0][0] = dot_product(r1,i1,size);
  dir_cos[1][0] = dot_product(r2,i1,size);
  dir_cos[2][0] = dot_product(r3,i1,size);
  dir_cos[0][1] = dot_product(r1,i2,size);
  dir_cos[1][1] = dot_product(r2,i2,size);
  dir_cos[2][1] = dot_product(r3,i2,size);
  dir_cos[0][2] = dot_product(r1,i3,size);
  dir_cos[1][2] = dot_product(r2,i3,size);
  dir_cos[2][2] = dot_product(r3,i3,size);

  xform_check = 1;

  return;

}

void lcm_lib::Transform::mat_vec(double vec1[], double vec2[], 
				 int size1, int size2)
{

  if(size2 != 3)return;

  if(xform_check == 0)
    {
      std::cout<<"ERROR: Transform::mat_vec: Cannot call this until directions cosines "<<
	"have been computed "<<std::endl;

      return;
    }

  for(int i=0 ; i<size1 ; ++i)
    {
      vec1[i] = 0.0;
      for(int j=0 ; j<size2 ; ++j)
	vec1[i] += dir_cos[i][j]*vec2[j];
    }

  return;

}

void lcm_lib::Transform::trans_mat_vec(double vec1[], double vec2[], 
				       int size1, int size2)
{

  if(size2 != 3)return;

  if(xform_check == 0)
    {
      std::cout<<"ERROR: Transform::trans_mat_vec: Cannot call this until directions cosines "<<
	"have been computed "<<std::endl;
      return;
    }

  for(int i=0 ; i<size1 ; ++i)
    {
      vec1[i] = 0.0;
      for(int j=0 ; j<size2 ; ++j)
	vec1[i] += dir_cos[j][i]*vec2[j];
    }

  return;

}

void lcm_lib::Transform::planar_coords(double& x, double& y, double x0, 
				       double y0, double z0)
{

  if(xform_check == 0)
    {
      std::cout<<"ERROR: Transform::planar_coords: Cannot call this until directions cosines "<<
	"have been computed "<<std::endl;
      return;
    }

  double point1[3],point2[3];

  point2[0] = x0;
  point2[1] = y0;
  point2[2] = z0;

  //Perform the tranlation to new origin

  for(int i=0 ; i<3 ; ++i)
    point2[i] -= translation[i];

  mat_vec(point1,point2,3,3);

  x = point1[0];
  y = point1[1];

  return;

}

double lcm_lib::tri_area(double x0,double y0,double x1,
			 double y1,double x2,double y2)
{

  double area=0.0;

  area = 0.5*(x0*y1 + y0*x2 + y2*x1
	      - y1*x2 - y0*x1 - x0*y2);

  return fabs(area);

}


bool lcm_lib::tri_check(double x, double y, double x0,double y0,
			double x1,double y1,double x2,double y2)
{

  //This function returns 0 if point lies inside triangle
  //else it returns 1


  double area,area1,area2,area3;
  double area_sum;

  area  = lcm_lib::tri_area(x0,y0,x1,y1,x2,y2);
  area1 = lcm_lib::tri_area(x,y,x1,y1,x2,y2);
  area2 = lcm_lib::tri_area(x0,y0,x,y,x2,y2);
  area3 = lcm_lib::tri_area(x0,y0,x1,y1,x,y);

  area_sum = area1+area2+area3;

  if(fabs(area_sum-area) <= 1.e-10)
    return true;
  else
    return false;

}


void lcm_lib::Transform::deplanar_coords(double x, double y, double& x0, 
					 double& y0, double& z0)
{

  if(xform_check == 0)
    {
      std::cout<<"ERROR: Transform::planar_coords: Cannot call this until directions cosines "<<
	"have been computed "<<std::endl;
      return;
    }

  double point1[3],point2[3];

  point2[0] = x;
  point2[1] = y;
  point2[2] = 0;

  trans_mat_vec(point1,point2,3,3);

  //Perform the tranlation to new origin

  for(int i=0 ; i<3 ; ++i)
    point1[i] += translation[i];

  x0 = point1[0];
  y0 = point1[1];
  z0 = point1[2];

  return;

}

