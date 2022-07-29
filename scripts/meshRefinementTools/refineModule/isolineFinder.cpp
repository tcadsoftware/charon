
#include "isolineFinder.hpp"
#include <vector>
#include <iostream>
#include <cmath>

std::vector<surfaceInfo> isolineFinder::createSurfaces(double levelValue)
{

  std::vector<surfaceInfo> localSurfs;


  //First, shift the function values to zero by levelValue

  for(std::size_t fv=0 ; fv<f.size() ; ++fv)
    f[fv] -= levelValue;

  int posCount = 0;
  int negCount = 0;
  int zeroCount = 0;

  for(std::size_t fv=0 ; fv<f.size() ; fv++)
    if( f[fv] <= 0.0)
      negCount++;

  for(std::size_t fv=0 ; fv<f.size() ; fv++)
    if( f[fv] >= 0.0)
      posCount++;

  for(std::size_t fv=0 ; fv<f.size() ; fv++)
    if( f[fv] == 0.0)
      zeroCount++;

  if(negCount==1)
    localSurfs = oneNodeSolution(-1);

  if(posCount==1)
    localSurfs = oneNodeSolution(1);

  if(negCount == 2)  //same solution if posCount == 2
    localSurfs = twoNodeSolution();




  return localSurfs;

}



std::vector<surfaceInfo> isolineFinder::twoNodeSolution()
{

  std::vector<surfaceInfo> localSurfs;
  surfaceInfo oneSurf;

  //First need to figure out if the two nodes of same sign are successive--This will most often be the case

  int nodes[] = {0,1,2,3,0};
  bool consecutive=true;
  int node1;  //This is the start of the first consecutive sequence
  int node2;  //This is the start of the second consecutive sequence
  double fVal=0.0;

  int sign=1;
  if(f[0] < 0.0)
    sign = -1;

  if(sign > 0 && f[2] > 0)
    consecutive = false;

  if(consecutive)  // This cell will have only one surface element
    {
      if(sign == -1)
	{
	  if(f[1] > 0)
	    {
	      node1 = 0;
	      node2 = 2;
	    }
	  else
	    {
	      node1 = 1;
	      node2 = 3;
	    }
	}
      else
	{
	  if(f[1] < 0)
	    {
	      node1 = 0;
	      node2 = 2;
	    }
	  else
	    {
	      node1 = 1;
	      node2 = 3;
	    }
	}
      
      //Find the x,y locations of te zero values
      double xloc,yloc;
      bool printStuff = false;

      if(sign == -1 && node1 == 1 && printStuff)
	{
	  std::cout<<"Finding intersections "<<node1<<"    "<<node2<<"    "<<nodes[node1]<<"    "<<nodes[node2]<<std::endl;
	  std::cout<<" of Values "<<f[0]<<"   "<<f[1]<<"    "<<f[2]<<"    "<<f[3]<<std::endl;
	  std::cout<<" of Values "<<f[nodes[node1+1]]<<"   "<<f[nodes[node1]]<<"    "<<f[nodes[node2+1]]<<"    "<<f[nodes[node2]]<<std::endl;
	}
      double fraction = (fVal - f[nodes[node1]])/(f[nodes[node1+1]] - f[nodes[node1]]);
      xloc = x[node1] + fraction*(x[nodes[node1+1]] - x[nodes[node1]]);
      yloc = y[node1] + fraction*(y[nodes[node1+1]] - y[nodes[node1]]);
      oneSurf.X.push_back(xloc);
      oneSurf.Y.push_back(yloc);
      oneSurf.Z.push_back(0);

      if(sign == -1 && node1 == 1 && printStuff)
	{
	  std::cout<<" Fraction 1  "<<fraction<<std::endl;
	}

      fraction = (fVal - f[nodes[node2]])/(f[nodes[node2+1]] - f[nodes[node2]]);
      xloc = x[node2] + fraction*(x[nodes[node2+1]] - x[nodes[node2]]);
      yloc = y[node2] + fraction*(y[nodes[node2+1]] - y[nodes[node2]]);
      oneSurf.X.push_back(xloc);
      oneSurf.Y.push_back(yloc);
      oneSurf.Z.push_back(0);

      if(sign == -1 && node1 == 1 && printStuff)
	{
	  std::cout<<" Fraction 2  "<<fraction<<std::endl;
	}

      oneSurf.xc = 0.5*(oneSurf.X[0] + oneSurf.X[1]);
      oneSurf.yc = 0.5*(oneSurf.Y[0] + oneSurf.Y[1]);
      oneSurf.zc = 0;
      
      oneSurf.normx =   oneSurf.Y[1] - oneSurf.Y[0];
      oneSurf.normy = -(oneSurf.X[1] - oneSurf.X[0]);
      oneSurf.normz = 0.0;

      //Make normal unit
      double normSize = sqrt(oneSurf.normx*oneSurf.normx + oneSurf.normy*oneSurf.normy + oneSurf.normz*oneSurf.normz);
      oneSurf.normx = oneSurf.normx/normSize;
      oneSurf.normy = oneSurf.normy/normSize;-(oneSurf.X[1] - oneSurf.X[0]);
      oneSurf.normz = oneSurf.normz/normSize;
      
      localSurfs.push_back(oneSurf);

    }    

  return localSurfs;

}

std::vector<surfaceInfo> isolineFinder::oneNodeSolution(int sign)
{

  std::vector<surfaceInfo> localSurfs;
  surfaceInfo oneSurf;

  //Find the odd node first--if sign=-1, odd node is negative, otherwise positive
  int oddNode;
  for(int in=0 ; in<4 ; ++in)
    {
      if(sign*f[in] > 0)
	oddNode = in;
    }

  int nodes[] = {3,0,1,2,3,0};
  double fVal=0.0;

  int node1 = nodes[oddNode+1];
  int node2 = nodes[oddNode+2];

  bool debug = false;

  if(debug)
    {
      std::cout<<oddNode<<"    "<<node1<<"    "<<node2<<std::endl;
      std::cout<<nodes[node1-1]<<"    "<<node1<<"    "<<nodes[node1+1]<<std::endl;
      std::cout<<f[0]<<"    "<<f[1]<<"    "<<f[2]<<"    "<<f[3]<<std::endl;
      std::cout<<x[0]<<"    "<<x[1]<<"    "<<x[2]<<"    "<<x[3]<<std::endl;
      std::cout<<y[0]<<"    "<<y[1]<<"    "<<y[2]<<"    "<<y[3]<<std::endl;
    }

  double fraction = (fVal - f[node1])/(f[node2] - f[node1]);
  double xloc = x[node1] + fraction*(x[node2] - x[node1]);
  double yloc = y[node1] + fraction*(y[node2] - y[node1]);
  oneSurf.X.push_back(xloc);
  oneSurf.Y.push_back(yloc);
  oneSurf.Z.push_back(0);


  node2 = nodes[oddNode];

  if(debug)
    {
      std::cout<<" Node1 "<<xloc<<"    "<<yloc<<std::endl;
      std::cout<<node1<<"    "<<node2<<std::endl;
    }

  fraction = (fVal - f[node1])/(f[node2] - f[node1]);
  xloc = x[node1] + fraction*(x[node2] - x[node1]);
  yloc = y[node1] + fraction*(y[node2] - y[node1]);
  oneSurf.X.push_back(xloc);
  oneSurf.Y.push_back(yloc);
  oneSurf.Z.push_back(0);
  
  if(debug)
    {
      std::cout<<" Node2 "<<xloc<<"    "<<yloc<<std::endl;
      std::cout<<std::endl;
    }

  oneSurf.xc = 0.5*(oneSurf.X[0] + oneSurf.X[1]);
  oneSurf.yc = 0.5*(oneSurf.Y[0] + oneSurf.Y[1]);
  oneSurf.zc = 0;
  
  oneSurf.normx =   oneSurf.Y[1] - oneSurf.Y[0];
  oneSurf.normy = -(oneSurf.X[1] - oneSurf.X[0]);
  oneSurf.normz = 0.0;

  //Make normal unit
  double normSize = sqrt(oneSurf.normx*oneSurf.normx + oneSurf.normy*oneSurf.normy + oneSurf.normz*oneSurf.normz);
  oneSurf.normx = oneSurf.normx/normSize;
  oneSurf.normy = oneSurf.normy/normSize;-(oneSurf.X[1] - oneSurf.X[0]);
  oneSurf.normz = oneSurf.normz/normSize;

  localSurfs.push_back(oneSurf);

  return localSurfs;

}
