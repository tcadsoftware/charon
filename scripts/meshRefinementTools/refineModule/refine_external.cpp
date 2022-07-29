

#include "refine_external.hpp"
#include "cgeometries.hpp"
#include "distanceFunctions.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

//#include "table.h"  
//#include "tet_table.h"  
//#include "vector_macros.h"
#include "cellCreator.hpp"  
#include "cells.hpp"
#include "lusolve.hpp"  
#include "surfaceFinder.hpp"
#include "depletionWidth.hpp"

//Constructor

meshRefine::meshRefine() :
  meshDimension(0),
  maxSurfaceRecursions(0),
  guaranteedSurfaceRecursions(0),
  BBxmin(0),
  BBymin(0),
  BBzmin(0),
  BBxmax(0),
  BBymax(0),
  BBzmax(0),
  numSurfs(0),
  refinementFactor(3),
  semiType(0),
  autoRefine(false)
{}


void meshRefine::helloWorld()
{

  std::cout<<" Hello World!"<<std::endl;

}

void meshRefine::setDimension(int dim)
{
  meshDimension = dim;
}

int meshRefine::getDimension()
{
  return meshDimension;
}

void meshRefine::setRefinementFactor(double rF)
{
  refinementFactor = rF;
}

void meshRefine::init(int mD)
{

  cellMinNSet = false;
  cellMinPSet = false;
  cellMinILSet = false;
  surfacesFileRead = false;
  writeJunctions = false;
  meshDimension = mD;

  if(mD == 3)
    {
      is3D = true;
    }
  else
    {
      is3D = false;
    }


}

void meshRefine::setSizeMeasure(std::string sM)
{

  if(sM == "average")
    sizeMeasure = 0; 
  else if (sM == "minimum")
    sizeMeasure = 1;
  else if (sM == "maximum")
    sizeMeasure = 2;
  else if (sM == "averagemax2")
    sizeMeasure = 3;
  else
    std::cout<<" Illegal size measure, "<<sM<<" specified."<<std::endl;

}

void meshRefine::setWriteJunctions(bool setJunctions)
{
  writeJunctions = setJunctions;
}

void meshRefine::readSurfaces(std::string filename)
{

  if(meshDimension == 2)
    {
      read2D(filename);
    }
  else if(meshDimension =3)
    {
      read3D(filename);
    }
  else
    std::cout<<"ERROR!!!: No mesh dimension specified."<<std::endl;
}

void meshRefine::read2D(std::string filename)
{

  std::ifstream surfFile(filename.c_str());

  std::cout<<" Starting to read surfaces "<<std::endl;

  if( !surfFile.is_open())
    {
      
      std::cout<<"Failed to open surfs file"<<std::endl;;
    }
  else
    std::cout<<"successfully opened surfs file."<<std::endl;
  
  //read in the surfs

  numSurfs=0;

  surfFile>>numSurfs;

  int numPoints;
  int surfType;
  double x,y,z;
  double nThick,pThick;

  surfs.resize(numSurfs);

  int surfCounter;

  for(surfCounter=0 ; surfCounter < numSurfs ; surfCounter++)
    {

      surfFile>>numPoints;
      surfFile>>surfType;
      surfFile>>nThick>>pThick;

      surfs[surfCounter].numNodes = numPoints;
      surfs[surfCounter].surfType = surfType;
      surfs[surfCounter].nThick = nThick;
      surfs[surfCounter].pThick = pThick;
      surfs[surfCounter].X.resize(numPoints);
      surfs[surfCounter].Y.resize(numPoints);
      surfs[surfCounter].Z.resize(numPoints);

      int node;
      double xc,yc,zc;
      xc = 0.0;
      yc = 0.0;
      zc = 0.0;

      for(node=0 ; node<numPoints ; node++)
	{

	  surfFile>>surfs[surfCounter].X[node]>>surfs[surfCounter].Y[node];

	  surfs[surfCounter].Z[node] = 0.0;

	  //fscanf(surfFile,"%lf %lf",&surfs[surfCounter].X[node],&surfs[surfCounter].Y[node]);

	  xc = xc + surfs[surfCounter].X[node];
	  yc = yc + surfs[surfCounter].Y[node];
	  zc = 0.0;

	}

      surfs[surfCounter].xc = xc/(double)numPoints;
      surfs[surfCounter].yc = yc/(double)numPoints;
      surfs[surfCounter].zc = zc/(double)numPoints;

      //Compute and store the normal to the surface
      //Note that the normal is set up by how the node are ordered in the linesegment

      surfs[surfCounter].normx =   surfs[surfCounter].X[1] - surfs[surfCounter].X[0];
      surfs[surfCounter].normy = -(surfs[surfCounter].Y[1] - surfs[surfCounter].Y[0]);
      surfs[surfCounter].normz =   0.0;

      double norm = sqrt(surfs[surfCounter].normx*surfs[surfCounter].normx + 
			 surfs[surfCounter].normy*surfs[surfCounter].normy + 
			 surfs[surfCounter].normz*surfs[surfCounter].normz);

      surfs[surfCounter].normx /= norm;
      surfs[surfCounter].normy /= norm;
      surfs[surfCounter].normz /= norm;

      surfs[surfCounter].size = norm;

    }

  surfFile.close();

}


void meshRefine::read3D(std::string filename)
{

  std::ifstream surfFile(filename.c_str());

  std::cout<<" Starting to read surfaces"<<std::endl;

  if( !surfFile.is_open())
    {

      std::cout<<"Failed to open surfs file surfaces.dat"<<std::endl;
    }
  else
    std::cout<<"successfully opened surfs file."<<std::endl;

  //read in the surfs

  numSurfs=0;

  surfFile>>numSurfs;

  int numPoints;
  int surfType;
  double x,y,z;
  double nThick,pThick;

  surfs.resize(numSurfs);

  int surfCounter;

  for(surfCounter=0 ; surfCounter < numSurfs ; surfCounter++)
    {

      surfFile>>numPoints;
      surfFile>>surfType;
      surfFile>>nThick>>pThick;


      surfs[surfCounter].numNodes = numPoints;
      surfs[surfCounter].surfType = surfType;
      surfs[surfCounter].nThick = nThick;
      surfs[surfCounter].pThick = pThick;
      surfs[surfCounter].X.resize(numPoints);
      surfs[surfCounter].Y.resize(numPoints);
      surfs[surfCounter].Z.resize(numPoints);

      int node;
      double xc,yc,zc;
      xc = 0.0;
      yc = 0.0;
      zc = 0.0;

      for(node=0 ; node<numPoints ; node++)
	{

	  surfFile>>surfs[surfCounter].X[node]>>surfs[surfCounter].Y[node]>>
	    surfs[surfCounter].Z[node];

	  xc = xc + surfs[surfCounter].X[node];
	  yc = yc + surfs[surfCounter].Y[node];
	  zc = zc + surfs[surfCounter].Z[node];
	}

      surfs[surfCounter].xc = xc/(double)numPoints;
      surfs[surfCounter].yc = yc/(double)numPoints;
      surfs[surfCounter].zc = zc/(double)numPoints;


      //Compute and store the normal to the surface

      double vecx[]={0,0},vecy[]={0,0},vecz[]={0,0};

      vecx[0] = surfs[surfCounter].X[1] - surfs[surfCounter].X[0];
      vecy[0] = surfs[surfCounter].Y[1] - surfs[surfCounter].Y[0];
      vecz[0] = surfs[surfCounter].Z[1] - surfs[surfCounter].Z[0];
      
      vecx[1] = surfs[surfCounter].X[2] - surfs[surfCounter].X[0];
      vecy[1] = surfs[surfCounter].Y[2] - surfs[surfCounter].Y[0];
      vecz[1] = surfs[surfCounter].Z[2] - surfs[surfCounter].Z[0];
  
      
      //Now create normal vector from cross product

      surfs[surfCounter].normx =   vecy[0]*vecz[1] - vecz[0]*vecy[1];
      surfs[surfCounter].normy = -(vecx[0]*vecz[1] - vecz[0]*vecx[1]);
      surfs[surfCounter].normz =   vecx[0]*vecy[1] - vecy[0]*vecx[1];

      double norm = sqrt(surfs[surfCounter].normx*surfs[surfCounter].normx + 
			 surfs[surfCounter].normy*surfs[surfCounter].normy + 
			 surfs[surfCounter].normz*surfs[surfCounter].normz);
      surfs[surfCounter].normx /= norm;
      surfs[surfCounter].normy /= norm;
      surfs[surfCounter].normz /= norm;

      int nodes[] = {0,1,2,0};
      double maxSide = 0;
      for(int i=0 ; i<3 ; ++i)
	{
	  double dist = distance(surfs[surfCounter].X[nodes[i+1]],surfs[surfCounter].Y[nodes[i+1]],
				 surfs[surfCounter].Z[nodes[i+1]],surfs[surfCounter].X[nodes[i]],
				 surfs[surfCounter].Y[nodes[i]], surfs[surfCounter].Z[nodes[i]]);

	  if(dist > maxSide)
	    maxSide = dist;

	  //maxSide += dist;
	}

      //maxSide /= 3.0;
      surfs[surfCounter].size = maxSide/2.0;

    }

  surfFile.close();

}


void meshRefine::printSurfs()
{

  int surfCounter;

  for(surfCounter=0 ; surfCounter < numSurfs ; surfCounter++)
    {

      std::cout<<"surf "<<surfCounter <<" has "<<surfs[surfCounter].numNodes<<" nodes."<<std::endl;

      int node;

      for(node=0 ; node<surfs[surfCounter].numNodes ; node++)
	{

	  std::cout<<"Node "<<node<<"   "<<surfs[surfCounter].X[node]<<"  "<<surfs[surfCounter].Y[node]
		   <<"  "<<surfs[surfCounter].Z[node]<<std::endl;
	}

      std::cout<<"Centroid at "<<surfs[surfCounter].xc<<"   "<<surfs[surfCounter].yc<<"    "
	       <<surfs[surfCounter].zc<<std::endl;

    }
}


/*  REMOVE THIS CODE
double * meshRefine::createDoublePtr(int size)
{

  double *x;
  x = malloc(size * sizeof(double));

  return x;
}
*/

void meshRefine::testX(std::vector<double> x)
{

  for(int i=0 ; i<10 ; ++i)
    x[i] = (double)i*i/(double)(i*i*i+1);

}


void meshRefine::printTestX(std::vector<double> x)
{

  for(int i=0 ; i<10 ; ++i)
    std::cout<<" x = "<<x[i]<<std::endl;

}

//----------------------------------------------------------------
//  set mionimum size of cell
//----------------------------------------------------------------

void meshRefine::setCellMinimum(double min)
{
  cellMinSize = min;
}


//----------------------------------------------------------------
//  set mionimum size of cell n side
//----------------------------------------------------------------

void meshRefine::setCellMinimumN(double min)
{
  cellMinSizeN = min;
  cellMinNSet = true;
}


//----------------------------------------------------------------
//  set mionimum size of cell p side
//----------------------------------------------------------------

void meshRefine::setCellMinimumP(double min)
{
  cellMinSizeP = min;
  cellMinPSet = true;
}


//----------------------------------------------------------------
//  set mionimum size of cell inversion layer
//----------------------------------------------------------------

void meshRefine::setCellMinimumIL(double min)
{
  cellMinSizeIL = min;
  cellMinILSet = true;
}


//----------------------------------------------------------------
//  Allocate memory for cell node location storage 
//----------------------------------------------------------------

void meshRefine::sizeCellNodes(int numCellNodes_)
{

  numCellNodes = numCellNodes_;

  cellNodeCoords.resize(numCellNodes);

  for(int i= 0 ; i<numCellNodes ; ++i)
    cellNodeCoords[i].resize(3);

  nodeCounter = 0;

}


//----------------------------------------------------------------
//  store cell node locations
//----------------------------------------------------------------

void meshRefine::fillCoordinates(std::vector<double> x_)
{

  int numDir;
  if(is3D)
    numDir=3;
  else
    {
      numDir=2;
      cellNodeCoords[nodeCounter][2] = 0.0;
    }

  for(int i=0 ; i<numDir ; i++)
    {
      cellNodeCoords[nodeCounter][i] = x_[i];
      //std::cout<<nodeCounter<<"    "<<i<<"    "<<x_[i]<<std::endl;
    }

  nodeCounter++;

}

//----------------------------------------------------------------
//  store cell node locations
//----------------------------------------------------------------

void meshRefine::resetNodeCounter()
{
  nodeCounter = 0;
}

void meshRefine::printCoords()
{

  for(int i=0 ; i<numCellNodes ; ++i)
    std::cout<<" Node "<<i<<" of "<<numCellNodes<<" coordinates are "<<
      cellNodeCoords[i][0]<<"   "<<cellNodeCoords[i][1]<<"    "<<cellNodeCoords[i][2]<<std::endl;
}


void meshRefine::freeCoords()
{

  for(int i=0 ; i<numCellNodes ; ++i)
    cellNodeCoords[i].clear();

  cellNodeCoords.clear();

  nodeCounter = 0;

}


void meshRefine::freeSurfs()
{

  for(int i=0 ; i<numSurfs ; ++i)
    {
      surfs[i].X.clear();
      surfs[i].Y.clear();
      surfs[i].Z.clear();
    }

  surfs.clear();

  numSurfs = 0;

}


bool meshRefine::doIRefine()
{

  if (cellNodeCoords[0][1] == 2)
    return true;
  else
    return false;


}

void meshRefine::setRefinementDistance(double dist_)
{
  refinementDistance = dist_;
}

void meshRefine::setRefinementDistanceN(double dist_)
{
  refinementDistanceN = dist_;
}

void meshRefine::setRefinementDistanceP(double dist_)
{
  refinementDistanceP = dist_;
}

void meshRefine::setRefinementDistanceIL(double dist_)
{
  refinementDistanceIL = dist_;
}


double meshRefine::getRefinementDistance()
{
  return refinementDistance;
}

int meshRefine::getNumSurfs()
{
  return numSurfs;
}

//----------------------------------------------------------------
//  base refinement on distance to an infinite y-z plane at xloc
//----------------------------------------------------------------


bool meshRefine::doIrefineXplane(double xloc)
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  double lengthMetric=getLengthMetric();


  if(!cellMinNSet)
    {
      std::cout<<"Must set n-side cell minimum size!"<<std::endl;
      exit(1);
    }

  if(!cellMinPSet)
    {
      std::cout  <<"Must set p-side cell minimum size!"<<std::endl;
      exit(1);
    }


  for(int i=0 ; i<numCellNodes ; ++i)
    {
      checkPSide=false;
      checkNSide=false;
      double rDist=0;

      if(cellNodeCoords[i][0] < xloc)
	rDist = refinementDistanceP;

      if(cellNodeCoords[i][0] >= xloc)
	rDist = refinementDistanceN;

      double rSide=0;
      if(cellNodeCoords[i][0] < xloc)
	{
	  rSide = cellMinSizeP;
	  checkNSide = true;
	}
      if(cellNodeCoords[i][0] >= xloc)
	{
	  rSide = cellMinSizeN;
	  checkPSide=true;
	}

      if(lengthMetric <= rSide)
	{
	  refine = false;
	}
      else
	if(fabs(cellNodeCoords[i][0] - xloc) < rDist)
	  {
	    refine = true;
	    if(checkPSide)
	      pRefine = true;
	    if(checkNSide)
	      nRefine = true;
	  }
    }


  if(pRefine || nRefine)
    {
      return true;
    }
  else
    {
      return false;
    }

}


//-------------------------------------------------------------------------------------
//  base refinement on distance to the centroid of the set of surfaces in memory
//-------------------------------------------------------------------------------------


bool meshRefine::doIrefineCentroid()
{

  bool refine=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  std::vector<double> point;

  point.resize(3);

  //  if(maxSide < cellMinSize)
  if(minSide < cellMinSize)
    {
      return refine;
    }

  for(int surf=0 ; surf < numSurfs ; surf++)
    {
      for(int i=0 ; i<numCellNodes ; ++i)
	{
	  if(distance(cellNodeCoords[i][0],cellNodeCoords[i][1],cellNodeCoords[i][2],
		      surfs[surf].xc,surfs[surf].yc,surfs[surf].zc)  < refinementDistance)
	    {
	      refine = true;
	      break;
	    }
	}
    }

  //In cases where the tet is large compared to the refinement distance, the surface may lie
  //inside the tet without any node within the refinement distance.  Need to check it.

  for(int surf=0 ; surf < numSurfs ; surf++)
    {
      //Point Check

      point[0] = surfs[surf].xc;;
      point[1] = surfs[surf].yc;
      point[2] = surfs[surf].zc;
      
      bool check1 = signedDistanceSameSideCheck(point,cellNodeCoords[0],cellNodeCoords[1],
						cellNodeCoords[2],cellNodeCoords[3]);
      bool check2 = signedDistanceSameSideCheck(point,cellNodeCoords[1],cellNodeCoords[0],
						cellNodeCoords[2],cellNodeCoords[3]);
      bool check3 = signedDistanceSameSideCheck(point,cellNodeCoords[2],cellNodeCoords[1],
						cellNodeCoords[0],cellNodeCoords[3]);
      bool check4 = signedDistanceSameSideCheck(point,cellNodeCoords[3],cellNodeCoords[1],
						cellNodeCoords[2],cellNodeCoords[0]);
      
      if(check1 && check2 && check3 && check4)
	{
	  refine = true;
	  break;
	}
    }

  return refine;

}


//-------------------------------------------------------------------------------------
//  base 3D refinement on distance to the surfaces in memory
//-------------------------------------------------------------------------------------


bool meshRefine::doIrefine3D()
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool iLRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  bool checkILSide=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  double aveSide = getAveSide(true);
  std::vector<double> point,localCoords;
  double cellSize;

  if (sizeMeasure==0)
    cellSize = aveSide;
  else if (sizeMeasure==1)
    cellSize = minSide;
  else if (sizeMeasure==2)
    cellSize = maxSide;
  else
    cellSize = aveSide;
  if(cellSize < cellMinSize)  //This has met size requirement
    {
      //std::cout<<" Met size requirement "<<cellSize<<"    "<<cellMinSize<<std::endl;
      //if(tetNum == 3106)
      //std::cout<<" Returning false for "<<tetNum<<std::endl;
      return false;
    }

  //Determine if this is a p cell or an n cell
  std::vector<double> centroidX(3,0.0);
  for (size_t i=0 ; i<cellNodeCoords.size() ; ++i)
    {
      centroidX[0] += cellNodeCoords[i][0];
      centroidX[1] += cellNodeCoords[i][1];
      centroidX[2] += cellNodeCoords[i][2];
    }
  centroidX[0] /= (double)cellNodeCoords.size();
  centroidX[1] /= (double)cellNodeCoords.size();
  centroidX[2] /= (double)cellNodeCoords.size();

  double centroidDoping = function.evaluateFunction(centroidX);

  if(centroidDoping >= 0) //ntype
    semiType = 1;
  if(centroidDoping < 0) // ptype
    semiType = 2;

  //Set up data structures
  point.resize(3);
  localCoords.resize(3);
  distanceFunctions dFun;


  //In cases where the tet is large compared to the refinement distance, the surface may lie
  //inside the tet without any node within the refinement distance.  Need to check it.
  //If the cell is greater than the minimum size, always refine it if it is transected by a junction

  for(int surf=0 ; surf < numSurfs ; surf++)
    {
      //Point Check
      
      surfaceInfo localSurf = surfs[surf];

      point[0] = localSurf.xc;
      point[1] = localSurf.yc;
      point[2] = localSurf.zc;
      

      double f0 = function.evaluateFunction(cellNodeCoords[0]);
      double f1 = function.evaluateFunction(cellNodeCoords[1]);
      double f2 = function.evaluateFunction(cellNodeCoords[2]);
      double f3 = function.evaluateFunction(cellNodeCoords[3]);
      if( f0*f1 < 0 or f0*f2 < 0 or f0*f3 < 0 or
	  f1*f2 < 0 or f1*f3 < 0 or f2*f3 < 0)
	{
	  //intersection 
	  return true;
	}

      //Set the junction refinement distance
      double junctionRefinement=0.0;
      junctionRefinement = refinementDistance;
      if (autoRefine)
	{
	  if (semiType == 1)
	    junctionRefinement = localSurf.nThick;
	  if (semiType == 2)
	    junctionRefinement = localSurf.pThick;
	}

      double xnew,ynew,znew;

      for(int i=0 ; i<numCellNodes ; ++i)
	{
	  localCoords[0] = cellNodeCoords[i][0]; 
	  localCoords[1] = cellNodeCoords[i][1]; 
	  localCoords[2] = cellNodeCoords[i][2];
	  if(localSurf.surfType == 2)  //Refine to a specified surface
	    {
	      //Get a true normal distance for this one
	      double distance = dFun.normalDistanceToSurface(localCoords[0],localCoords[1],localCoords[2],
							  xnew,ynew,znew,localSurf);

	      if(distance <= localSurf.refinementDistance)
		{
		  return true;
		}

	      if(distance <= 2*localSurf.refinementDistance and cellSize > refinementFactor*cellMinSize)
		{
		  return true;
		}

	      if(distance <= 3*localSurf.refinementDistance and cellSize > refinementFactor*refinementFactor*cellMinSize)
		{
		  return true;
		}
	    }

	  if(localSurf.surfType == 0)  //Refine to a junction
	    {
	      //Get a normal distance plus ends for this one
	      double distance = dFun.distanceToSurface(localCoords[0],localCoords[1],localCoords[2],
						       xnew,ynew,znew,localSurf);

	      bool localReturn=false;


	      if(distance <= junctionRefinement)
		{
		  localReturn = true;
		}

	      if(distance <= 2*junctionRefinement and cellSize > refinementFactor*cellMinSize)
		{
		  localReturn = true;
		}

	      if(distance <= 3*junctionRefinement and cellSize > refinementFactor*refinementFactor*cellMinSize)
		{
		  localReturn = true;
		}
	      if(localReturn)
		{
		  return true;
		}

	    }
	}
    }


  return refine;

}


//-------------------------------------------------------------------------------------
//  base refinement on distance to the centroid of the set of surfaces in memory
//-------------------------------------------------------------------------------------


bool meshRefine::doIrefineCentroidSided()
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool iLRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  bool checkILSide=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  double aveSide = getAveSide(true);
  std::vector<double> point;

  point.resize(3);

  //find a subset of surfaces to check by centroid
  std::vector<int> subSurfList;
  int numSubs=10;

  subSurfList.resize(numSubs);

  findClosest(numSubs,subSurfList);

  flagTet=false;

  /*
    if( fabs(cxc-1700) < 75 &&
    fabs(cyc-0) < 35 &&
    fabs(czc-300) < 35)
    flagTet = true;
  */
  

  if(!cellMinNSet)
    {
      std::cout  <<"Must set n-side cell minimum size!"<<std::endl;
      exit(1);
    }
  
  if(!cellMinPSet)
    {
      std::cout  <<"Must set p-side cell minimum size!"<<std::endl;
      exit(1);
    }
  
  double localCoords[] = {cxc,cyc,czc};

  //for(int surf=0 ; surf < numSurfs ; surf++)
  for(int surf=0 ; surf < numSubs ; surf++)
    {
      for(int i=0 ; i<numCellNodes ; ++i)
	{
	  
	  localCoords[0] = cellNodeCoords[i][0]; 
	  localCoords[1] = cellNodeCoords[i][1]; 
	  localCoords[2] = cellNodeCoords[i][2]; 
	  
	  checkPSide=false;
	  checkNSide=false;
	  checkILSide=false;
	  surfaceInfo localSurf = surfs[subSurfList[surf]];
	  //DEBUG std::cout<<i<<"  Local surf "<<subSurfList[surf]<<"    "<<surf<<std::endl;  //DEBUG 
	  //Determine which type of surf we're looking at
	  //	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],surfs[surf]);
	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],localSurf);
	  double rDist = 0;
	  double rSize=0;
	  if(type == 0)  //p type
	    {
	      rDist = refinementDistanceP;
	      rSize = cellMinSizeP;
	      rDist = localSurf.pThick;
	      rSize = localSurf.pThick/10;
	      checkPSide=true;
	    }
	  if(type == 1)  //n type
	    {
	      rDist = refinementDistanceN;
	      rSize = cellMinSizeN;
	      rDist = localSurf.nThick;
	      rSize = localSurf.nThick/10;
	      checkNSide=true;
	    }
	  if(type == 2)  //inversion layer
	    {
	      rDist = refinementDistanceIL;
	      rSize = cellMinSizeIL;
	      checkILSide=true;
	    }
	  
	  
	  if(aveSide < rSize)
	    return false;

	  double angle = 0;
	  //This is not elegant, but if it's N check, flip the normal
	  if(checkPSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
				      localSurf.normy,localSurf.normz);  //cosine is even function
	  if(checkNSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
				      -localSurf.normy,-localSurf.normz);  //cosine is even function
	  if(checkILSide)
	    {
	      angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
					localSurf.normy,localSurf.normz);  //cosine is even function
	      if(cos(angle) < 0)
		angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					  localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
					  -localSurf.normy,-localSurf.normz);  //cosine is even function
	    }
	  

	  if(angle<0)
	    std::cout<<" SOMETHING HAS GONE WRONG"<<std::endl;
	  
	  //double checkDistance = rDist*cos(angle) + localSurf.size*sin(angle);
	  
	  double checkDistance = sqrt((localCoords[0]-localSurf.xc)*(localCoords[0]-localSurf.xc) +
				      (localCoords[1]-localSurf.yc)*(localCoords[1]-localSurf.yc) +
				      (localCoords[2]-localSurf.zc)*(localCoords[2]-localSurf.zc))*cos(angle);
	  
	  bool refineOrNot=false;

	  if(checkDistance <= rDist)
	    refineOrNot = true;
	  
	  //Special consideration for oxide adjacent cells
	  if(checkILSide && cos(angle) < 1e-10 )
	    {
	      if(distance(localCoords[0],localCoords[1],localCoords[2],
			  localSurf.xc,localSurf.yc,localSurf.zc) < localSurf.size)
		{
		  // This node sits on the oxide;
		  //std::cout<<"Checking against inversion layer surface,  %f   %f \n",angle,cos(angle));
		  refineOrNot = true;
		  checkDistance = 0;
		}
	    }
	  
	  //If it's close enough, check if it's in the cylinder
	  if(refineOrNot)
	    {
	      double xcomp = localCoords[0] - (localSurf.xc + localSurf.normx*checkDistance);
	      double ycomp = localCoords[1] - (localSurf.yc + localSurf.normy*checkDistance);
	      double zcomp = localCoords[2] - (localSurf.zc + localSurf.normz*checkDistance);
	      
	      double cylinder = sqrt(xcomp*xcomp + ycomp*ycomp + zcomp*zcomp);
	      
	      /*
		if(checkILSide)
		std::cout<<"Checking against IL %d  %f  %f  \n",surf,cylinder,localSurf.size);
	      */
	      
	      //std::cout<<"cylinder check "<<localSurf.size<<"   "<<cylinder<<std::endl;
	      //std::cout<<"cylinder check "<< localCoords[0]<<"    "<<localCoords[1]<<"    "<<localCoords[2]<<std::endl;


	      //Outside the cylinder; shut down refinement
	      if(cylinder > localSurf.size)
		refineOrNot = false;
	    }
	  
	  //On the inversion layer, just use the straight distance
	  if(checkILSide)checkDistance = rDist;
	  
	  //if(distance(localCoords[0],localCoords[1],localCoords[2],
	  //	      localSurf.xc,localSurf.yc,localSurf.zc)  < rDist)
	  //if(distance(localCoords[0],localCoords[1],localCoords[2],
	  //	      localSurf.xc,localSurf.yc,localSurf.zc)  < checkDistance)
	  if(refineOrNot)
	    {
	      
	      refine = true;
	      if(checkPSide)
		pRefine = true;
	      if(checkNSide)
		nRefine = true;
	      if(checkILSide)
		iLRefine = true;
	      
	      if(flagTet)
		{
		  std::cout  <<" Made a positive check at local surf "<<surf<<"  distance "<<checkDistance
			     <<" angle "<<angle<<std::endl;
		  std::cout  <<"    Node located at: "<<localCoords[0]<<"    "
			     <<localCoords[1]<<"    "<<localCoords[2]<<std::endl;
		  std::cout  <<" surface located at: "<<localSurf.xc<<"    "<<localSurf.yc<<"    "
			     <<localSurf.zc<<std::endl;
		  std::cout  <<"          with norm: "<<localSurf.normx<<"    "<<localSurf.normy<<"    "
			     <<localSurf.normz<<std::endl;
		  
		  std::cout<<" checkComponents: "<<rDist<<"    "<<cos(angle)<<"    "<<localSurf.size<<"    "
			   <<sin(angle)<<std::endl;
		  if(checkNSide)
		    std::cout<< "N Check "<<distance(localCoords[0],localCoords[1],localCoords[2],
						     localSurf.xc,localSurf.yc,localSurf.zc)
			     <<"    "<<rDist
			     <<"    "<<aveSide<<"    "<<rSize<<std::endl;
		  
		  if(checkPSide)
		    std::cout<<" P  check is true"<<std::endl<<std::endl;
		  if(checkNSide)
		    std::cout<<" N  check is true"<<std::endl<<std::endl;
		  if(checkILSide)
		    std::cout<<" IL check is true " <<std::endl<<std::endl;
		}
	    }
	}
    }
  
  //In cases where the tet is large compared to the refinement distance, the surface may lie
  //inside the tet without any node within the refinement distance.  Need to check it.
  
  //If the decision to refine has already been made, skip this 
  if(!refine)
    {
      //for(int surf=0 ; surf < numSurfs ; surf++)
      for(int surf=0 ; surf < numSubs ; surf++)
	{
	  //Point Check

	  surfaceInfo localSurf = surfs[subSurfList[surf]];

	  point[0] = localSurf.xc;;
	  point[1] = localSurf.yc;
	  point[2] = localSurf.zc;
      
	  bool check1 = signedDistanceSameSideCheck(point,cellNodeCoords[0],cellNodeCoords[1],
						    cellNodeCoords[2],cellNodeCoords[3]);
	  bool check2 = signedDistanceSameSideCheck(point,cellNodeCoords[1],cellNodeCoords[0],
						    cellNodeCoords[2],cellNodeCoords[3]);
	  bool check3 = signedDistanceSameSideCheck(point,cellNodeCoords[2],cellNodeCoords[1],
						    cellNodeCoords[0],cellNodeCoords[3]);
	  bool check4 = signedDistanceSameSideCheck(point,cellNodeCoords[3],cellNodeCoords[1],
						    cellNodeCoords[2],cellNodeCoords[0]);
      
	  double sizeCheck = 0;
	  if(cellMinSizeN < cellMinSizeP)
	    sizeCheck = cellMinSizeN;
	  else
	    sizeCheck = cellMinSizeP;

	  if(check1 && check2 && check3 && check4 && (aveSide > sizeCheck))
	    {
	      refine = true;
	    }
	  if(refine)break;
	}
    }

  subSurfList.clear();

  return refine;

}



//-------------------------------------------------------------------------------------
//  base refinement on distance to the set of surfaces in memory
//-------------------------------------------------------------------------------------

bool meshRefine::doIrefine2D()
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool iLRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  bool checkILSide=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  double aveSide = getAveSide(true);
  double norm[3];

  //Determine if this is a p cell or an n cell
  std::vector<double> centroidX(3,0.0);
  for (size_t i=0 ; i<cellNodeCoords.size() ; ++i)
    {
      centroidX[0] += cellNodeCoords[i][0];
      centroidX[1] += cellNodeCoords[i][1];
    }
  centroidX[0] /= (double)cellNodeCoords.size();
  centroidX[1] /= (double)cellNodeCoords.size();

  double centroidDoping = function.evaluateFunction(centroidX);

  if(centroidDoping >= 0) //ntype
    semiType = 1;
  if(centroidDoping < 0) // ptype
    semiType = 2;

  //If this cell has hit mimimum size threshold, return false
  bool hitMinimumSize=false;
  //if(aveSide < cellMinSize)
  if(maxSide < cellMinSize)
    hitMinimumSize = true;
  if(hitMinimumSize)return false;
	  
  std::vector<double> localCoords;
  localCoords.resize(2);

  distanceFunctions dFun;

  for(size_t surf=0 ; surf < surfs.size() ; surf++)
    {
      //First check if the surface transects or lies inside the cell
      //If it does select it for refinement provided cell is larger than minimum size
      surfaceInfo localSurf = surfs[surf];
      int nodes[4]={0,1,2,0};

      double axmin,axmax,aymin,aymax;

      if(localSurf.X[0] <= localSurf.X[1])
	{
	  axmin = localSurf.X[0];
	  axmax = localSurf.X[1];
	}
      else
	{
	  axmin = localSurf.X[1];
	  axmax = localSurf.X[0];
	}
      if(localSurf.Y[0] <= localSurf.Y[1])
	{
	  aymin = localSurf.Y[0];
	  aymax = localSurf.Y[1];
	}
      else
	{
	  aymin = localSurf.Y[1];
	  aymax = localSurf.Y[0];
	}

      double sizeax = axmax-axmin;
      double sizeay = aymax-aymin;
      double posax  = axmin;
      double posay  = aymin;
  
      for(int curve=0 ; curve < 3 ; ++curve)
	{
	  double bxmin,bxmax,bymin,bymax;

	  if(cellNodeCoords[nodes[curve]][0] <= cellNodeCoords[nodes[curve+1]][0])
	    {
	      bxmin = cellNodeCoords[nodes[curve]][0];
	      bxmax = cellNodeCoords[nodes[curve+1]][0];
	    }
	  else
	    {
	      bxmin = cellNodeCoords[nodes[curve+1]][0];
	      bxmax = cellNodeCoords[nodes[curve]][0];
	    }
	  if(cellNodeCoords[nodes[curve]][1] <= cellNodeCoords[nodes[curve+1]][1])
	    {
	      bymin = cellNodeCoords[nodes[curve]][1];
	      bymax = cellNodeCoords[nodes[curve+1]][1];
	    }
	  else
	    {
	      bymin = cellNodeCoords[nodes[curve+1]][1];
	      bymax = cellNodeCoords[nodes[curve]][1];
	    }

	  double sizebx = bxmax-bxmin;
	  double sizeby = bymax-bymin;
	  double posbx  = bxmin;
	  double posby  = bymin;
  
	  if ((posax < posbx + sizebx) &&
	      (posay < posby + sizeby) &&
	      (posbx < posax + sizeax) &&
	      (posby < posay + sizeay))
	    {
	      return true;
	    }
	}
      
      //Set the junction refinement distance
      double junctionRefinement=0.0;
      junctionRefinement = refinementDistance;
      if (autoRefine)
	{
	  if (semiType == 1)
	    junctionRefinement = localSurf.nThick;
	  if (semiType == 2)
	    junctionRefinement = localSurf.pThick;
	}

      double xnew,ynew;
      std::vector<double> surfXCoords=localSurf.X;
      std::vector<double> surfYCoords=localSurf.Y;
      for(int i=0 ; i<numCellNodes ; ++i)
	{
	  localCoords[0] = cellNodeCoords[i][0]; 
	  localCoords[1] = cellNodeCoords[i][1];
	  double distance = 0; 
	  switch(localSurf.surfType)
	    {

	    case 0:
	      //Get a normal distance plus ends for this one
	      //double distance = dFun.normalDistanceToLine(localCoords[0],localCoords[1],
	      //					  xnew,ynew,localSurf.X,localSurf.Y);
	      distance = dFun.distanceToLine(localCoords[0],localCoords[1],
					     xnew,ynew,localSurf.X,localSurf.Y);

	      //Check for refinement--grade away to 3x distance increasing minimum distance by 10x each time
	      if(distance <= junctionRefinement and maxSide > cellMinSizeN)
		{
		  return true;
		}

	      if(distance <= 2*junctionRefinement and maxSide > refinementFactor*cellMinSizeN)
		{
		  return true;
		}

	      if(distance <= 3*junctionRefinement and maxSide > refinementFactor*refinementFactor*cellMinSizeN)
		{
		  return true;
		}

	      break;

	    case 2:

	      //Get a true normal distance for this one
	      distance = dFun.normalDistanceToLine(localCoords[0],localCoords[1],
						   xnew,ynew,localSurf.X,localSurf.Y);
	      if(distance <= localSurf.refinementDistance and maxSide > cellMinSizeIL)
		{
		  return true;
		}

	      if(distance <= 2*localSurf.refinementDistance and maxSide > refinementFactor*cellMinSizeIL)
		{
		  return true;
		}

	      if(distance <= 3*localSurf.refinementDistance and maxSide > refinementFactor*refinementFactor*cellMinSizeIL)
		{
		  return true;
		}

	      break;

	    default:
	      std::cout<<" Unknown surfType found "<<localSurf.surfType<<std::endl;
	      exit;

	    }
	}
    }

  return false;

}
//-------------------------------------------------------------------------------------
//  base refinement on distance to the centroid of the set of surfaces in memory
//-------------------------------------------------------------------------------------


bool meshRefine::doIrefine2DDepricated()
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool iLRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  bool checkILSide=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  double aveSide = getAveSide(true);
  std::vector<double> point;
  double norm[3];

  point.resize(3);

  //find a subset of surfaces to check by centroid
  std::vector<int> subSurfList;
  int numSubs=10;

  numSubs = numSurfs;

  subSurfList.resize(numSubs);

  for(int i=0 ; i<numSubs ; ++i)
    subSurfList[i] = i;

  //findClosest(numSubs,subSurfList);

  flagTet=true;

  /*
  for(int i=0 ; i<3 ; ++i)
    {
      if(cellNodeCoords[i][0] > 1.0 && cellNodeCoords[i][0] < 2.0)
	if(cellNodeCoords[i][1] > -0.2)
	  {
	    flagTet = true;
	  }
    }
  */

  if(!cellMinNSet)
    {
      std::cout  <<"Must set n-side cell minimum size!"<<std::endl;
      exit(1);
    }

  if(!cellMinPSet)
    {
      std::cout  <<"Must set p-side cell minimum size!"<<std::endl;
      exit(1);
    }

  double localCoords[] = {cxc,cyc,czc};

  //for(int surf=0 ; surf < numSurfs ; surf++)
  for(int surf=0 ; surf < numSubs ; surf++)
    {
      //std::cout<<" SURFS "<<surf<<"    "<<surfs[surf].surfType<<std::endl;
      //if(surf < 20)
	{
	  //std::cout<<" SURFS "<<surf<<"     "<<surfs[surf].size<<std::endl;
	  //std::cout<<" SURFS "<<surf<<"     "<<surfs[surf].pThick<<std::endl;
	  //std::cout<<" SURFS "<<surf<<"    "<<surfs[surf].X[0]<<"     "<<surfs[surf].X[1]<<std::endl;
	  //std::cout<<" SURFS "<<surf<<"    "<<surfs[surf].Y[0]<<"     "<<surfs[surf].Y[1]<<std::endl;
	}
      for(int i=0 ; i<numCellNodes ; ++i)
	{

	  localCoords[0] = cellNodeCoords[i][0]; 
	  localCoords[1] = cellNodeCoords[i][1]; 
	  localCoords[2] = cellNodeCoords[i][2]; 

	  if(fabs(localCoords[2]) > 1.e-8)
	    std::cout<<" bad Z coordinate "<<localCoords[2]<<std::endl;

	  checkPSide=false;
	  checkNSide=false;
	  checkILSide=false;
	  surfaceInfo localSurf = surfs[subSurfList[surf]];
	  //Determine which type of surf we're looking at
	  //	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],surfs[surf]);

	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],localSurf);
	  double rDist = 0;
	  double rSize=0;
	  if(type == 0)  //p type
	    {
	      //rDist = refinementDistanceP;
	      //rSize = cellMinSizeP;
	      rDist = localSurf.pThick;
	      rSize = localSurf.pThick/10;
	      checkPSide=true;
	    }
	  if(type == 1)  //n type
	    {
	      //rDist = refinementDistanceN;
	      //rSize = cellMinSizeN;
	      rDist = localSurf.nThick;
	      rSize = localSurf.nThick/10;
	      checkNSide=true;
	    }
	  if(type == 2)  //inversion layer
	    {
	      rDist = refinementDistanceIL;
	      rSize = cellMinSizeIL;
	      checkILSide=true;
	    }


	  double angle = 0;
	  //This is not elegant, but if it's N check, flip the normal
	  if(checkPSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
				      localSurf.normy,localSurf.normz);  //cosine is even function
	  if(checkNSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
				      -localSurf.normy,-localSurf.normz);  //cosine is even function
	  if(checkILSide)
	    {
	      angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
					localSurf.normy,localSurf.normz);  //cosine is even function
	      if(cos(angle) < 0)
		angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					  localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
					  -localSurf.normy,-localSurf.normz);  //cosine is even function
	    }

	  //if(angle<0)
	  if(angle>1.571)
	    std::cout<<" SOMETHING HAS GONE WRONG"<<std::endl;

	  double checkDistance = sqrt((localCoords[0]-localSurf.xc)*(localCoords[0]-localSurf.xc) +
				      (localCoords[1]-localSurf.yc)*(localCoords[1]-localSurf.yc) +
				      (localCoords[2]-localSurf.zc)*(localCoords[2]-localSurf.zc))*cos(angle);


	  double checkDistanceSin = sqrt((localCoords[0]-localSurf.xc)*(localCoords[0]-localSurf.xc) +
					 (localCoords[1]-localSurf.yc)*(localCoords[1]-localSurf.yc) +
					 (localCoords[2]-localSurf.zc)*(localCoords[2]-localSurf.zc))*sin(angle);



	  bool refineOrNot=false;

	  //std::cout<<" checkDistance = %f %f %f \n",checkDistance,rDist,cos(angle));

	  if(aveSide<rSize)
	    continue;
	  // return false;

	  if(checkDistanceSin <= rDist && checkDistance < localSurf.size/2 && localCoords[0]<1.0 )
	    std::cout<<"setting true based on size "<<localCoords[0]<<std::endl;

	  if(checkDistanceSin <= rDist && checkDistance < localSurf.size/2)
	    return true;

	  //Special consideration for oxide adjacent cells
	  if(checkILSide && cos(angle) < 1e-10 )
	    {
	      if(distance(localCoords[0],localCoords[1],localCoords[2],
			  localSurf.xc,localSurf.yc,localSurf.zc) < localSurf.size)
		{
		  // This node sits on the oxide;
		  //std::cout<<"Checking against inversion layer surface,  %f   %f \n",angle,cos(angle));
		  refineOrNot = true;
		  checkDistance = 0;
		  if(localCoords[0]<1.0 )
		    std::cout<<"Not even near oxide "<<localCoords[0]<<std::endl;
		}
	    }

	  //If it's close enough, check if it's in the cylinder
	  if(refineOrNot)
	    {
	      double xcomp = localCoords[0] - (localSurf.xc + localSurf.normx*checkDistance);
	      double ycomp = localCoords[1] - (localSurf.yc + localSurf.normy*checkDistance);
	      double zcomp = localCoords[2] - (localSurf.zc + localSurf.normz*checkDistance);

	      double cylinder = sqrt(xcomp*xcomp + ycomp*ycomp + zcomp*zcomp);
	      
	      /*
		if(checkILSide)
		std::cout<<"Checking against IL %d  %f  %f  \n",surf,cylinder,localSurf.size);
	      */

	      //std::cout<<" checking cylinder.  %f  %f  \n",cylinder,localSurf.size);

	      if(surf<20)
		std::cout<<" SURFCHECK "<<cylinder<<"    "<<localSurf.size<<"    "<<localCoords[0]<<std::endl;

	      //Outside the cylinder; shut down refinement
	      if(cylinder > localSurf.size)
		refineOrNot = false;
	    }
	  
	  //On the inversion layer, just use the straight distance
	  if(checkILSide)checkDistance = rDist;

	  if(refineOrNot)
	    {
	      
	      refine = true;
	      if(checkPSide)
		pRefine = true;
	      if(checkNSide)
		nRefine = true;
	      if(checkILSide)
		iLRefine = true;

	      if(flagTet)
		{
		  std::cout  <<" Made a positive check at local surf "<<surf<<"  distance "<<checkDistance
			     <<" angle "<<angle<<std::endl;
		  std::cout  <<"    Node located at: "<<localCoords[0]<<"    "
			     <<localCoords[1]<<"    "<<localCoords[2]<<std::endl;
		  std::cout  <<" surface located at: "<<localSurf.xc<<"    "<<localSurf.yc<<"    "
			     <<localSurf.zc<<std::endl;
		  std::cout  <<"          with norm: "<<localSurf.normx<<"    "<<localSurf.normy
			     <<"    "<<localSurf.normz<<std::endl;
		  
		  std::cout<<" checkComponents: "<<rDist<<"    "<<cos(angle)<<"    "<<localSurf.size
			   <<"    "<<sin(angle)<<std::endl;
		  if(checkNSide)
		    std::cout<< "N Check "
			     <<distance(localCoords[0],localCoords[1],localCoords[2],
					localSurf.xc,localSurf.yc,localSurf.zc)
			     <<"    "<<rDist<<"    "<<aveSide<<"    "<<rSize<<std::endl;;
		  
		  if(checkPSide)
		    std::cout<<" P  check is true "<<std::endl<<std::endl;
		  if(checkNSide)
		    std::cout<<" N  check is true "<<std::endl;
		  if(checkILSide)
		    std::cout<<" IL check is true "<<std::endl<<std::endl;
		}
	    }
	}
    }
  
  //In cases where the tet is large compared to the refinement distance, the surface may lie
  //inside the tet without any node within the refinement distance.  Need to check it.

  //If the decision to refine has already been made, skip this
  //Revisit this for 2D

  /*
    if(flagTet)
    {
    if(refine) std::cout<<" REFINE IT\n");
    if(!refine)std::cout<<"DON'T REFINE \n");
    }
  */
  
  if(!refine)
    {
      //for(int surf=0 ; surf < numSurfs ; surf++)
      for(int surf=0 ; surf < numSubs ; surf++)
	{
	  //Point Check
	  
	  surfaceInfo localSurf = surfs[subSurfList[surf]];
	  
	  int nodes[4]={0,1,2,0};
	  
	  for(int curve=0 ; curve < 3 ; ++curve)
	    {

	      double axmin;
	      double axmax;
	      double aymin;
	      double aymax;

	      if(localSurf.X[0] <= localSurf.X[1])
		{
		  axmin = localSurf.X[0];
		  axmax = localSurf.X[1];
		}
	      else
		{
		  axmin = localSurf.X[1];
		  axmax = localSurf.X[0];
		}
	      if(localSurf.Y[0] <= localSurf.Y[1])
		{
		  aymin = localSurf.Y[0];
		  aymax = localSurf.Y[1];
		}
	      else
		{
		  aymin = localSurf.Y[1];
		  aymax = localSurf.Y[0];
		}

	      double bxmin;
	      double bxmax;
	      double bymin;
	      double bymax;


	      if(cellNodeCoords[nodes[curve]][0] <= cellNodeCoords[nodes[curve+1]][0])
		{
		  bxmin = cellNodeCoords[nodes[curve]][0];
		  bxmax = cellNodeCoords[nodes[curve+1]][0];
		}
	      else
		{
		  bxmin = cellNodeCoords[nodes[curve+1]][0];
		  bxmax = cellNodeCoords[nodes[curve]][0];
		}
	      if(cellNodeCoords[nodes[curve]][1] <= cellNodeCoords[nodes[curve+1]][1])
		{
		  bymin = cellNodeCoords[nodes[curve]][1];
		  bymax = cellNodeCoords[nodes[curve+1]][1];
		}
	      else
		{
		  bymin = cellNodeCoords[nodes[curve+1]][1];
		  bymax = cellNodeCoords[nodes[curve]][1];
		}

	      if( (axmin <= bxmax)  &&
		  (axmax >= bxmin)  &&
	          (aymin <= bymax)  &&
		  (aymax >= bymin))
		{
		  //std::cout<<" INTERSECTION "<<std::endl;
		  return true;
		}

	    }


	  /*

	  point[0] = localSurf.xc;;
	  point[1] = localSurf.yc;
	  point[2] = localSurf.zc;
	  norm[0] = localSurf.normx;
	  norm[1] = localSurf.normy;
	  norm[2] = 0.0;

	  double vec1[2],vec2[2],vec3[2];

	  vec1[0] = cellNodeCoords[0][0] - point[0];
	  vec1[1] = cellNodeCoords[0][1] - point[1];
       
	  vec2[0] = cellNodeCoords[1][0] - point[0];
	  vec2[1] = cellNodeCoords[1][1] - point[1];
       
	  vec3[0] = cellNodeCoords[2][0] - point[0];
	  vec3[1] = cellNodeCoords[2][1] - point[1];
       
	  int check1,check2,check3;

	  double DP1 = vec1[0]*vec2[0] + vec1[1]*vec2[1];
	  double DP2 = vec1[0]*vec3[0] + vec1[1]*vec3[1];
	  double DP3 = vec3[0]*vec2[0] + vec3[1]*vec2[1];

	  if(DP1 == 0 || DP2 == 0 || DP3 == 0)
	    return true;


	  if(DP1 > 0 )
	    check1 = 1;
	  else
	    check1 = -1;

	  if(DP2 > 0 )
	    check2 = 1;
	  else
	    check2 = -1;

	  if(DP3 > 0 )
	    check3 = 1;
	  else
	    check3 = -1;

	  if(DP1==0.0 || DP2==0.0 || DP3==0.0)
	    {
	      return true;
	    }

	  if(cellNodeCoords[0][0] < 1)
	    flagTet = true;

	  if(flagTet)
	    {
	      std::cout<<"\n POINT      %f %f \n ",point[0],point[1]);
	      std::cout<<" FLAG  vec1 %f %f  %f\n",vec1[0],vec1[1],DP1);
	      std::cout<<" FLAG  vec2 %f %f  %f\n",vec2[0],vec2[1],DP2);
	      std::cout<<" FLAG  vec2 %f %f  %f\n",vec3[0],vec3[1],DP3);
	      std::cout<<"\n");
	    }

	  //CheckSide returns either 1 or -1.  If not ALL returns are the same, the surface
	  //must lie inside the traingle.  If so, refine it.

	  double sizeCheck = 0;
	  if(cellMinSizeN < cellMinSizeP)
	    sizeCheck = cellMinSizeN;
	  else
	    sizeCheck = cellMinSizeP;

	  if(  (check1 > 0 && check2 > 0 && check3 > 0)  ||
	       (check1 < 0 && check2 < 0 && check3 < 0))
	    refine = false;
	  else
	    if(aveSide > sizeCheck)
	      {
		refine = true;
	      }

	  if(flagTet)
	    {
	    if(refine)
	      {
		std::cout<<"refine TRUE %f %f %f  %d  %d  %d\n",point[0],point[1],point[2],check1,check2,check3);
		std::cout<<"Coords 1  %f %f %f \n",cellNodeCoords[0][0],cellNodeCoords[0][1],cellNodeCoords[0][2]);
		std::cout<<"Coords 2  %f %f %f \n",cellNodeCoords[1][0],cellNodeCoords[1][1],cellNodeCoords[1][2]);
		std::cout<<"Coords 3  %f %f %f \n",cellNodeCoords[2][0],cellNodeCoords[2][1],cellNodeCoords[2][2]);
	      }
	    else
	      {
		std::cout<<"refine false %f %f \n",point[0],point[1],check1,check2,check3);
		std::cout<<"Coords 1  %f %f %f \n",cellNodeCoords[0][0],cellNodeCoords[0][1],cellNodeCoords[0][2]);
		std::cout<<"Coords 2  %f %f %f \n",cellNodeCoords[1][0],cellNodeCoords[1][1],cellNodeCoords[1][2]);
		std::cout<<"Coords 3  %f %f %f \n",cellNodeCoords[2][0],cellNodeCoords[2][1],cellNodeCoords[2][2]);
	      }

	    }
	  flagTet = false;


	  if(refine)return refine;
	  //if(refine)break;
	  */
	}
    }

  subSurfList.clear();
  return refine;

}

//-------------------------------------------------------------------------------------
//  base refinement on distance to the centroid of the set of surfaces in memory
//-------------------------------------------------------------------------------------


double meshRefine::doIrefineCentroidSidedFinishMetric()
{

  bool refine=false;
  bool pRefine=false;
  bool nRefine=false;
  bool iLRefine=false;
  bool checkPSide=false;
  bool checkNSide=false;
  bool checkILSide=false;
  double maxSide = getMaxSide();
  double minSide = getMinSide();
  double aveSide = getAveSide(true);
  std::vector<double> point;
  double metric=0.0;
  //find a subset of surfaces to check by centroid
  std::vector<int> subSurfList;
  int numSubs=10;

  point.resize(3);

  subSurfList.resize(numSubs);

  findClosest(numSubs,subSurfList);

  flagTet=false;

  if(!cellMinNSet)
    {
      std::cout  <<"Must set n-side cell minimum size!"<<std::endl;
      exit(1);
    }

  if(!cellMinPSet)
    {
      std::cout  <<"Must set p-side cell minimum size!"<<std::endl;
      exit(1);
    }

  double localCoords[] = {cxc,cyc,czc};

  //for(int surf=0 ; surf < numSurfs ; surf++)
 for(int surf=0 ; surf < numSubs ; surf++)
    {
      for(int i=0 ; i<numCellNodes ; ++i)
	{

	  localCoords[0] = cellNodeCoords[i][0]; 
	  localCoords[1] = cellNodeCoords[i][1]; 
	  localCoords[2] = cellNodeCoords[i][2]; 

	  checkPSide=false;
	  checkNSide=false;
	  checkILSide=false;
	  surfaceInfo localSurf = surfs[subSurfList[surf]];
	  //Determine which type of surf we're looking at
	  //	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],surfs[surf]);
	  int type = getSurfType(localCoords[0],localCoords[1],localCoords[2],localSurf);
	  double rDist = 0;
	  double rSize=0;
	  if(type == 0)  //p type
	    {
	      rDist = refinementDistanceP;
	      rSize = cellMinSizeP;
	      rDist = localSurf.pThick;
	      rSize = localSurf.pThick/10;
	      checkPSide=true;
	    }
	  if(type == 1)  //n type
	    {
	      rDist = refinementDistanceN;
	      rSize = cellMinSizeN;
	      rDist = localSurf.nThick;
	      rSize = localSurf.nThick/10;
	      checkNSide=true;
	    }
	  if(type == 2)  //inversion layer
	    {
	      rDist = refinementDistanceIL;
	      rSize = cellMinSizeIL;
	      checkILSide=true;
	    }


	  double angle = 0;
	  //This is not elegant, but if it's N check, flip the normal
	  if(checkPSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
				      localSurf.normy,localSurf.normz);  //cosine is even function
	  if(checkNSide)
	    angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
				      localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
				      -localSurf.normy,-localSurf.normz);  //cosine is even function
	  if(checkILSide)
	    {
	      angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					localSurf.xc,localSurf.yc,localSurf.zc,localSurf.normx,
					localSurf.normy,localSurf.normz);  //cosine is even function
	      if(cos(angle) < 0)
		angle = getVectorRadAngle(localCoords[0],localCoords[1],localCoords[2],
					  localSurf.xc,localSurf.yc,localSurf.zc,-localSurf.normx,
					  -localSurf.normy,-localSurf.normz);  //cosine is even function
	    }

	  if(angle<0)
	    std::cout<<" SOMETHING HAS GONE WRONG "<<std::endl;

	  //double checkDistance = rDist*cos(angle) + localSurf.size*sin(angle);

	  double checkDistance = sqrt((localCoords[0]-localSurf.xc)*(localCoords[0]-localSurf.xc) +
				      (localCoords[1]-localSurf.yc)*(localCoords[1]-localSurf.yc) +
				      (localCoords[2]-localSurf.zc)*(localCoords[2]-localSurf.zc))*cos(angle);

	  bool refineOrNot=false;

	  if(checkDistance <= rDist)
	    refineOrNot = true;

	  //Special consideration for oxide adjacent cells
	  if(checkILSide && cos(angle) < 1e-10 )
	    {
	      if(distance(localCoords[0],localCoords[1],localCoords[2],
			  localSurf.xc,localSurf.yc,localSurf.zc) < localSurf.size)
		{
		  // This node sits on the oxide;
		  //std::cout<<"Checking against inversion layer surface,  %f   %f \n",angle,cos(angle));
		  refineOrNot = true;
		  checkDistance = 0;
		}
	    }

	  //If it's close enough, check if it's in the cylinder
	  if(refineOrNot)
	    {
	      double xcomp = localCoords[0] - (localSurf.xc + localSurf.normx*checkDistance);
	      double ycomp = localCoords[1] - (localSurf.yc + localSurf.normy*checkDistance);
	      double zcomp = localCoords[2] - (localSurf.zc + localSurf.normz*checkDistance);

	      double cylinder = sqrt(xcomp*xcomp + ycomp*ycomp + zcomp*zcomp);

	      //Outside the cylinder; shut down refinement
	      if(cylinder > localSurf.size)
		refineOrNot = false;
	    }

	  if(refineOrNot)
	    {

	      refine = true;
	      metric = aveSide/rSize;
	    }
	  if(refine)break;
	}
      if(refine)break;
    }

  //In cases where the tet is large compared to the refinement distance, the surface may lie
  //inside the tet without any node within the refinement distance.  Need to check it.

  //If the decision to refine has already been made, skip this
  if(!refine)
    {
      //for(int surf=0 ; surf < numSurfs ; surf++)
      for(int surf=0 ; surf < numSubs ; surf++)
	{
	  //Point Check

	  surfaceInfo localSurf = surfs[subSurfList[surf]];

	  point[0] = localSurf.xc;;
	  point[1] = localSurf.yc;
	  point[2] = localSurf.zc;
      
	  bool check1 = signedDistanceSameSideCheck(point,cellNodeCoords[0],cellNodeCoords[1],
						    cellNodeCoords[2],cellNodeCoords[3]);
	  bool check2 = signedDistanceSameSideCheck(point,cellNodeCoords[1],cellNodeCoords[0],
						    cellNodeCoords[2],cellNodeCoords[3]);
	  bool check3 = signedDistanceSameSideCheck(point,cellNodeCoords[2],cellNodeCoords[1],
						    cellNodeCoords[0],cellNodeCoords[3]);
	  bool check4 = signedDistanceSameSideCheck(point,cellNodeCoords[3],cellNodeCoords[1],
						    cellNodeCoords[2],cellNodeCoords[0]);
      
	  double sizeCheck = 0;
	  if(cellMinSizeN < cellMinSizeP)
	    sizeCheck = cellMinSizeN;
	  else
	    sizeCheck = cellMinSizeP;

	  if(check1 && check2 && check3 && check4)
	    {
	      refine = true;
	      metric = aveSide/sizeCheck;
	    }
	  if(refine) break;
	}
    }

  subSurfList.clear();

  if(!refine)metric = 0.0;

  return metric;

}


//-------------------------------------------------------------------------------------
//  calculate the distance between two points.
//-------------------------------------------------------------------------------------

double meshRefine::distance(double x1,double y1,double z1,double x2,double y2,double z2)
{
  double dist = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
  return dist;
}

//-------------------------------------------------------------------------------------
//  return the length metric to be used in refinement decision making
//-------------------------------------------------------------------------------------

double meshRefine::getLengthMetric()
{

  double metric;

  switch(sizeMeasure)
    {
    case 0:
      metric = getAveSide(true);
      break;

    case 1:
      metric = getMinSide();
      break;

    case 2:
      metric = getMaxSide();
      break;

    case 3:
      metric = getAveSide(false);
      break;

    }

  return metric;

}

//-------------------------------------------------------------------------------------
//  calculate the distance between two points.
//-------------------------------------------------------------------------------------

double meshRefine::getMaxSide()
{

  double maxDist = -1.0e20;
  double dist;

  int numSides;
  int iNum;
  int jNum;

  if(is3D)
    {
      numSides = 6;
      iNum=3;
      jNum = 4;
    }
  else
    {
      numSides=3;
      iNum=2;
      jNum=3;
    }

  for(int i=0 ; i<iNum ; ++i)
    for(int j=i+1 ; j<jNum ; ++j)
      {
	dist = distance(cellNodeCoords[i][0],cellNodeCoords[i][1],cellNodeCoords[i][2],
			cellNodeCoords[j][0],cellNodeCoords[j][1],cellNodeCoords[j][2]);

	if(dist > maxDist)
	  maxDist = dist;
      }

  return maxDist;

}

//-------------------------------------------------------------------------------------
//  calculate the distance between two points.
//-------------------------------------------------------------------------------------

double meshRefine::getAveSide(bool getAll)
{

  double aveDist = 0.0;
  double dist=0;
  int numSides;
  int iNum;
  int jNum;
  int numAve;
  std::vector<double> sideLength;

  if(is3D)
    {
      numSides = 6;
      iNum=3;
      jNum = 4;
    }
  else
    {
      numSides=3;
      iNum=2;
      jNum=3;
    }

  if(getAll)
    numAve = numSides;
  else
    numAve = 2;

  for(int i=0 ; i<iNum ; ++i)
    for(int j=i+1 ; j<jNum ; ++j)
      {
	dist = distance(cellNodeCoords[i][0],cellNodeCoords[i][1],cellNodeCoords[i][2],
			cellNodeCoords[j][0],cellNodeCoords[j][1],cellNodeCoords[j][2]);

	sideLength.push_back(dist);

      }

  std::sort(sideLength.begin(), sideLength.end());

  for(size_t i=0 ; i<numSides ; ++i)
    aveDist += sideLength[numSides-i-1];

  return aveDist/(double)numAve;

}



//-------------------------------------------------------------------------------------
//  calculate the distance between two points.
//-------------------------------------------------------------------------------------

double meshRefine::getMinSide()
{

  double minDist = 1.0e20;
  double dist=0;
  int numSides;
  int iNum;
  int jNum;

  if(is3D)
    {
      numSides = 6;
      iNum=3;
      jNum = 4;
    }
  else
    {
      numSides=3;
      iNum=2;
      jNum=3;
    }


  for(int i=0 ; i<iNum ; ++i)
    for(int j=i+1 ; j<jNum ; ++j)
      {
	dist = distance(cellNodeCoords[i][0],cellNodeCoords[i][1],cellNodeCoords[i][2],
			cellNodeCoords[j][0],cellNodeCoords[j][1],cellNodeCoords[j][2]);

	if(dist < minDist)
	  minDist = dist;
      }

  return minDist;

}

//-------------------------------------------------------------------------------------
// Check which side a node is on.
//-------------------------------------------------------------------------------------

int meshRefine::checkSide(std::vector<double> norm, std::vector<double> point, std::vector<double> coord)
{

  double vector[3];

  vector[0] = coord[0] - point[0];
  vector[1] = coord[1] - point[1];
  vector[2] = coord[2] - point[2];

  double dotProduct = norm[0]*vector[0] + norm[1]*vector[1] + norm[2]*vector[2];

  if(dotProduct > 0)
    return 1;
  else
    return -1;

}


//-------------------------------------------------------------------------------------
//  Check to see if a point and a node are on the same side of a tet face
//-------------------------------------------------------------------------------------

bool meshRefine::signedDistanceSameSideCheck(std::vector<double> pointm1, std::vector<double> point0, std::vector<double> point1, std::vector<double> point2, std::vector<double> point3)
{
  
  bool sameSide=false;
  double t1[3],t2[3],n1[3],cent[3];

  //Compute the point in the center of triangle

  cent[0] = (point1[0] + point2[0] + point3[0])/3.0;
  cent[1] = (point1[1] + point2[1] + point3[1])/3.0;
  cent[2] = (point1[2] + point2[2] + point3[2])/3.0;

  //Compute two tangent vectors

  t1[0] = point1[0] - point2[0];
  t1[1] = point1[1] - point2[1];
  t1[2] = point1[2] - point2[2];

  t2[0] = point1[0] - point3[0];
  t2[1] = point1[1] - point3[1];
  t2[2] = point1[2] - point3[2];

  //now compute the cross product to get the normal vector

  n1[0] =   t1[1]*t2[2] - t2[1]*t1[2];
  n1[1] = -(t1[0]*t2[2] - t2[0]*t1[2]);
  n1[2] =   t1[0]*t2[1] - t2[0]*t1[1];

  //Do a quick sanity check

  double normalLength2 = n1[0]*n1[0] + n1[1]*n1[1] +n1[2]*n1[2]; 

  if(normalLength2 < 1.0e-16)
    std::cout<<" Invalid description of plane in signed distance "<<normalLength2<<std::endl;

  //If I really cared about the signed distance, I'd normalize n1.
  //But I really ony care about signs in this case

  //Check pointm1

  double vector[3];

  vector[0] = pointm1[0] - cent[0];
  vector[1] = pointm1[1] - cent[1];
  vector[2] = pointm1[2] - cent[2];

  double dist1 = vector[0]*n1[0] + vector[1]*n1[1] + vector[2]*n1[2];

  //check point0

  vector[0] = point0[0] - cent[0];
  vector[1] = point0[1] - cent[1];
  vector[2] = point0[2] - cent[2];

  double dist2 = vector[0]*n1[0] + vector[1]*n1[1] + vector[2]*n1[2];

  if(dist1*dist2 >= 0)
    sameSide = true;
  else
    sameSide = false;

  return sameSide;

}

int meshRefine::getSurfType(double x, double y, double z, surfaceInfo surf)
{

  if(surf.surfType == 1) //This is inversion layer surface
    return 2;

  double dotProduct = surf.normx*(x-surf.xc) + surf.normy*(y-surf.yc) + surf.normz*(z-surf.zc);

  if(dotProduct>=0) 
    return 0;  //p type
  else
    return 1;  //n type
}

//----------------------------------------------------------
// Find closest surfs to tet by centroid
//----------------------------------------------------------

void meshRefine::findClosest(int numSubs,std::vector<int> &subSurfList)
{

  double xp=0,yp=0,zp=0,dist;
  double surfdist[numSubs];
  int maxi=0;

  int numCellNodes;
  if(is3D)
    numCellNodes = 4;
  else
    numCellNodes = 3;

  for(int i=0 ; i<numSubs ; ++i)
    surfdist[i] = 1e10;

  for(int i=0 ; i<numCellNodes ; ++i)
    {
      xp += cellNodeCoords[i][0];
      yp += cellNodeCoords[i][1];
      zp += cellNodeCoords[i][2];
    }

  xp /= (double)numCellNodes;
  yp /= (double)numCellNodes;
  zp /= (double)numCellNodes;

  //Temporary
  cxc = xp;
  cyc = yp;
  czc = zp;


  for(int surf_count=0 ; surf_count < numSurfs ; ++surf_count)
    {
      dist = distance(surfs[surf_count].xc,surfs[surf_count].yc,surfs[surf_count].zc,xp,yp,zp);

     for (int i=0 ; i<numSubs ; ++i)
	{
	  if(dist <= surfdist[i])
	    {
	      ++maxi;
	      for(int j=numSubs-1 ; j>i ; --j)
		{
		  surfdist[j] = surfdist[j-1];
		  subSurfList[j] = subSurfList[j-1];
		}
	      surfdist[i] = dist;
	      subSurfList[i] = surf_count;
	      break;
	    }
	}
    }


  return;

}


//-------------------------------------------------------------------
//return the angle between two vectors
//-------------------------------------------------------------------
double meshRefine::getVectorRadAngle(double x1, double y1, double z1, 
			 double x2, double y2, double z2,
			 double normx, double normy, double normz)
{

  double vecSize1 = distance(x1, y1,z1,x2,y2,z2);

  double vec1x = (x1-x2)/vecSize1;
  double vec1y = (y1-y2)/vecSize1;
  double vec1z = (z1-z2)/vecSize1;

  double vecSize2 = sqrt(normx*normx + normy*normy +normz*normz);

  double vec2x = normx/vecSize2;
  double vec2y = normy/vecSize2;
  double vec2z = normz/vecSize2;

  double dotProduct = vec1x*vec2x + vec1y*vec2y + vec1z*vec2z;

  double radAngle = acos(dotProduct);

  return radAngle;

}


//-------------------------------------------------------------------
//call methods to create srufaces from analytically defined doping
//-------------------------------------------------------------------
void meshRefine::createSurfaces()
{

  if(meshDimension == 2)
    createSurfaces2D();
  else
    createSurfaces3D();


}


//-------------------------------------------------------------------
//return a set of generated surfaces in 2D (line segments
//-------------------------------------------------------------------
void meshRefine::createSurfaces2D()
{

  double xlo =  BBxmin;
  double xhi =  BBxmax;
  double ylo =  BBymin;
  double yhi =  BBymax;

  cellCreator cC;
  if(maxSurfaceRecursions > 0)
    cC.setMaxLevels(maxSurfaceRecursions);
  if(guaranteedSurfaceRecursions > 0)
    {
      if(guaranteedSurfaceRecursions > maxSurfaceRecursions)
	cC.setGuaranteedRecursions(maxSurfaceRecursions);
      else
	cC.setGuaranteedRecursions(guaranteedSurfaceRecursions);
    }

  std::vector<cells> cellData;
  cellData = cC.createCells(xlo,xhi,ylo,yhi,function);

  isolineFinder isoF;

  std::vector<surfaceInfo> localSurfs;
  std::vector<double> x,y,f;
  x.resize(4);
  y.resize(4);
  f.resize(4);

  int vertexCount = 1;

  bool journalCells=false;
  //bool journalSurfaces=true;
  bool journalSurfaces = writeJunctions;

  std::ofstream surfaceJournalFile;
  if(journalSurfaces)
    surfaceJournalFile.open("junctionSurfaces.jou");

  for (size_t icell=0 ; icell<cellData.size() ; ++icell)
    {

      x[0] = cellData[icell].vertices[0][0];
      x[1] = cellData[icell].vertices[1][0];
      x[2] = cellData[icell].vertices[2][0];
      x[3] = cellData[icell].vertices[3][0];

      y[0] = cellData[icell].vertices[0][1];
      y[1] = cellData[icell].vertices[1][1];
      y[2] = cellData[icell].vertices[2][1];
      y[3] = cellData[icell].vertices[3][1];

      f[0] = cellData[icell].values[0];
      f[1] = cellData[icell].values[1];
      f[2] = cellData[icell].values[2];
      f[3] = cellData[icell].values[3];

      /*
      if(y[0] > -2)
	std::cout<<" 0  "<<x[0]<<"    "<<y[0]<<"    "<<f[0]<<std::endl;
      if(y[1] > -2)
	std::cout<<" 1  "<<x[1]<<"    "<<y[1]<<"    "<<f[1]<<std::endl;
      if(y[2] > -2)
	std::cout<<" 2  "<<x[2]<<"    "<<y[2]<<"    "<<f[2]<<std::endl;
      if(y[3] > -2)
	std::cout<<" 3  "<<x[3]<<"    "<<y[3]<<"    "<<f[3]<<std::endl;
      */

      if(journalCells)
	{
	  std::cout<<"create vertex "<<x[0]<<"    "<<y[0]<<std::endl;
	  std::cout<<"create vertex "<<x[1]<<"    "<<y[1]<<std::endl;
	  std::cout<<"create vertex "<<x[2]<<"    "<<y[2]<<std::endl;
	  std::cout<<"create vertex "<<x[3]<<"    "<<y[3]<<std::endl;
	  //std::cout<<"create surface vertex "<<vertexCount<<"    "<<vertexCount+1<<"    "<<vertexCount+2<<"    "<<vertexCount+3<<std::endl;
	  std::cout<<"create curve vertex "<<vertexCount<<"    "<<vertexCount+1<<std::endl;
	  std::cout<<"create curve vertex "<<vertexCount+1<<"    "<<vertexCount+2<<std::endl;
	  std::cout<<"create curve vertex "<<vertexCount+2<<"    "<<vertexCount+3<<std::endl;
	  std::cout<<"create curve vertex "<<vertexCount+3<<"    "<<vertexCount<<std::endl;
	  vertexCount += 8;
	}

      isoF.setCell(x,y,f);

      localSurfs = isoF.createSurfaces(0.0);

      for(size_t ls=0 ; ls<localSurfs.size() ; ++ls)
	{
	  localSurfs[ls].surfType = 0;  //junction surface
	  surfs.push_back(localSurfs[ls]);
	}

    }

  if(journalSurfaces)
    {
      for(size_t sf=0 ; sf<surfs.size() ; ++sf)
	{
	  surfaceJournalFile<<"create vertex "<<surfs[sf].X[0]<<"    "<<surfs[sf].Y[0]<<std::endl;
	  surfaceJournalFile<<"create vertex "<<surfs[sf].X[1]<<"    "<<surfs[sf].Y[1]<<std::endl;
	  surfaceJournalFile<<"create curve vertex "<<vertexCount<<"    "<<vertexCount+1<<std::endl;
	  vertexCount += 2;
	}
    }

  numSurfs = surfs.size();
  std::cout<<" I created "<<surfs.size()<<" surfaces to refine to."<<std::endl;


  if(journalSurfaces)
    {
      surfaceJournalFile<<"export acis \"jS.sat\" overwrite"<<std::endl;
      surfaceJournalFile<<"exit"<<std::endl;
      surfaceJournalFile.close();
    }

  depletionWidth dW;
  double stepOff = 1e-2;
  double nDoping,pDoping;
  std::vector<double> xOff(3,0.0),norm(3,0.0);

  for (size_t surf=0 ; surf<surfs.size() ; ++surf)
    {
      xOff[0] = surfs[surf].xc + stepOff*surfs[surf].normx;
      xOff[1] = surfs[surf].yc + stepOff*surfs[surf].normy;
      double dope = function.evaluateFunction(xOff);
      if (dope >= 0)
	nDoping = dope;
      if (dope < 0)
	pDoping = -dope;

      xOff[0] = surfs[surf].xc - stepOff*surfs[surf].normx;
      xOff[1] = surfs[surf].yc - stepOff*surfs[surf].normy;
      dope = function.evaluateFunction(xOff);
      if (dope >= 0)
	nDoping = dope;
      if (dope < 0)
	pDoping = -dope;

      double thick,nThick,pThick;

      dW.calculateDepletionWidths(pDoping,nDoping,thick,pThick,nThick);

      surfs[surf].nThick = 1.0e4*nThick;  //convert from cm to um
      surfs[surf].pThick = 1.0e4*pThick;  //convert from cm to um
    }
  

}


//-------------------------------------------------------------------
//create surfaces for specified line segment location
//-------------------------------------------------------------------

void meshRefine::addRefineToLine(double xmin, double ymin, double xmax, double ymax, double ilThick)
{

  double sLength;
  int addedSegments = 1;

  double totalLength = sqrt((xmax - xmin)*(xmax - xmin) + (ymax - ymin)*(ymax - ymin));

  double xIncrements = (xmax-xmin)/(double)addedSegments;
  double yIncrements = (ymax-ymin)/(double)addedSegments;

  double x=xmin,y=ymin;

  for (int ls=0 ; ls<addedSegments ; ++ls)
    {

      surfaceInfo localSurf;
      localSurf.X.push_back(x);
      localSurf.Y.push_back(y);
      localSurf.Z.push_back(0);

      x += xIncrements;
      y += yIncrements;

      localSurf.X.push_back(x);
      localSurf.Y.push_back(y);
      localSurf.Z.push_back(0);

      localSurf.xc = 0.5*(localSurf.X[0] + localSurf.X[1]);
      localSurf.yc = 0.5*(localSurf.Y[0] + localSurf.Y[1]);
      localSurf.zc = 0;

      localSurf.normx =  0;
      localSurf.normy = -1;
      localSurf.normz =  0;

      localSurf.pThick = 0.1;
      localSurf.nThick = 0.1;

      localSurf.surfType = 2;

      localSurf.refinementDistance = ilThick;

      localSurf.size = sqrt((localSurf.X[0] - localSurf.X[1])*(localSurf.X[0] - localSurf.X[1]) 
			    + (localSurf.Y[0] - localSurf.Y[1])*(localSurf.Y[0] - localSurf.Y[1]));

      surfs.push_back(localSurf);

    }

  numSurfs += addedSegments;
  std::cout<<" I created "<<surfs.size()<<" SPECIFIED surfaces to refine to."<<std::endl;
}


//-------------------------------------------------------------------
//create surfaces for refinement
//-------------------------------------------------------------------

void meshRefine::addRefineToSurface(double x0, double y0, double z0, double x1, double y1, double z1, 
				    double x2, double y2, double z2, double x3, double y3, double z3,
				    double ilThick)
{

  surfaceInfo localSurf0, localSurf1;
  localSurf0.X.push_back(x0);
  localSurf0.Y.push_back(y0);
  localSurf0.Z.push_back(z0);

  localSurf0.X.push_back(x1);
  localSurf0.Y.push_back(y1);
  localSurf0.Z.push_back(z1);

  localSurf0.X.push_back(x2);
  localSurf0.Y.push_back(y2);
  localSurf0.Z.push_back(z2);

  localSurf1.X.push_back(x0);
  localSurf1.Y.push_back(y0);
  localSurf1.Z.push_back(z0);

  localSurf1.X.push_back(x2);
  localSurf1.Y.push_back(y2);
  localSurf1.Z.push_back(z2);

  localSurf1.X.push_back(x3);
  localSurf1.Y.push_back(y3);
  localSurf1.Z.push_back(z3);

  localSurf0.normx =  0;
  localSurf0.normy = -1;
  localSurf0.normz =  0;

  localSurf1.normx =  0;
  localSurf1.normy = -1;
  localSurf1.normz =  0;

  localSurf0.pThick = ilThick;
  localSurf0.nThick = ilThick;

  localSurf1.pThick = ilThick;
  localSurf1.nThick = ilThick;

  localSurf0.surfType = 2;
  localSurf1.surfType = 2;

  localSurf0.refinementDistance = ilThick;
  localSurf1.refinementDistance = ilThick;

  surfs.push_back(localSurf0);
  surfs.push_back(localSurf1);

  std::cout<<" I created "<<2<<" SPECIFIED surfaces to refine to."<<std::endl;
}



//-------------------------------------------------------------------
//return a set of generated surfaces in 3D
//-------------------------------------------------------------------
void meshRefine::createSurfaces3D()
{

  //  scalarFunction function;

  surfaceFinder surface;

  int node_number[8]={0,3,2,1,4,7,6,5};

  double vertices[8][3],rvertices[8][3];

  double values[8];

  double threshold=0.0; 
  double epsilon=1e-8;
  int *total_verts;
  edgeStruct **edge_list=NULL;
  int *total_tris;
  int **triangle_list=NULL;

  edge_list =  new edgeStruct*;
  triangle_list = new int*;
  *edge_list = NULL;
  *triangle_list = NULL;

  total_verts = new int;
  total_tris = new int;

  double xlo = BBxmin;
  double xhi = BBxmax;
  double ylo = BBymin;
  double yhi = BBymax;
  double zlo = BBzmin;
  double zhi = BBzmax;

  cellCreator cC;
  if(maxSurfaceRecursions > 0)
    cC.setMaxLevels(maxSurfaceRecursions);
  if(guaranteedSurfaceRecursions > 0)
    {
      if(guaranteedSurfaceRecursions > maxSurfaceRecursions)
	cC.setGuaranteedRecursions(maxSurfaceRecursions);
      else
	cC.setGuaranteedRecursions(guaranteedSurfaceRecursions);
    }

  std::vector<cells> cellData;
  cellData = cC.createCells(xlo,xhi,ylo,yhi,zlo,zhi,function);

  std::cout<<" I created "<<cellData.size()<<" 3D cells."<<std::endl;

  int vertex_counter=1;

  bool journalSurfaces = writeJunctions;

  std::ofstream surfaceJournalFile;
  if(journalSurfaces)
    surfaceJournalFile.open("junctionSurfaces.jou");


  for (size_t icell=0 ; icell<cellData.size() ; ++icell)
    {

      surface.draw_cell(HEX,cellData[icell].vertices,cellData[icell].values,threshold,epsilon,total_verts,
			edge_list,total_tris,triangle_list);

      int *tptr;
      tptr = *triangle_list;

      for(int i=0 ; i< *total_tris ; ++i)
	{

	  surfaceInfo localSurf;
	  localSurf.numNodes = 3;
	  localSurf.xc = 0.0;
	  localSurf.yc = 0.0;
	  localSurf.zc = 0.0;
	  localSurf.X.resize(3);
	  localSurf.Y.resize(3);
	  localSurf.Z.resize(3);

	  for(int inode=0 ; inode<3 ; ++inode)
	    {
	      localSurf.xc += (*edge_list)[(*tptr)].x[0];
	      localSurf.yc += (*edge_list)[(*tptr)].x[1];
	      localSurf.zc += (*edge_list)[(*tptr)].x[2];
	      localSurf.X[inode] = (*edge_list)[(*tptr)].x[0];
	      localSurf.Y[inode] = (*edge_list)[(*tptr)].x[1];
	      localSurf.Z[inode] = (*edge_list)[(*tptr)].x[2];
	      tptr++;

	      if(journalSurfaces)
		surfaceJournalFile<<"create vertex "<<localSurf.X[inode]<<"    "<<localSurf.Y[inode]<<"    "<<localSurf.Z[inode]<<std::endl;
	    }
	  if(journalSurfaces)
	    surfaceJournalFile<<"create surface vertex "<<vertex_counter<<"    "<<vertex_counter+1<<"    "<<vertex_counter+2<<std::endl;

	  vertex_counter += 3;

	  localSurf.xc /= 3.0;
	  localSurf.yc /= 3.0;
	  localSurf.zc /= 3.0;
	  localSurf.nThick = 0.1;
	  localSurf.pThick = 0.1;
	  localSurf.surfType = 0;

	  //Compute and store the normal to the surface

	  double vecx[]={0,0},vecy[]={0,0},vecz[]={0,0};

	  vecx[0] = localSurf.X[1] - localSurf.X[0];
	  vecy[0] = localSurf.Y[1] - localSurf.Y[0];
	  vecz[0] = localSurf.Z[1] - localSurf.Z[0];
      
	  vecx[1] = localSurf.X[2] - localSurf.X[0];
	  vecy[1] = localSurf.Y[2] - localSurf.Y[0];
	  vecz[1] = localSurf.Z[2] - localSurf.Z[0];
  
      
	  //Now create normal vector from cross product

	  localSurf.normx =   vecy[0]*vecz[1] - vecz[0]*vecy[1];
	  localSurf.normy = -(vecx[0]*vecz[1] - vecz[0]*vecx[1]);
	  localSurf.normz =   vecx[0]*vecy[1] - vecy[0]*vecx[1];

	  double norm = sqrt(localSurf.normx*localSurf.normx + 
			     localSurf.normy*localSurf.normy + 
			     localSurf.normz*localSurf.normz);
	  localSurf.normx /= norm;
	  localSurf.normy /= norm;
	  localSurf.normz /= norm;
	  
	  //Check orientation of the norm

	  for(int inode=0 ; inode<numCellNodes ; ++inode)
	    {
	      bool checked = false;
	      bool switchsign = false;


	      if(cellData[icell].values[inode] < 0)
		{
		  double vec1x = cellData[icell].vertices[inode][0] - localSurf.xc;
		  double vec1y = cellData[icell].vertices[inode][1] - localSurf.yc;
		  double vec1z = cellData[icell].vertices[inode][2] - localSurf.zc;

		  double check = vec1x*localSurf.normx + vec1y*localSurf.normy + vec1z*localSurf.normz;

		  checked = true;
		  if(check < 0.0)switchsign = true;
		}

	      if(switchsign)
		{
		  localSurf.normx *= -1;
		  localSurf.normy *= -1;
		  localSurf.normz *= -1;
		}
	      if(checked) break;
	    }


	  localSurf.surfType = 0; //junction surface

	  surfs.push_back(localSurf);


	}

    }

  if(journalSurfaces)
    {
      //surfaceJournalFile<<"export acis \"junctionSurfaces.sat\" overwrite"<<std::endl;
      surfaceJournalFile<<"unite body all"<<std::endl;
      surfaceJournalFile<<"export acis \"jS.sat\" overwrite"<<std::endl;
      surfaceJournalFile<<"exit"<<std::endl;
      surfaceJournalFile.close();
    }

  depletionWidth dW;
  double stepOff = 1e-2;
  double nDoping,pDoping;
  std::vector<double> xOff(3,0.0),norm(3,0.0);

  for (size_t surf=0 ; surf<surfs.size() ; ++surf)
    {
      xOff[0] = surfs[surf].xc + stepOff*surfs[surf].normx;
      xOff[1] = surfs[surf].yc + stepOff*surfs[surf].normy;
      xOff[2] = surfs[surf].zc + stepOff*surfs[surf].normz;
      double dope = function.evaluateFunction(xOff);
      if (dope >= 0)
	nDoping = dope;
      if (dope < 0)
	pDoping = -dope;

      xOff[0] = surfs[surf].xc - stepOff*surfs[surf].normx;
      xOff[1] = surfs[surf].yc - stepOff*surfs[surf].normy;
      xOff[2] = surfs[surf].zc + stepOff*surfs[surf].normz;
      dope = function.evaluateFunction(xOff);
      if (dope >= 0)
	nDoping = dope;
      if (dope < 0)
	pDoping = -dope;

      double thick,nThick,pThick;

      dW.calculateDepletionWidths(pDoping,nDoping,thick,pThick,nThick);

      surfs[surf].nThick = 1.0e4*nThick;  //convert from cm to um
      surfs[surf].pThick = 1.0e4*pThick;  //convert from cm to um
    }
  
  numSurfs = surfs.size();
  std::cout<<" I created "<<surfs.size()<<" surfaces to refine to."<<std::endl;

}


//-------------------------------------------------------------------
//create a journal file of the surfaces
//-------------------------------------------------------------------
void meshRefine::journalSurfaces()
{

  std::ofstream ofile("surfs.jou");
  int vertexCounter=0;

  for(int isurf=0 ; isurf<surfs.size() ; ++isurf)
    {
      
      for(int inode = 0; inode<3 ; ++inode)
	ofile<<"create vertex  "<<surfs[isurf].X[inode]<<"    "<<surfs[isurf].Y[inode]<<"    "<<surfs[isurf].Z[inode]<<std::endl;
      ofile<<"create surface vertex "<<vertexCounter+1<<"     "<<vertexCounter+2<<"    "<<vertexCounter+3<<std::endl;

      vertexCounter += 3;

    }

  ofile.close();

}



//-----------------------------------------------------------------------
// create the doping function
//-----------------------------------------------------------------------

void meshRefine::setDopingFunction(std::string name, std::string functionType, std::string dopingType)
{
  function.setDopingFunction(name, functionType, dopingType);
}

//-----------------------------------------------------------------------
// set the minimum and maximum geometric bounds
//-----------------------------------------------------------------------
void meshRefine::setXBounds(std::string name, std::string axis, double min, double max)
{
  function.setXBounds(name, axis, min, max);
}

//-----------------------------------------------------------------------
// set the minimum and maximum doping concentrations for this function
//-----------------------------------------------------------------------
void meshRefine::setDopingBounds(std::string name, double min, double max)
{
  function.setDopingBounds(name, min,  max);
}



//-----------------------------------------------------------------------
// set the coordinate direction (+/-) of the function
//-----------------------------------------------------------------------

void meshRefine::setDirection(std::string name, std::string direction)
{
  function.setDirection(name,direction);
}



//-----------------------------------------------------------------------
// set function location--a translation from the origin--3D
//-----------------------------------------------------------------------
void meshRefine::setDopingLocation(std::string name, double x, double y, double z)
{
  function.setLocation(name, x, y, z);
}


//-----------------------------------------------------------------------
// set function location--a translation from the origin--2D
//-----------------------------------------------------------------------
void meshRefine::setDopingLocation(std::string name, double x, double y)
{
  function.setLocation(name, x, y);
}


//-----------------------------------------------------------------------
// set function width
//-----------------------------------------------------------------------
void meshRefine::setDopingWidth(std::string name, double x, double y, double z)
{
  function.setWidth(name, x, y, z);
}


//-----------------------------------------------------------------------
// set function width
//-----------------------------------------------------------------------
void meshRefine::setDopingWidth(std::string name, double x, double y)
{
  function.setWidth(name, x, y);
}


//-----------------------------------------------------------------------
// set function width
//-----------------------------------------------------------------------
void meshRefine::listFunctions()
{
  function.listFunctions();
}

//-----------------------------------------------------------------------
// print doping function XML equivalents
//-----------------------------------------------------------------------
void meshRefine::printXMLDopingFunctions()
{
  function.printXMLFunctions();
}

