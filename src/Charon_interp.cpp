
///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <algorithm>

// Charon
#include "Charon_interp.hpp"
#include "Charon_Vector.hpp"

//-----------------------------------------------------------------------------
//Methods for the cluster interpolator class
//-----------------------------------------------------------------------------

//------------------------------------------
//Set interpolant method
//------------------------------------------
bool charon::clusterInterpolator::setMethod(std::string methodName, float shepardPneg)
{

  bool success = interpolantMethodFactory(methodName, shepardPneg);

  return success;

}

//------------------------------------------
//Set influence radius
//------------------------------------------
bool charon::clusterInterpolator::setInfluenceRadius(double iR)
{

  influenceRadius = iR;

  if(iR > 0.0)
    influenceRadiusSet = true;

  return true;

}

//------------------------------------------
//set data file type
//------------------------------------------
bool charon::clusterInterpolator::setFileType(std::string fileTypeName)
{

  bool success = readFileTypeFactory(fileTypeName);

  return success;

}

//------------------------------------------
//Factory to create the interpolant object
//------------------------------------------
bool charon::clusterInterpolator::interpolantMethodFactory(std::string methodName,
                                                           const float& shepardPneg)
{

  bool success=false;

  if(methodName == "SHEPARDS_1D")
  {
    num_dim = 1;
    ipMethod = Teuchos::rcp(new shepardsMethod(num_dim,shepardPneg));
    success = true;
  }
  if(methodName == "SHEPARDS_2D")
  {
    num_dim = 2;
    ipMethod = Teuchos::rcp(new shepardsMethod(num_dim,shepardPneg));
    success = true;
  }
  if(methodName == "SHEPARDS_3D")
  {
    num_dim = 3;
    ipMethod = Teuchos::rcp(new shepardsMethod(num_dim,shepardPneg));
    success = true;
  }
  if(methodName == "ONED_LINEAR_INTERPOLATION")
  {
    num_dim = 1;
    ipMethod = Teuchos::rcp(new oneDLinearInterpolationMethod());
    success = true;
  }

  return success;

}

//------------------------------------------
//Factory to create the file reader object
//------------------------------------------
bool charon::clusterInterpolator::readFileTypeFactory(std::string /* methodName */)
{

  bool success=false;

  rfType = Teuchos::rcp(new clusterFiles(num_dim));
  success = true;

  return success;

}


//------------------------------------------
//interpolate at a point from the clusters
//------------------------------------------
void charon::clusterInterpolator::interpolateToPoint(double xn, double yn,
                                                     double zn, double tn,
                                                     double& pn)
{

  //If the influence radius has not been set, need to calculate it.
  if(!influenceRadiusSet)
    {
      influenceRadius = calculateInfluenceRadius();
    }


  ipMethod->interpolateToPoint(dPS, xn, yn, zn, tn, pn, influenceRadius);

  return;

}



bool charon::clusterInterpolator::readFiles(std::vector<std::string> fileNames)
{

  bool success=false;

  dPS.clear();

  success = rfType->readFiles(fileNames, dPS);

  return success;

}


void charon::clusterInterpolator::InitializeClusterBCVectors(charon::Vector<std::string> clusterNames_)
{

  clusterNames = clusterNames_;
  clusterElectron.resize(clusterNames_.size());
  clusterHole.resize(clusterNames_.size());
  clusterAcceptor.resize(clusterNames_.size());
  clusterDonor.resize(clusterNames_.size());
  clusterFound.resize(clusterNames_.size());

}

//------------------------------------------
//calculate the radius of influence
//------------------------------------------

double charon::clusterInterpolator::calculateInfluenceRadius()
{

  double RofI=-1.0;

  //The calculated influence radius is defined to be twice the greatest minimum distance
  //between any pair of cluster locations.


  //When coupled to xyce models of a cluster, there will be no cluster information early on
  //deal with it.
  if(dPS.size() == 0)
    return RofI;

  RofI = 0.0;

  for(size_t iPS=0 ; iPS<dPS.size() ; ++iPS)
    {

      double minDistForI=1.0e10;

      for(size_t jPS=0 ; jPS<dPS.size() ; ++jPS)
	{
	  double x1 = dPS[iPS].xp;
	  double y1 = dPS[iPS].yp;
	  double z1 = dPS[iPS].zp;
	  double x2 = dPS[jPS].xp;
	  double y2 = dPS[jPS].yp;
	  double z2 = dPS[jPS].zp;

	  //Without loss of generality, can set y, z = 0 if in lower dimensions
	  switch(num_dim)
	    {
	    case 1:
	      y1 = z1 = y2 = z2 = 0.0;
	      break;
	    case 2 :
	      z1 = z2 = 0.0;
	      break;
	    default: //setting distance to zero means RofI is infinite
	      x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
	      break;
	    }

	  double localDist = distance(x1,y1,z1,x2,y2,z2);

	  if(iPS != jPS)
	    {
	      if(localDist < minDistForI)
		minDistForI = localDist;
	    }
	}
      if(minDistForI >  RofI)
	RofI = minDistForI;
      //RofI += minDistForI;
    }


  /*
  if(dPS.size() == 0)
    RofI = -1.0;
  else
    RofI = 2.0*RofI/(double)dPS.size();
  */

  RofI *= 2.0;

  influenceRadiusSet = true;

  return RofI;

}


//------------------------------------------
//calculate the distance between two points
//------------------------------------------

double charon::clusterInterpolator::distance(double x1,double y1,double z1,
					     double x2,double y2,double z2)
{

  double distance = sqrt((x2-x1)*(x2-x1) + 
			 (y2-y1)*(y2-y1) + 
			 (z2-z1)*(z2-z1));

  return distance;

}


//-----------------------------------------------------------------------------
//Methods for the Shepards Method interpolant class
//-----------------------------------------------------------------------------

//------------------------------------------
//interpolate at a point from the clusters
//by Shepard's method
//------------------------------------------
void charon::shepardsMethod::interpolateToPoint(std::vector<dataPointSet> dPS_,
                                                 double xn, double yn,
                                                 double zn, double tn,
						double& pn, double RofI)
{

  //Given a point, loop over the data points and interpolate to the given point

  double dist_to_point,total_distance,x_dist,y_dist,z_dist;

  total_distance = 0.0;
  dist_to_point = 0.0;

  double pfunc, w;

  pfunc = 0.0;
  w = 0.0;

  int timeIndexes[2];

  std::vector<double> timeStepVector;

  //Depart if outside clusters
  if(dPS_.size() > 0)
    {
      if(xn > dPS_[dPS_.size()-1].xp  ||
	 xn < dPS_[0].xp)
	{
	  pn = 0.0;
	  return;
	}
    }
  else
    {
      pn = 0.0;
      return;
    }


  for(size_t ipoint=0 ; ipoint<dPS_.size() ; ++ipoint)
  {

    x_dist = dPS_[ipoint].xp - xn;
    y_dist = dPS_[ipoint].yp - yn;
    z_dist = dPS_[ipoint].zp - zn;

    switch(num_dim)
    {
      case 1:
        dist_to_point = sqrt(x_dist*x_dist); //1D
        break;
      case 2:
        dist_to_point = sqrt(x_dist*x_dist + y_dist*y_dist); //2D
        break;
      case 3:
        dist_to_point = sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist); //3D
        break;
      default:
        std::cout << "Error with number of dimensions in Shepards method" << std::endl;
    }

    //Find the time indexes in between which to interpolate the sink

    timeStepVector = dPS_[ipoint].time;
    // If the time is less or equal to the last specified cluster time
    // then do the interpolation, otherwise if the time is greater than
    // the last cluster time then just use the last value.
    double sink;
    size_t last_index = timeStepVector.size()-1;
    if (tn <= timeStepVector[last_index])
    {
      findTimeIndexes(tn,timeStepVector,timeIndexes);

      double sinkLo = dPS_[ipoint].sinkLifetime[timeIndexes[0]];
      double sinkHi = dPS_[ipoint].sinkLifetime[timeIndexes[1]];

      sink = sinkLo + (tn-dPS_[ipoint].time[timeIndexes[0]])*(sinkHi-sinkLo)
        /(dPS_[ipoint].time[timeIndexes[1]]-dPS_[ipoint].time[timeIndexes[0]]);
    }
    else
    {
      sink = dPS_[ipoint].sinkLifetime[last_index];
    }

    //If the point is close enough to a cluster point, use its value

    if(dist_to_point > 1.e-12)
    {
      w = pow((std::max(0.0,RofI-dist_to_point)/(RofI*dist_to_point)), p);

      pfunc += w*sink;

      total_distance += w;

	/*
      total_distance += (w = pow(dist_to_point,-p));
      pfunc += w * sink;
	*/
    }
    else
    {
      pfunc = sink;
      total_distance = 1.0;
      break;
    }
  }


  //This will happen if no points are used in the interpolation:
  //Localization causes this.  Simply set the interpolation to 0.0
  if(total_distance < 1.e-16)
  {
    pfunc = 0.0;
    total_distance = 1.0;
  }

  pn = pfunc/total_distance;

  return;

}


//------------------------------------------
//Find the time indexes
//------------------------------------------
void charon::interpolantMethod::findTimeIndexes(double localTime,
                                                std::vector<double>& timeStepVector,
                                                       int* indexes)
{


  indexes[0] = 0;
  indexes[1] = 1;

  size_t tSVsize = timeStepVector.size();

  if(localTime <= timeStepVector[0])return;

  for(size_t i=0 ; i < tSVsize-1 ; ++i)
  {

    if(localTime >= timeStepVector[i] &&
       localTime <= timeStepVector[i+1])
      return;
    else
    {
      indexes[0] += 1;
      indexes[1] += 1;
    }
  }

  if(localTime >= timeStepVector[tSVsize-1])
  {
    indexes[0] = tSVsize-2;
    indexes[1] = tSVsize-1;
    return;
  }

  indexes[0] = -1;
  indexes[1] = -1;

  return;

}



//-----------------------------------------------------------------------------
//Methods for the oneD linear interpolation Method interpolant class
//-----------------------------------------------------------------------------

//------------------------------------------
//interpolate at a point from the clusters
//by linear interpolation
//------------------------------------------
void charon::oneDLinearInterpolationMethod::interpolateToPoint(std::vector<dataPointSet> dPS_,
                                                               double xn, double /* yn */,
                                                               double /* zn */, double tn,
                                                               double& pn, double /* RofI */)
{

  //Given a point, loop over the data points and interpolate to the given point

  double x_dist;

  int timeIndexes[2];
  int spaceIndexes[2]={0,0};

  std::vector<double> timeStepVector;

  double lowSideDist=-1e101,highSideDist=1e101;
  bool lowSideFound=false,highSideFound=false;

  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);

  for(size_t ipoint=0 ; ipoint < dPS_.size() ; ++ipoint)
  {

    //Need to add support if model aligned with y axis.

    x_dist = dPS_[ipoint].xp - xn;

    if(x_dist < 0.0 && x_dist > lowSideDist)
    {
      lowSideDist = x_dist;
      spaceIndexes[0] = ipoint;
      lowSideFound = true;
    }
    if(x_dist > 0.0 && x_dist < highSideDist)
    {
      highSideDist = x_dist;
      spaceIndexes[1] = ipoint;
      highSideFound = true;
    }
  }


  if(!lowSideFound || !highSideFound)
  {
    pn = 0.0;
    return;
  }



  //Find the time indexes and Calculate the func on the low side

  timeStepVector = dPS_[spaceIndexes[0]].time;

  double sinkLo;
  double sinkHi;
  double lowSideSink;
  double highSideSink;

  // If the time is less or equal to the last specified cluster time
  // then do the interpolation, otherwise if the time is greater than
  // the last cluster time then just use the last value.
  size_t last_index = timeStepVector.size()-1;
  if (tn <= timeStepVector[last_index])
  {
    findTimeIndexes(tn,timeStepVector,timeIndexes);

    sinkLo = dPS_[spaceIndexes[0]].sinkLifetime[timeIndexes[0]];
    sinkHi = dPS_[spaceIndexes[0]].sinkLifetime[timeIndexes[1]];

    lowSideSink = sinkLo + (tn-dPS_[spaceIndexes[0]].time[timeIndexes[0]])*(sinkHi-sinkLo)
      /(dPS_[spaceIndexes[0]].time[timeIndexes[1]]-dPS_[spaceIndexes[0]].time[timeIndexes[0]]);
  }
  else
  {
    lowSideSink = dPS_[spaceIndexes[0]].sinkLifetime[last_index];
  }

  //Find the time indexes and Calculate the func on the high side

  timeStepVector = dPS_[spaceIndexes[1]].time;
  last_index = timeStepVector.size()-1;
  if (tn <= timeStepVector[last_index])
  {
    findTimeIndexes(tn,timeStepVector,timeIndexes);

    sinkLo = dPS_[spaceIndexes[1]].sinkLifetime[timeIndexes[0]];
    sinkHi = dPS_[spaceIndexes[1]].sinkLifetime[timeIndexes[1]];

    highSideSink = sinkLo + (tn-dPS_[spaceIndexes[1]].time[timeIndexes[0]])*(sinkHi-sinkLo)
      /(dPS_[spaceIndexes[1]].time[timeIndexes[1]]-dPS_[spaceIndexes[1]].time[timeIndexes[0]]);
  }
  else
  {
    highSideSink = dPS_[spaceIndexes[1]].sinkLifetime[last_index];
  }

  //Interpolate the value in between clusters

  double sink = lowSideSink + (xn-dPS_[spaceIndexes[0]].xp)*(highSideSink-lowSideSink)
    /(dPS_[spaceIndexes[1]].xp-dPS_[spaceIndexes[0]].xp);


  pn = sink;

  return;

}


//-----------------------------------------------------------------------------
//Methods for reading Cluster files
//-----------------------------------------------------------------------------
bool charon::clusterFiles::readFiles(std::vector<std::string> fileNames,
                                       std::vector<dataPointSet>& dPS_)
{

  bool success=false;

  for(size_t ifile=0 ; ifile<fileNames.size() ; ++ifile)
  {

    infile.open(fileNames[ifile].c_str());

    assert(!infile.fail());

    if(infile.fail())
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "OPEN FILE FAILED!!! ");
      return false;
    }

    // Read first line of Cluster file for x,y,z
    // and error if dimensions don't match
    std::getline(infile,line);
    std::istringstream firstline(line);

    switch(num_dim)
    {
      case 1:
        firstline>>xn;
        yn = 0.0;
        zn = 0.0;
        if(firstline>>trap && trap != 0.0)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Dimensions don't match in cluster file!!! ");
          return false;
        }
        break;

      case 2:
        firstline>>xn;
        if(!(firstline>>yn))
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Y value missing from cluster file!!! ");
          return false;
        }
        zn = 0.0;
        if(firstline>>trap && trap != 0.0)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Dimensions don't match in cluster file!!! ");
          return false;
        }
        break;

      case 3:
        firstline>>xn;
        if(!(firstline>>yn))
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Y value missing from cluster file!!! ");
          return false;
        }
        if(!(firstline>>zn))
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Z value missing from cluster file!!! ");
          return false;
        }
        if(firstline>>trap && trap != 0.0)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Dimensions don't match in cluster file!!! ");
          return false;
        }
        break;

      default:
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
              "Dimensions don't match in cluster file!!! ");
          return false;
    }


    dPS_temp.xp = xn;
    dPS_temp.yp = yn;
    dPS_temp.zp = zn;


    while(!infile.eof())
    {
      infile>>tn>>ps;
      // This avoids problems at the end of the file, otherwise
      // there will be an extra data point appended that is a
      // duplicate of the last data point
      if (infile)
      {
        dPS_temp.time.push_back(tn);
        dPS_temp.sinkLifetime.push_back(ps);
      }
    }

    dPS_.push_back(dPS_temp);

    dPS_temp.clear();

    infile.close();


  }

  success = true;

  return success;

}
