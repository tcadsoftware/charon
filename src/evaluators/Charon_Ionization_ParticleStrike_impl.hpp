
#ifndef CHARON_IONIZATION_PARTICLE_STRIKE_IMPL_HPP
#define CHARON_IONIZATION_PARTICLE_STRIKE_IMPL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <typeinfo>

// Charon
#include "Charon_Ionization_ParticleStrike.hpp"
#include "Charon_Names.hpp"

// Panzer
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Ionization_ParticleStrike<EvalT, Traits>::
Ionization_ParticleStrike(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  int_rule_degree = ir->cubature_degree;
  


  //We always want the basis for this evaluator
  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  basis_name = basis->name();

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);
  num_dim = ir->dl_vector->dimension(2);

  // fields
  ionization_particle_strike_rate = MDField<ScalarT,Cell,Point>(n.field.ionization_particle_strike_rate,scalar);

  intrin_conc = MDField<const ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  C0 = scaleParams->scale_params.C0;

  generationRate=0.0;
  totalCharge = 0.0;

  //working things
  startX = p.get<double>("Start Point X");
  startY = p.get<double>("Start Point Y");
  startZ = p.get<double>("Start Point Z");
  endX = p.get<double>("End Point X");
  endY = p.get<double>("End Point Y");
  endZ = p.get<double>("End Point Z");
  strikeRadius = p.get<double>("Strike Radius");
  if(p.isParameter("Generation Rate"))
     generationRate = p.get<double>("Generation Rate");
  if(p.isParameter("Total Charge"))
     totalCharge = p.get<double>("Total Charge");
  startTime = p.get<double>("Start Time");
  endTime = p.get<double>("End Time");
  temporalWaveform = p.get<std::string>("Temporal Waveform");

  //Check for charge specification error

  if(fabs(generationRate) < 1.0e-20)
    generationRate = 0.0;
  if(fabs(totalCharge) < 1.0e-20)
    totalCharge = 0.0;

  if(fabs(generationRate) < 1.0e-20 && fabs(totalCharge) < 1.0e-20)
    {
      //Throw an error
      std::string msg = "Error in the parameter prescription for particle strike; generation rate and total charge cannot both be zero.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }

  if(fabs(generationRate) > 1.0e-20 && fabs(totalCharge) > 1.0e-20)
    {
      //Throw an error
      std::string msg = "Error in the parameter prescription for particle strike; cannot specify both a generation rate and total charge.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }

  if(startTime < 0)
    {
      //Throw an error
      std::string msg = "Error in the parameter prescription for particle strike; the start time must be greater than or equal to zero.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }

  if(startTime > endTime)
    {
      //Throw an error
      std::string msg = "Error in the parameter prescription for particle strike; the end time must be strictly greater than the start time.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }

  //Calculate the generation rate if a total charge was requested
  //assume the charged region is cylindrical

  double regionLength = sqrt((endX-startX)*(endX-startX) + (endY-startY)*(endY-startY)
                             + (endZ-startZ)*(endZ-startZ));

  //The following scale factor is present for two reasons
  //It converts coulombs to numbers of electrons for the generation rate
  //It converts cubit microns to cubic cm.
  double scaleFactor =  6.2415091e+18;

  double volumeConversion;

  if(num_dim==2)
    volumeConversion = 1.0e8;  //convert square microns to square cm
  else
    volumeConversion = 1.0e12;  //convert cubic microns to cubic cm

  double pi = acos(-1.0);

  double regionCrossSectionalArea = 0.0;
  if(num_dim == 2)
    regionCrossSectionalArea = 2.0*strikeRadius;
  else
    regionCrossSectionalArea = 2.0*pi*strikeRadius*strikeRadius;

  double strikeVolume = regionLength*regionCrossSectionalArea/volumeConversion;

  if(fabs(totalCharge) > 1.0e-20)
    {
      if(temporalWaveform == "Square")
        {
          generationRate = scaleFactor*totalCharge/(strikeVolume*(endTime-startTime));
        }
      else if(temporalWaveform == "Gaussian")
        {
          generationRate = scaleFactor*totalCharge/(strikeVolume);
        }
    }

  //Functions

  // Evaluated fields
  this->addEvaluatedField(ionization_particle_strike_rate);


  std::string name = "Ionization_Particle_Strike_Rate";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Ionization_ParticleStrike<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);

  //Always want the basis for this evaluator
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Ionization_ParticleStrike<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // get the current (present) time in [s]
  double curr_time = workset.time * t0;

  ScalarT scale_const = 1.0;

  if(C0 != 0) scale_const = t0/C0;

  std::vector<double> nodalX, nodalY, nodalZ;
  std::vector<int> nodesInThePath;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    nodalX.clear();
    nodalY.clear();
    nodalZ.clear();

    //Get the nodal coordinates
    for(int point=0 ; point < num_points ; ++point)
      {
        nodalX.push_back((workset.bases[basis_index])->basis_coordinates(cell,point,0));
        nodalY.push_back((workset.bases[basis_index])->basis_coordinates(cell,point,1));
	if(num_dim == 2)
	  nodalZ.push_back(0.0);
	else 
	  nodalZ.push_back((workset.bases[basis_index])->basis_coordinates(cell,point,2));
      }

    nodesInThePath.clear();
    nodesInThePath = amILitUp(nodalX, nodalY, nodalZ);

    if(nodesInThePath.size() == 0)  //If there are no nodes in the ionization trail, move on
      {
        for (int point = 0; point < num_points; ++point)
          ionization_particle_strike_rate(cell,point) = 0.0;
        continue;
      }

    if( curr_time < startTime
        || curr_time > endTime)
      {
        for (int point = 0; point < num_points; ++point)
          ionization_particle_strike_rate(cell,point) = 0.0;
        continue;
     }

    double timeFactor = 1.0;
    if(temporalWaveform != "Square")
      timeFactor = getTimeFactor(curr_time);


    for (int point = 0; point < num_points; ++point)
      {
        ionization_particle_strike_rate(cell,point) = -timeFactor*generationRate*scale_const;
      }

  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Ionization_ParticleStrike<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  //p->set<bool>("Is IP Set", false);

  p->set<double>("Start Point X",0.0);
  p->set<double>("Start Point Y",0.0);
  p->set<double>("Start Point Z",0.0);
  p->set<double>("End Point X",0.0);
  p->set<double>("End Point Y",0.0);
  p->set<double>("End Point Z",0.0);
  p->set<double>("Strike Radius",0.0);
  p->set<double>("Generation Rate",0.0);
  p->set<double>("Total Charge",0.0);
  p->set<double>("Start Time",0.0);
  p->set<double>("End Time",0.0);
  p->set<std::string>("Temporal Waveform","Square");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

///////////////////////////////////////////////////////////////////////////////
//
//  amILitUp()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<int>
Ionization_ParticleStrike<EvalT, Traits>::amILitUp(
   std::vector<double> nodalX, std::vector<double> nodalY, std::vector<double> nodalZ)
{

  std::vector<int> localNodes;

  if(num_dim == 3) {
    for(size_t inode=0 ; inode<nodalX.size() ; ++inode) {
      double lengthsq = (endX-startX)*(endX-startX) + (endY-startY)*(endY-startY) + 
	                (endZ-startZ)*(endZ-startZ);
      double radius_sq = strikeRadius*strikeRadius;
      double distsq = PointInCylinderTest(startX, startY, startZ,
                          endX, endY, endZ, lengthsq, radius_sq,
			  nodalX[inode], nodalY[inode], nodalZ[inode]);
      if(distsq >= 0.0) localNodes.push_back(inode);
     
    }
  } else { // 2D
    double x1,y1,x2,y2;  //,z1,z2;

    x1 = startX;
    y1 = startY;
    //z1 = startZ;
    x2 = endX;
    y2 = endY;
    //z2 = endZ;

    for(size_t inode=0 ; inode<nodalX.size() ; ++inode)
      {
	
	double xmin,xmax,ymin,ymax;
	double xnew=0.0,ynew=0.0,znew=0.0;
	double slope;

        /*
	cases:
	   horizontal (y1 = y2)
	   vertical (x1 = x2)
	   general:
	   compute perpendicular point to center of path
	   if point is off either end of path, then distance is to endpoint
	   else perpendicular point is the closest point on path
        */

	if (y1 == y2) {
	  xmin = std::min(x1,x2);
	  xmax = std::max(x1,x2);
	  if (nodalX[inode]< xmin)
	    {
	      xnew = xmin;
	      znew = 0.0;
	      if(xnew == x1)
		ynew = y1;
	      else
		ynew = y2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else if (nodalX[inode]> xmax)
	    {
	      xnew = xmax;
	      znew = 0.0;
	      if(xnew == x1)
		ynew = y1;
	      else
		ynew = y2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else
	    {
	      xnew = nodalX[inode];
	      ynew = y1;
	      znew = 0.0;
	      if(fabs(nodalY[inode]-y1) < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	} else if (x1 == x2) {
	  ymin = std::min(y1,y2);
	  ymax = std::max(y1,y2);
	  if (nodalY[inode]< ymin)
	    {
	      ynew = ymin;
	      znew = 0.0;
	      if(ynew == y1)
		xnew = x1;
	      else
		xnew = x2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else if (nodalY[inode]> ymax)
	    {
	      ynew = ymax;
	      znew = 0.0;
	      if(ynew == y1)
		xnew = x1;
	      else
		xnew = x2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else
	    {
	      xnew = x1;
	      ynew = nodalY[inode];
	      znew = 0.0;
	      if(fabs(nodalX[inode]-x1) < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	} else {
	  slope = (y2-y1) / (x2-x1);

	  double A,B,C,lam;
	  A = slope;
	  B = -1;
	  C = y1 - slope*x1;
	  lam = -(2*(A*nodalX[inode] + B*nodalY[inode] + C))/(A*A+B*B);
	  
	  xnew = 0.5*lam*A + nodalX[inode];
	  ynew = 0.5*lam*B + nodalY[inode];
	  znew = 0.0;

	  xmin = std::min(x1,x2);
	  xmax = std::max(x1,x2);
	  if (xnew < xmin)
	    {
	      xnew = xmin;
	      if(xnew == x1)
		ynew = y1;
	      else
		ynew = y2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else if (xnew > xmax)
	    {
	      xnew = xmax;
	      if(xnew == x1)
		ynew = y1;
	      else
		ynew = y2;
	      if(distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode])
		 < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	  else
	    {
	      double dist = distance(xnew,ynew,znew,nodalX[inode],nodalY[inode],nodalZ[inode]);
	      if(dist < strikeRadius)  //in the killzone!
		localNodes.push_back(inode);
	    }
	}
      }
  }

  return localNodes;
}


//-----------------------------------------------------
//distance calculator
//-----------------------------------------------------

template<typename EvalT, typename Traits>
double
Ionization_ParticleStrike<EvalT, Traits>::distance(double x1, double y1, double z1,
                                                   double x2, double y2, double z2)
{

  double dist = (x1-x2)*(x1-x2)
    + (y1-y2)*(y1-y2)
    + (z1-z2)*(z1-z2);

  return sqrt(dist);

}

  

template <typename EvalT, typename Traits> double
Ionization_ParticleStrike<EvalT, Traits>::PointInCylinderTest(
     double pt1_x, double pt1_y, double pt1_z, // first point to define cyl axis
     double pt2_x, double pt2_y, double pt2_z, // second point to define cyl axis
     double lengthsq, // cyl axis length squared
     double radius_sq, // cyl radius squared
     double testpt_x, double testpt_y, double testpt_z // test point
    ) 
{
 
  // vector d from point pt1 to point pt2
  double dx = pt2_x - pt1_x;	
  double dy = pt2_y - pt1_y;     
  double dz = pt2_z - pt1_z;

  // vector pd from point pt1 to test point.
  double pdx = testpt_x - pt1_x;		
  double pdy = testpt_y - pt1_y;
  double pdz = testpt_z - pt1_z;
  
  // Dot the d and pd vectors to see if point lies behind the 
  // cylinder cap at pt1
  double dot = pdx * dx + pdy * dy + pdz * dz;
  // If dot is less than zero the point is behind the pt1 cap.
  // If greater than the cylinder axis line segment length squared
  // then the point is outside the other end cap at pt2.
  if(dot < 0.0 || dot > lengthsq) 
    return -1.0; // test point outside the cylinder caps

  // calculate distance squared to the cylinder axis
  // dsq = pd*pd - (pd*cos(pd,d))^2 = pd*pd - dot*dot/ lengthsq
  double dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot/lengthsq;
  if( dsq > radius_sq )
    return -1.0;
  else 
    return dsq; // return distance squared to axis
}



//-----------------------------------------------------
//calculate time factor for non-square waveforms
//-----------------------------------------------------

template<typename EvalT, typename Traits>
double
Ionization_ParticleStrike<EvalT, Traits>::getTimeFactor(double currentTime)
{

  double timeFactor=1.0;

  if(temporalWaveform == "Gaussian")
    {

      //Assume that start time and end time of the temporal pulse waveform
      //spans 6 standard deviations of a normal waveform.
      //Assume that the width is +/- 3 standard deviations

      double width = endTime-startTime;
      double stdDev = width/6;

      double midPt = 0.5*(startTime+endTime);

      if(startTime < 0)
        {
          //Throw an error
          std::string msg = "Error in temporal Gaussian pulse prescription for particle strike; the pulse starts before time=0.\n";
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
        }

      double pi = acos(-1.0);
      double prefactor = 1.0/(stdDev*sqrt(2.0*pi));
      double argument = (currentTime-midPt)*(currentTime-midPt)/(2.0*stdDev*stdDev);

      timeFactor = prefactor*exp(-argument);
    }


  return timeFactor;

}

}

#endif //CHARON_IONIZATION_PARTICLE_STRIKE_IMPL_HPP
