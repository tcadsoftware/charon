
#include <fstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "Charon_empiricalConvolution.hpp"

//-----------------------------------------------------------------------------
//Methods for the empirical Convolution class
//-----------------------------------------------------------------------------

//--------------------------------------------------------
// Name:      empiricalConvolution::empiricalConvolution
// Purpose:   Constructor
// Notes:     This constructor takes a filename that
//            points mu tables and takes parameters for
//            a square wave pulse or Gaussian pulse
//
// Author: Lawrence C Musson
// Date: 5 February 2015
//--------------------------------------------------------

charon::empiricalConvolution::empiricalConvolution(charon::PulseDamage_Spec const& inp_damageSpec,
                                                   std::string muFilename, bool pIR) :
  pulseCount(0),
  runningVoltage(-100),
  timeOfRecord(0.0),
  NfpOfRecord(0.0),
  damageSpec(inp_damageSpec),
  pulseIsRate(pIR)
{
  muCalculator = Teuchos::rcp(new charon::muData(muFilename));
  convolution = Teuchos::rcp(new charon::charonConvolute());
}


//--------------------------------------------------------
// Name:      empiricalConvolution::copmuteNfpMu
// Purpose:   compute the product of the pulse and the
//            empirical factor for recombination
//
// Author: Lawrence C Musson
// Date: 5 February 2015
//--------------------------------------------------------


double charon::empiricalConvolution::computeNfpMu(const double time, const double timeStep, double voltage)
{

  //This is just stupid.
  if(timeStep==0.0)
    return NfpOfRecord;

  //If this calculation has already been done at this time, return the same NfpMu

  if(fabs(time-timeOfRecord) < 1.e-14) //already done here
      return NfpOfRecord;

  //If the timestep failed, the current time will be less than time of record
  //revert to an old mu and start again.

  if(time < timeOfRecord)
    restoreOldMu();

  //Save the old mu vector to handle a timestep failure

  storeOldMu();

  //Query the pulse to see if we need to add pulses.

  std::vector<double>::size_type countEmUp = damageSpec.grabPulses(time);

  //If we haven't yet hit the pulse, there's no point in going on.
  if(countEmUp == 0)
    return 0.0;

  if(countEmUp > pulseCount)
    addMu(countEmUp);

  //Check to see if we need to create a new spline
  if(voltage != runningVoltage)
    {
      muCalculator->createSpline(voltage);
      runningVoltage = voltage;
    }

  //Calculate the new mus
  calculateMu(time, timeStep);

  double NfpMu=0.0;

  //It makes no sense to compute a convolution integral if there is only a single
  //delta pulse.  Check that and operate as appropriate.
  if (damageSpec.shape() == charon::PulseDamage_Spec::Delta ||
      (damageSpec.shape() == charon::PulseDamage_Spec::File && damageSpec.numberOfPoints() == 1))
  {
    NfpMu = damageSpec.values()[0] * mu[0] * muCalculator->getMuZero();
  }
  else
  {

    NfpMu = convolution->convolve(pulseIsRate, time, muCalculator->getMuZero(),
                                  damageSpec.values(), damageSpec.times(), mu);
  }

  pulseCount = countEmUp;

  NfpOfRecord = NfpMu;
  timeOfRecord = time;

  return NfpMu;
}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------

//--------------------------------------------------------
// Name:      empiricalConvolution::addMu
// Purpose:   add new mus to the vector of mus
// Notes:
//
// Author: Lawrence C Musson
// Date: 5 February 2015
//--------------------------------------------------------


void charon::empiricalConvolution::addMu(int countEmUp)
{

  for(std::vector<double>::size_type newMu=0 ; newMu < countEmUp-pulseCount ; ++newMu)
    mu.push_back(1.0);

  //also add mutStar
  for(std::vector<double>::size_type newtStar=0 ; newtStar<countEmUp-pulseCount ; ++newtStar)
    mutStar.push_back(0.0);

}


//--------------------------------------------------------
// Name:      empiricalConvolution::storeOldMu
// Purpose:   store old mu as a hedge against time integration failure
// Notes:
//
// Author: Lawrence C Musson
// Date: 2 June 2015
//--------------------------------------------------------


void charon::empiricalConvolution::storeOldMu()
{

  muOld.clear();

  muOld = mu;
  //  for(int i=0 ; i<mu.size() ; ++i)
  //muOld.push_back(mu[i]);

}

//--------------------------------------------------------
// Name:      empiricalConvolution::restoreOldMu
// Purpose:   revert to old value of mu if time step failed
// Notes:
//
// Author: Lawrence C Musson
// Date: 2 June 2015
//--------------------------------------------------------


void charon::empiricalConvolution::restoreOldMu()
{

  mu.clear();

  mu=muOld;

  //  for(int i=0 ; i<mu.size() ; ++i)
  //muOld.push_back(mu[i]);

}


//--------------------------------------------------------
// Name:      empiricalConvolution::addMu
// Purpose:   add new mus to the vector of mus
// Notes:
//
// Author: Lawrence C Musson
// Date: 5 February 2015
//--------------------------------------------------------


void charon::empiricalConvolution::calculateMu(double currentTime, double timeStep)
{

  //This calculates the new value of mu*
  double logDt = log10(currentTime) - log10(currentTime-timeStep);

  for(std::vector<double>::size_type evolveMu=0 ; evolveMu < mu.size() ; ++evolveMu)
  {
    //If the current time minus the time the pulse hits is less than the
    //first time in the mu array, no evolution
    double timeSincePulse = currentTime - damageSpec.onsetTime(evolveMu);

    if(timeSincePulse < muCalculator->getMuStartTime())
    {
      mutStar[evolveMu] = 0.0; //Doesn't really matter what this is.
      mu[evolveMu] = 1.0;
      break;
    }

    double localMu = mu[evolveMu];
    mutStar[evolveMu] = muCalculator->getTime(localMu);
    double dMudt = muCalculator->getDMuDt(mutStar[evolveMu]);

    mu[evolveMu] += dMudt * logDt;
  }

}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//Methods for the muData class
//-----------------------------------------------------------------------------


//--------------------------------------------------------
// Name:      muData::muData
// Purpose:   Constructor
// Notes:     This constructor takes a filename to read a table of mu data
//
// Author: Lawrence C Musson
// Date: 5 February 2015
//--------------------------------------------------------

charon::muData::muData(std::string filename)
{

  std::ifstream muFile(filename);

  if(muFile.fail())
    {
      std::string msg = "Listen, buster.  Either give me a mu data filename I can open or move on.\n";
      msg += filename + " either doesn't exist or is not readable.";
      TEUCHOS_TEST_FOR_EXCEPTION(!muFile.is_open(), std::logic_error, msg);
    }

  //The first row in the file MUST have a name for each column.  No spaces in the name please.

  //The mu data will be arranged with current / voltage across the top row starting with
  //the first column will be time.  Subsequent columns will be the mu values.
  //This assumes that all points in the measured mu will be at identical times.  The table
  //can have any number of columns and rows, but it must be rectangular.

  //Read the top row as tokens so that we can count columns

  std::string firstLine;

  std::getline(muFile,firstLine);

  boost::char_separator<char> sep(" ");

  boost::tokenizer<boost::char_separator<char> > firstLineTokens(firstLine, sep);

  std::vector<std::string> tableColumnNames;

  for (const auto& t : firstLineTokens)
    {
      tableColumnNames.push_back(t);
    }

  //Now read the data
  std::vector<double> readVector;
  bool gottaReadCurrents=true;
  while(!muFile.eof())
    {
      if(gottaReadCurrents)  //First row after names is the current or voltage row
        {
          for(std::vector<std::string>::size_type columns=0 ; columns< tableColumnNames.size()-1 ; ++columns)
            {
              double value;
              muFile>>value;
              voltage.push_back(value);
            }
          gottaReadCurrents = false;
        }
      else
        {
          for(std::vector<std::string>::size_type columns=0 ; columns<tableColumnNames.size() ; ++columns)
            {
              double value;
              muFile>>value;
              if(columns==0)
                time.push_back(value);
              else
                readVector.push_back(value);
            }
          muTable.push_back(readVector);
          readVector.clear();
        }
    }

  //For some completely idiotic reason, the preceding code adds an extra
  //line to the time vector and mu tables.  I hate parsing files.  Delete the extra now.

  time.erase(time.end()-1);
  muTable.erase(muTable.end()-1);

  //if time[0] == 0, get the muZeroVector

  if(time[0] == 0.0)
    {
      for(size_t v=0 ; v<voltage.size() ; ++v)
        muZeroVector.push_back(muTable[0][v]);
    }


  muSpline = Teuchos::rcp(new charon::charonSpline());
  muInverseSpline = Teuchos::rcp(new charon::charonSpline());

}


//--------------------------------------------------------
// Name:      muData::createSpline
// Purpose:   builds splines to evaluate mu at various times for a given voltage
//
// Author: Lawrence C Musson
// Date: 11 February 2015
//--------------------------------------------------------

void charon::muData::createSpline(double volt)
{

  //First, call a method to get the mu trace and store the values in a vector
  getTrace(volt);

  //Create the spline on which all interpolations will be done.
  muSpline->createSpline(muTableTime,muTrace);
  muInverseSpline->createSpline(muTrace,muTableTime);

}


//--------------------------------------------------------
// Name:      muData::getTrace
// Purpose:   create trace at constant voltage in the mu table
// Notes:
//
// Author: Lawrence C Musson
// Date: 11 February 2015
//--------------------------------------------------------

void charon::muData::getTrace(double volt)
{

  std::vector<double> workingTime(time.size(),0.0);
  std::vector<double> workingMu(time.size(),0.0);

  //find the indexes of the columns in the mu table and compute a fraction...
  int start=0,end=0;
  for(size_t muColumn=0 ;  muColumn<voltage.size()-1 ; ++muColumn)
    {
    if(volt >= voltage[muColumn] && volt <= voltage[muColumn+1])
      {
        start = muColumn;
        end = muColumn+1;
      }
    }

  //The voltage falls outside the data supplied in the mu table.
  bool notInTheTable=false;
  if(start==0 && end==0)
    notInTheTable = true;
  if(notInTheTable)
  {
    std::ostringstream msg;
    msg << "ERROR: A voltage has been requested in the emprical model that is not spanned by the supplied data. ";
    msg << "The value was " << volt;
    TEUCHOS_TEST_FOR_EXCEPTION(notInTheTable, std::logic_error, msg.str());
  }


  //Grab the trace interpolating between voltage columns in the mu table
  //We want to modify this trace so that it starts at the first node away from
  //time=0 where the slope becomes non-zero.

  double voltFraction = (volt-voltage[start])/(voltage[end] - voltage[start]);

  //INterpolate the muZeroValue for this trace;
  muZeroValue = muZeroVector[start] + voltFraction * (muZeroVector[end] - muZeroVector[start]);

  //Extract the trace into working vectors
  for(size_t traceIndex=0 ; traceIndex<time.size() ; ++traceIndex)
    {

      workingTime[traceIndex] = time[traceIndex];

      double tempMu = muTable[traceIndex][start] +
        voltFraction * (muTable[traceIndex][end] - muTable[traceIndex][start]);

      workingMu[traceIndex] = tempMu;
    }

  int traceStart=0;
  //Find an acceptable place to start the trace;
  for(size_t findStart=0 ; findStart<workingTime.size() ; ++findStart)
    if(workingMu[findStart] != muZeroValue)
      {
        traceStart = findStart-1;
        break;
      }

  //This trace needs to be normalized to muZero

  //Check the smallest time, if it's zero, we can't take the log
  if(time[0] <= 0.0)
    time[0] = time[1]/2.0;

  for(size_t fillTrace=traceStart ; fillTrace<workingTime.size() ; ++fillTrace)
    {
      muTableTime.push_back(log10(time[fillTrace]));
      muTrace.push_back(workingMu[fillTrace]/muZeroValue);
    }

  muStartTime = time[traceStart];
}


//--------------------------------------------------------
// Name:      muData::getTime
// Purpose:   Use the inverse spline to get time from mu
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 February 2015
//--------------------------------------------------------

double charon::muData::getTime(double inputMu)
{

  if(inputMu > muTrace[0])
    return 0.0;

  double newTime = muInverseSpline->evaluateSpline(inputMu);

  return newTime;

}


//--------------------------------------------------------
// Name:      muData::getMu
// Purpose:   Use the spline to get mu at time inputTime
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 February 2015
//--------------------------------------------------------

double charon::muData::getMu(double inputTime)
{

  if(inputTime < muTableTime[0])
    return 1.0;

  double newMu = muSpline->evaluateSpline(inputTime);

  return newMu;

}


//--------------------------------------------------------
// Name:      muData::getDMuDt
// Purpose:   Use the spline to get dMuDt at time inputTime
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 February 2015
//--------------------------------------------------------

double charon::muData::getDMuDt(double inputTime)
{

  if(inputTime < muTableTime[0]) //prior to the first time in the spline, mu is constant in time
    return 0.0;

  double newDMuDt = muSpline->evaluateSplineDerivative(inputTime);

  return newDMuDt;

}


//-----------------------------------------------------------------------------
//Methods for the charonSpline class
//-----------------------------------------------------------------------------


//--------------------------------------------------------
// Name:      charonSpline::charonSpline
// Purpose:   Constructor
// Notes:     This constructor takes a filename to read a table of mu data
//
// Author: Lawrence C Musson
// Date: 10 February 2015
//--------------------------------------------------------

charon::charonSpline::charonSpline(std::vector<double> x, std::vector<double> y)
{
  createSpline(x,y);
}

//--------------------------------------------------------
// Name:      charonSpline::createSpline
// Purpose:   Constructor
// Notes:     This constructor takes a filename to read a table of mu data
//
// Author: Lawrence C Musson
// Date: 10 February 2015
//--------------------------------------------------------

bool charon::charonSpline::createSpline(std::vector<double> x, std::vector<double> y)
{

  //create local working vectors
  std::vector<double> c,l,u,z,h;

  x_.clear();
  y_.clear();

  int vectorSize = x.size()-1;

  x_.reserve(vectorSize);
  copy(x.begin(),x.end(),back_inserter(x_));

  y_.reserve(vectorSize);
  copy(y.begin(),y.end(),back_inserter(y_));

  //Spline parameters
  a_.resize(vectorSize);
  b_.resize(vectorSize);
  c_.resize(vectorSize+1);
  d_.resize(vectorSize);

  //local working vectors
  l.resize(vectorSize+1);
  u.resize(vectorSize+1);
  z.resize(vectorSize+1);
  h.resize(vectorSize+1);


  //Start the forward pass

  l[0] = 1.0;
  u[0] = 0;
  z[0] = 0;
  h[0] = x[1]-x[0];

  for(int i=1 ; i < vectorSize ; ++i)
    {
      h[i] = x[i+1]-x[i];
      l[i] = 2.0*(x[i+1]-x[i-1]) - h[i-1]*u[i-1];
      u[i] = h[i]/l[i];
      a_[i] = 3.0/h[i]*(y[i+1]-y[i]) - 3.0/h[i-1]*(y[i]-y[i-1]);
      z[i] = (a_[i] - h[i-1]*z[i-1])/l[i];
    }

  //Do the reverse pass

  l[vectorSize] = 1.0;
  z[vectorSize] = 0.0;
  c_[vectorSize] = 0.0;
  for(int j=vectorSize-1 ; j >= 0 ; --j)
    {
      c_[j] = z[j]-u[j]*c_[j+1];
      b_[j] = (y[j+1]-y[j])/h[j] -
               h[j]*(c_[j+1]+2.0*c_[j])/3.0;
      d_[j] = (c_[j+1]-c_[j])/(3.0*h[j]);
    }

  for(int i=0 ; i < vectorSize ; ++i)
    a_[i] = y[i];

  return true;

}


//--------------------------------------------------------
// Name:      charonSpline::evaluateSpline
// Purpose:   Constructor
// Notes:     Given x, return y
//
// Author: Lawrence C Musson
// Date: 10 February 2015
//--------------------------------------------------------

double charon::charonSpline::evaluateSpline(double xloc)
{

  int index=binarySearch(xloc);

  double xtmp = xloc - x_[index];

  double yReturn = a_[index]
    + b_[index] * xtmp
    + c_[index] * xtmp * xtmp
    + d_[index] * xtmp * xtmp * xtmp;

  return yReturn;

}

//--------------------------------------------------------
// Name:      charonSpline::evaluateSplineDerivative
// Purpose:   Constructor
// Notes:     Given x, return dy/dx
//
// Author: Lawrence C Musson
// Date: 10 February 2015
//--------------------------------------------------------

double charon::charonSpline::evaluateSplineDerivative(double xloc)
{

  int index=binarySearch(xloc);

  double xtmp = xloc - x_[index];

  double yReturn = b_[index]
    + 2.0 * c_[index] * xtmp
    + 3.0 * d_[index] * xtmp * xtmp;

  return yReturn;

}

//--------------------------------------------------------
// Name:      charonSpline::reverseEvaluate
// Purpose:   Constructor
// Notes:     Given Y, return x.
//
// Author: Lawrence C Musson
// Date: 12 February 2015
//--------------------------------------------------------

double charon::charonSpline::reverseEvaluateSpline(double yloc)
{

  //Find the closest index
  int index=0;
  bool foundOne=false;

  for(size_t splineIndex=1 ; splineIndex<y_.size() ; ++splineIndex)
    {
      //First check if the point is located right on the node.
      if(yloc == y_[splineIndex])
        {
          index = splineIndex;
          return x_[splineIndex];
        }
      if((yloc > y_[splineIndex-1] && yloc < y_[splineIndex]) ||
         (yloc < y_[splineIndex-1] && yloc > y_[splineIndex]))
        {
          index = splineIndex;
          foundOne=true;
          break;
        }
      index = y_.size()-1;
    }

  if(!foundOne) //Go with the nearest
    {

      double distance=1e10;;
      int minIndex=-1;
      for(size_t splineIndex=0 ; splineIndex<y_.size() ; ++splineIndex)
        {
          if(abs(yloc-y_[splineIndex]) < distance)
            {
              distance = abs(yloc-y_[splineIndex]);
              minIndex = splineIndex;
            }
        }
      index = minIndex;
    }

  double xReturn=0.0;

  //For now, just return the linearly interpolated location between nodes.

  double fraction = (yloc - y_[index-1])/(y_[index]-y_[index-1]);

  xReturn = x_[index-1] + fraction * (x_[index] - x_[index-1]);

  return xReturn;

}

//--------------------------------------------------------
// Name:      charonSpline::printSpline
// Purpose:   Constructor
// Notes:     print the values of a spline to the screen
//
// Author: Lawrence C Musson
// Date: 12 August 2017
//--------------------------------------------------------

void charon::charonSpline::printSpline()
{


  std::cout<<"-----------------"<<std::endl;
  std::cout<<"   X            Y"<<std::endl;
  for(size_t index=0 ; index<x_.size() ; ++index)
    std::cout<<std::scientific<<std::setprecision(6)<<x_[index]<<"      "<<y_[index]<<std::endl;

}

//--------------------------------------------------------
// Name:      charonSpline::binarySearch
// Purpose:   Constructor
// Notes:     Find an index of the spline given a x location
//
// Author: Lawrence C Musson
// Date: 10 February 2015
//--------------------------------------------------------

int charon::charonSpline::binarySearch(double xloc)
{

  int vectorSize = x_.size();
  int xmin,xmax;
  bool forward=true;
  if(x_[0] < x_[vectorSize-1])
    {
      xmin = 0;
      xmax = vectorSize-1;
    }
  else
    {
      xmin = vectorSize-1;
      xmax = 0;
      forward=false;
    }

  //if not on the curve, it's really kind of invalid for the spline.
  //return the closest endpoint.  it's up to user to handle this case.
  if(xloc <= x_[xmin])return xmin;
  if(xloc >= x_[xmax])return xmax;

  int imin = 0;
  int imax = x_.size()-1;

  while(imax > imin + 1)
    {
      int imid = imin + (imax - imin)/2;
      if(forward)
        if(x_[imid] > xloc)
          imax = imid;
        else
          imin = imid;
      else
        if(x_[imid] < xloc)
          imax = imid;
        else
          imin = imid;
    }
    return imin;

}

//--------------------------------------------------------
// Name:      charonConvolute::convolve
// Purpose:   does convolution
// Notes:
//
// Author: Lawrence C Musson
// Date: 17 February 2015
//--------------------------------------------------------

double charon::charonConvolute::convolve(bool pulseIsRate, double currentTime, double factor,
                                  std::vector<double> vector1, std::vector<double> time,
                                  std::vector<double> vector2)
{

  double convolution = 0.0;
  double timeFactor=0.0;
  int totalPoints=0;

  // Find out how many pulses have been seen to this point in time,
  // currentTime.
  if(currentTime >= time[time.size()-1])
    totalPoints = time.size();
  else
  {
    while (currentTime > time[totalPoints])
      totalPoints++;
  }

  if(totalPoints == 0) return 0.0;

  for(int convInt=1 ; convInt<totalPoints ; ++convInt)
    {

      //If the pulse waveform is a rate, multiply by timeSpan squared
      //This is to get a concentration of Frenkel pairs and the factor
      //for integration over the pulse by trapezoidal rule.
      //If the waveform is concentration, the factor is just the timespan
      if(pulseIsRate)
        timeFactor=(time[convInt] - time[convInt-1])*(time[convInt] - time[convInt-1]);
      else
        timeFactor=(time[convInt] - time[convInt-1]);

       convolution += 0.5*factor*(vector1[convInt-1]*vector2[convInt-1]
                                 + vector1[convInt]*vector2[convInt]) * timeFactor;


    }

  //final piece
  if(currentTime < time[time.size()-1])
    {
      double timeContractionFactor = (currentTime - time[totalPoints-1]) /
        (time[totalPoints] - time[totalPoints-1]);

       //If the pulse waveform is a rate, multiply by timeSpan squared
      //This is to get a concentration of Frenkel pairs and the factor
      //for integration over the pulse by trapezoidal rule.
      //If the waveform is concentration, the factor is just the timespan
      if(pulseIsRate)
        timeFactor=(currentTime - time[totalPoints-1])*(currentTime - time[totalPoints-1]);
      else
        timeFactor=(currentTime - time[totalPoints-1]);

      convolution += 0.5*factor*timeContractionFactor *
        (vector1[totalPoints-1] + vector1[totalPoints]) *
        (vector2[totalPoints-1]) * timeFactor;

    }


  return convolution;

}
