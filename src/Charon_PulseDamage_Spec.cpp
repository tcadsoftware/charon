
#include <boost/algorithm/string.hpp>

#include "Panzer_Workset.hpp"

#include "Teuchos_Assert.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

#include "Charon_PulseDamage_Spec.hpp"

#include "Charon_ASCII_FileReader.hpp"

/****************************************************************************************************/
/****************************************************************************************************/
charon::PulseDamage_Spec::PulseDamage_Spec(double const inpTimeSF,
                                           Teuchos::ParameterList const& inpParams) :
  time_sf(inpTimeSF),
  inp_PL(inpParams)
{
  Teuchos::RCP<Teuchos::ParameterList> validPL = getValidParameters();

  inp_PL.validateParametersAndSetDefaults(*validPL);
}

/****************************************************************************************************/
/****************************************************************************************************/
Teuchos::RCP<Teuchos::ParameterList>
charon::PulseDamage_Spec::getValidParameters() const
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

  p->set<std::string>("material name","all");

  p->set<double>("eb x low", -std::numeric_limits<double>::max());
  p->set<double>("eb x high",  std::numeric_limits<double>::max());
  p->set<double>("eb y low", -std::numeric_limits<double>::max());
  p->set<double>("eb y high",  std::numeric_limits<double>::max());
  p->set<double>("eb z low", -std::numeric_limits<double>::max());
  p->set<double>("eb z high",  std::numeric_limits<double>::max());

  p->set<double>("thermal velocity", 0.0);
  p->set<double>("cross section", 0.0);

  p->set<std::string>("pulse data file", "");
  p->set<std::string>("mu data file", "");
  p->set<std::string>("pulse type", "");
  p->set<bool>("Is IP Set", false);

  p->set<double>("pulse start", 0.0);
  p->set<double>("pulse end", 0.0);
  p->set<double>("pulse magnitude", 0.0);

  p->set<int>("pulse resolution", 0);
  p->set<bool>("pulse is rate",false);

  p->set<double>("eb voltage override", -1.0);
  p->set<bool>("eb voltage override bool",false);

  p->set<double>("cb voltage override", -1.0);
  p->set<bool>("cb voltage override bool",false);

  p->set<std::string>("file pulse sampling scheme","all");


  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

/****************************************************************************************************/
/****************************************************************************************************/
charon::PulseDamage_Spec::Shape
charon::PulseDamage_Spec::shape(std::string const& pulseName)
{
  if (boost::iequals(pulseName, "delta"))
  {
    return charon::PulseDamage_Spec::Delta;
  }
  else if (boost::iequals(pulseName, "square"))
  {
    return charon::PulseDamage_Spec::Square;
  }
  else if (boost::iequals(pulseName, "gaussian"))
  {
    return charon::PulseDamage_Spec::Gaussian;
  }
  else if (boost::iequals(pulseName, "gaussianlog"))
  {
    return charon::PulseDamage_Spec::GaussianLog;
  }
  else if (boost::iequals(pulseName, "file"))
  {
    return charon::PulseDamage_Spec::File;
  }

  std::ostringstream err_msg;
  err_msg << "Unknown empirical damage pulse shape \"" << pulseName << "\" specified in input file";
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, err_msg.str());

  // Should never get to this point
  return charon::PulseDamage_Spec::Unknown;
}

/****************************************************************************************************/
/****************************************************************************************************/
std::string charon::PulseDamage_Spec::rythmosBPlist()
{

   TEUCHOS_TEST_FOR_EXCEPTION(pulse_times.size() == 0, std::runtime_error,
                             "This shouldn't happen in charon::PulseDamage_Spec::rythmosBPlist()");

  std::ostringstream bp_list_;

  bp_list_ << pulse_times[0]/time_sf;
  for (std::vector<double>::size_type i=1; i < pulse_times.size(); ++i)
  {
    bp_list_ << ", " << pulse_times[i]/time_sf;
  }

  return bp_list_.str();
}

/****************************************************************************************************/
/****************************************************************************************************/
std::vector<double>::size_type charon::PulseDamage_Spec::grabPulses(double const& currentTime)
{

  // If the current time is past the point of the last pulse, the number
  // of pulses is time.size
  if(currentTime > pulse_times[pulse_times.size()-1])
    return pulse_times.size();

  std::vector<double>::size_type pulses=0;
  for(std::vector<double>::size_type ipulses=0 ; ipulses < pulse_times.size() ; ++ipulses)
    if(currentTime >= pulse_times[pulses])
      ++pulses;

  TEUCHOS_ASSERT(pulses <= pulse_times.size());

  return pulses;

}

/****************************************************************************************************/
/****************************************************************************************************/
void charon::PulseDamage_Spec::checkRequiredParameters(char const* pulseName,
                                                       Teuchos::ParameterList const& inpParamList,
                                                       std::vector<std::string> const& paramList)
{
  std::vector<std::string>::const_iterator it = paramList.begin();

  std::ostringstream err_msg;
  unsigned int error_count = 0;
  for (;it != paramList.end(); ++it)
  {
    if (!inpParamList.isParameter(*it))
    {
      error_count++;
      err_msg << "ERROR[" << error_count << "]: You must specify \"" << *it
              << "\" in the \"Empirical Defect Recombination\" section of the "
              << "input file for a " <<  pulseName << " damage pulse" << std::endl;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION((error_count > 0), std::logic_error, err_msg.str());
}

/****************************************************************************************************/
/****************************************************************************************************/
charon::Delta_PulseDamage_Spec::Delta_PulseDamage_Spec(Teuchos::ParameterList const& emParamList,
                                                       double timeScaleFactor) :
  PulseDamage_Spec(timeScaleFactor,
                   emParamList),
  requiredParamNames({
      "pulse start",
      "pulse magnitude"
    })
{

  checkRequiredParameters("delta", inp_PL, requiredParamNames);

  double timeOn = inp_PL.get<double>(requiredParamNames[0]);
  double mag = inp_PL.get<double>(requiredParamNames[1]);

  if(timeOn < 0.0)
  {
    //Throw an error
    std::string msg = "Nice try, Einstein.  Can't start a pulse before time begins.\n";
    msg += " Try again. Time at pulse = "+std::to_string(timeOn);
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  pulse_times.push_back(timeOn);
  pulse_values.push_back(mag);

  pulse_shape = PulseDamage_Spec::shape("delta");
}

/****************************************************************************************************/
/****************************************************************************************************/
charon::Square_PulseDamage_Spec::Square_PulseDamage_Spec(Teuchos::ParameterList const& emParamList,
                                                         double timeScaleFactor) :
  PulseDamage_Spec(timeScaleFactor,
                   emParamList),
  requiredParamNames({"pulse start",
        "pulse end",
        "pulse magnitude",
        "pulse resolution"})

{
  checkRequiredParameters("square", inp_PL, requiredParamNames);

  double timeOn = inp_PL.get<double>(requiredParamNames[0]);
  double timeOff = inp_PL.get<double>(requiredParamNames[1]);
  double magnitude = inp_PL.get<double>(requiredParamNames[2]);
  int pulsePoints = inp_PL.get<int>(requiredParamNames[3]);

  if(timeOn >= timeOff)
  {
    //Throw an error
    std::string msg = "Error in square pulse prescription; there is a negative or zero-width pulse.\n";
    msg += " Try again. Time on = "+std::to_string(timeOn)+" ; time off = "+std::to_string(timeOff);
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  if(pulsePoints==1)
  {
    pulse_times.push_back(0.5*(timeOn+timeOff));
    pulse_values.push_back(magnitude);
  }
  else
  {
    double timeIncrement = (timeOff-timeOn)/(double)(pulsePoints-1);
    double currentTime = timeOn;
    for(int pulseTime=0 ; pulseTime<pulsePoints ; ++pulseTime)
    {
      pulse_times.push_back(currentTime);
      pulse_values.push_back(magnitude);
      currentTime += timeIncrement;
    }
  }

  pulse_shape = PulseDamage_Spec::shape("square");

}

/****************************************************************************************************/
/****************************************************************************************************/
charon::Gaussian_PulseDamage_Spec::Gaussian_PulseDamage_Spec(Teuchos::ParameterList const& emParamList,
                                                             double timeScaleFactor) :
  PulseDamage_Spec(timeScaleFactor,
                   emParamList),
  requiredParamNames({"pulse start",
        "pulse end",
        "pulse magnitude",
        "pulse resolution"})
{

  //Assume that the width is +/- 3 standard deviations
  checkRequiredParameters("gaussian", inp_PL, requiredParamNames);

  double timeOn = inp_PL.get<double>(requiredParamNames[0]);
  double timeOff = inp_PL.get<double>(requiredParamNames[1]);
  double magnitude = inp_PL.get<double>(requiredParamNames[2]);
  int pulsePoints = inp_PL.get<int>(requiredParamNames[3]);

  double width = timeOff-timeOn;
  double stdDev = width/6;

  double midPt = 0.5*(timeOn+timeOff);

  if(timeOn < 0)
  {
    //Throw an error
    std::string msg = "Error in Gaussian pulse prescription; the pulse starts before time=0.\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  if(pulsePoints <= 0)
  {
    //Throw an error
    std::string msg = "You must specify an \"pulse resolution\" greater than or equal to 1 for a gaussian damage pulse. \n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  if(pulsePoints == 1)
  {
    std::cout<<"WARNING:: Only one point specified for a gaussian pulse.  Pulse will be a delta at midpoint of wave."<<std::endl;
    pulse_times.push_back(midPt);
    pulse_values.push_back(magnitude);
  }
  else
  {
    double timeIncrement = (timeOff-timeOn)/(double)(pulsePoints-1);

    double currentTime = timeOn;

    for(int i=0 ; i<pulsePoints ; ++i)
    {
      double prefactor = 1.0;  ///(stdDev*sqrt(2.0*pi));
      double argument = (currentTime-midPt)*(currentTime-midPt)/(2.0*stdDev*stdDev);

      double waveform = prefactor*exp(-argument);

      pulse_times.push_back(currentTime);
      pulse_values.push_back(magnitude*waveform);

      currentTime += timeIncrement;
    }
  }

  pulse_shape = PulseDamage_Spec::shape("gaussian");

}

/****************************************************************************************************/
/****************************************************************************************************/
charon::GaussianLog_PulseDamage_Spec::GaussianLog_PulseDamage_Spec(Teuchos::ParameterList const& emParamList,
                                                                   double timeScaleFactor) :
  PulseDamage_Spec(timeScaleFactor,
                   emParamList),
  requiredParamNames({"pulse start",
        "pulse end",
        "pulse magnitude",
        "pulse resolution"})
{

  checkRequiredParameters("gauss log", inp_PL, requiredParamNames);

  double timeOn = inp_PL.get<double>(requiredParamNames[0]);
  double timeOff = inp_PL.get<double>(requiredParamNames[1]);
  double magnitude = inp_PL.get<double>(requiredParamNames[2]);
  int pulsePoints = inp_PL.get<int>(requiredParamNames[3]);

  //Assume that the width is +/- 3 standard deviations

  double width = log10(timeOff)-log10(timeOn);
  double stdDev = width/6;

  double midPt = 0.5*(log10(timeOn)+log10(timeOff));

  if(pulsePoints == 1)
  {
    std::cout<<"WARNING:: Only one point specified for a gaussian pulse.  Pulse will be a delta at midpoint of wave."<<std::endl;
    pulse_times.push_back(midPt);
    pulse_values.push_back(magnitude);
  }
  else
  {
    double timeIncrement = width/(double)(pulsePoints-1);

    double currentTime = log10(timeOn);

    for(int i=0 ; i<pulsePoints ; ++i)
    {
      double prefactor = 1.0;
      double argument = (currentTime-midPt)*(currentTime-midPt)/(2.0*stdDev*stdDev);

      double waveform = prefactor*exp(-argument);

      pulse_times.push_back(pow(10.0,currentTime));
      pulse_values.push_back(magnitude*waveform);

      currentTime += timeIncrement;
    }
  }

  pulse_shape = PulseDamage_Spec::shape("gaussianlog");

}

/****************************************************************************************************/
/****************************************************************************************************/
charon::File_PulseDamage_Spec::File_PulseDamage_Spec(Teuchos::ParameterList const& emParamList,
                                                     double timeScaleFactor) :
  PulseDamage_Spec(timeScaleFactor,
                   emParamList),
  requiredParamNames({"pulse data file"})
{
  checkRequiredParameters("pulse data file", inp_PL, requiredParamNames);

  std::string filename = inp_PL.get<std::string>(requiredParamNames[0]);

  std::string fileDataPulses = inp_PL.get<std::string>("file pulse sampling scheme");

  int pulsePoints = inp_PL.get<int>("pulse resolution");

  // This assumes a very simple format with time in the first column and
  // the pulse magnitude in the second. Anything else will likely break
  // with ugliness.

  ASCII_FileReader file_reader_(filename);

  std::vector<double> tmp_times,tmp_mags;

  for(size_t c=0; c < file_reader_.rowCount(); ++c)
  {
    tmp_times.push_back(file_reader_.dataValues(0)[c]);
    tmp_mags.push_back(file_reader_.dataValues(1)[c]);
  }

  if(fileDataPulses == "all")
    {
      pulse_times = tmp_times;
      pulse_values = tmp_mags;
    }
  else if(fileDataPulses == "linear")
    {
      double timeIncrement = (tmp_times[tmp_times.size()-1] - tmp_times[0])/(double)(pulsePoints-1);
      //Grab the first
      pulse_times.push_back(tmp_times[0]);
      pulse_values.push_back(tmp_mags[0]);

      double time = tmp_times[0];
      for(int pP=1 ; pP < pulsePoints-1 ; ++pP)
        {
          time += timeIncrement;
          pulse_times.push_back(time);
          pulse_values.push_back(getInterpolatedPulse(time,tmp_times,tmp_mags,fileDataPulses));
        }

      //Grab the last
      pulse_times.push_back(tmp_times[tmp_times.size()-1]);
      pulse_values.push_back(tmp_mags[tmp_mags.size()-1]);
    }
  else if(fileDataPulses == "log")
    {
      //Need to safeguard against taking the log of a t=0 initial time.  It's unlikely that a
      //pulse will ever start before a femtosecond.  Use that in place of zero.

      double logTime = log10(1.0e-12);
      if(tmp_times[0] > 0.0)
        logTime = log10(tmp_times[0]);

      double timeIncrement = (log10(tmp_times[tmp_times.size()-1]) - logTime)/(double)(pulsePoints-1);
      //Grab the first
      pulse_times.push_back(tmp_times[0]);
      pulse_values.push_back(tmp_mags[0]);

      for(int pP=1 ; pP < pulsePoints-1 ; ++pP)
        {
          logTime += timeIncrement;
          pulse_times.push_back(std::pow(10,logTime));
          pulse_values.push_back(getInterpolatedPulse(logTime,tmp_times,tmp_mags,fileDataPulses));
        }

      //Grab the last
      pulse_times.push_back(tmp_times[tmp_times.size()-1]);
      pulse_values.push_back(tmp_mags[tmp_mags.size()-1]);
    }
  else
    {
      //Throw an error
      std::string msg = " I don't have a valid file interpolant type in File_PulseDamage_Spec. \n";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
     }

  pulse_shape = PulseDamage_Spec::shape("file");

}


/****************************************************************************************************/
/****************************************************************************************************/

double charon::PulseDamage_Spec::getInterpolatedPulse(double time, std::vector<double> times,
                                                           std::vector<double> pulses, std::string type)
  {

    for(size_t index=1 ; index < times.size() ; index++)
      {

        if(type == "Linear")
          {
            if(time >= times[index-1] &&
               time <= times[index])
              {
                double offset = (time - times[index-1])/
                  (times[index] - times[index-1]);
                double pulseMag = pulses[index-1] + offset*(pulses[index] - pulses[index-1]);
                return pulseMag;
              }
          }

        if(type == "Log")
          {
            double timeim1;
            if(index-1 == 0)
              timeim1 = 1.0e-15;
            else
              timeim1 = times[index-1];

            if(time >= log10(timeim1) &&
               time <= log10(times[index]))
              {
                double offset = (time - log10(timeim1))/
                  (log10(times[index]) - log10(timeim1));
                double pulseMag = pulses[index-1] + offset*(pulses[index] - pulses[index-1]);
                return pulseMag;
              }
          }

      }

    //Throw an error
    std::string msg = "I reached the end of the list when trying to intperolate across tabulated pulses.  \n There must be an errir in  file_PulseDamageSpec\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);

  }

