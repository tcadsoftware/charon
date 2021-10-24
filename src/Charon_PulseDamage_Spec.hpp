
#ifndef CHARON_PULSEDAMAGE_SPEC
#define CHARON_PULSEDAMAGE_SPEC

#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "Teuchos_ParameterList.hpp"

namespace charon {

class PulseDamage_Spec
{

public:

  enum Shape {Delta, Square, Gaussian, GaussianLog, File, Unknown};

  /** Default c'tor
   */
  PulseDamage_Spec(double const inpTimeSF,
                   Teuchos::ParameterList const& inpPL);

  static Shape shape(std::string const& pulseName);

  /** Return the time values of the delta spikes/pulse samples.
   */
  std::vector<double> const& times() {return pulse_times;}

  /** Return the time magnitude of the delta spikes/pulse samples.
   */
  std::vector<double> const& values() {return pulse_values;}

  /** Return the value of the pulse/sample at a specific index
   */
  double valueAtIndex(std::vector<double>::size_type index) {return pulse_values[index];}

  /** Return the number of delta spikes/samples representing the pulse
   */
  std::vector<double>::size_type numberOfPoints() {return pulse_times.size();}

  /** Return the pulse shape
   */
  Shape shape() {return pulse_shape;}

  /** Construct and return a string of break point times suitable for
   * the Rythmos parameter list.
   */
  std::string rythmosBPlist();

  /** return time of the pulse
   */
  double onsetTime(std::vector<double>::size_type pulseNumber){return pulse_times[pulseNumber];}

  /** Get the number of delta pulses/samples that have occurred up to
   * the input time
   */
  std::vector<double>::size_type grabPulses(double const& currenTime);


  /** Clear data storage
   */
  void clear(void)
  {
    pulse_times.clear();
    pulse_values.clear();
  }


protected:

  void checkRequiredParameters(char const* pulseType,
                               Teuchos::ParameterList const& inpParamList,
                               std::vector<std::string> const& paramNames);

  double getInterpolatedPulse(double time, std::vector<double> times, std::vector<double> pulses, std::string type);


  double time_sf;

  std::vector<double> pulse_times;
  std::vector<double> pulse_values;

  std::vector<std::string> const requiredParamNames;

  Shape pulse_shape;

  Teuchos::ParameterList inp_PL;

private:

  /** Disallow default ctor
   */
  PulseDamage_Spec() {;}

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

};

//----------------------------------------------------------
// delta pulse
//----------------------------------------------------------
class Delta_PulseDamage_Spec : public PulseDamage_Spec
{

public:

  ///Constructs a single delta pulse at a given time
  Delta_PulseDamage_Spec(Teuchos::ParameterList const& emParamList, double timeScaleFactor);

private:
  std::vector<std::string> const requiredParamNames;

};

//----------------------------------------------------------
/// square pulse
//----------------------------------------------------------
class Square_PulseDamage_Spec : public PulseDamage_Spec
{

public:

  ///Constructs a square wave pulse as a set of delta pulses
  Square_PulseDamage_Spec(Teuchos::ParameterList const& emParamList, double timeScaleFactor);

private:
  std::vector<std::string> const requiredParamNames;

};

//----------------------------------------------------------
/// gaussian pulse
//----------------------------------------------------------
class Gaussian_PulseDamage_Spec : public PulseDamage_Spec
{

public:

  ///constructs a gaussian shaped pulse with a set of delta pulses
  Gaussian_PulseDamage_Spec(Teuchos::ParameterList const& emParamList, double timeScaleFactor);

private:
  std::vector<std::string> const requiredParamNames;

};

//----------------------------------------------------------
/// gaussian log pulse
//----------------------------------------------------------
class GaussianLog_PulseDamage_Spec : public PulseDamage_Spec
{

public:

  ///constructs a gaussian shaped pulse in log time with a set of delta pulses
  GaussianLog_PulseDamage_Spec(Teuchos::ParameterList const& emParamList, double timeScaleFactor);

private:
  std::vector<std::string> const requiredParamNames;

};

//----------------------------------------------------------
/// tabulated pulse
//----------------------------------------------------------
class File_PulseDamage_Spec : public PulseDamage_Spec
{

public:

  ///reads in a set of pulses from a user-provided file
  File_PulseDamage_Spec(Teuchos::ParameterList const& emParamList, double timeScaleFactor);

private:
  std::vector<std::string> const requiredParamNames;

};

}

#endif
