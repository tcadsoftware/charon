

#ifndef CHARON_EMPIRICAL_CONVOLUTION
#define CHARON_EMPIRICAL_CONVOLUTION

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Charon_PulseDamage_Spec.hpp"

//Forward Declarations

namespace charon {
  class charonSpline;
  class charonConvolute;
}

namespace charon {

  /**
   * \brief calculates a recombination rate based on empirical data
   */

  //------------------------------------------------------
  ///struct containing data for mudata
  //------------------------------------------------------
  class muData
  {

  public:

    ///constructor
    muData(std::string);

    std::vector <double>time,current,voltage;
    std::vector<std::vector <double> > muTable;
    ///creates a spline from a series of points
    void createSpline(double voltage);
    ///extracts a line segment of a surface
    void getTrace(double volt);
    ///gets the t* from mu table of empirical data
    double getTime(double mu);
    ///gets mu from a table of empirical data
    double getMu(double time);
    ///calculates the slope of the mu surface at a particular time
    double getDMuDt(double time);

    ///gets the time zero value of mu
    double getMuZero(){return muZeroValue;}

    ///gets the start time of mu where recombination starts
    double getMuStartTime(){return muStartTime;}

    ///Clears all mu data from the object
    void clear()
    {
      time.clear();
      for(std::vector<double>::size_type i=0 ; i<muTable.size() ; ++i)
        muTable[i].clear();
      muTable.clear();
      current.clear();
      voltage.clear();
    }

    ///destructor
    ~muData(){}

  private:

    std::vector<std::string> tableColumnNames;
    Teuchos::RCP<charonSpline> muSpline;
    Teuchos::RCP<charonSpline> muInverseSpline;

    std::vector<double> muTableTime,muTrace;
    std::vector<double> muZeroVector;
    double muZeroValue,muStartTime;

  };


  //--------------------------------------------------------
  /// The "main" class for computing empirical recombination
  /// coefficients from experimental data
  //--------------------------------------------------------
  class empiricalConvolution
  {


  public:

    /// Construct the object that calculates the convolution integral of recombination data
    empiricalConvolution(charon::PulseDamage_Spec const& damageSpec, std::string muFilename, bool pulseIsRate);

    empiricalConvolution(std::string muFilename);

    ///destructs the object that calculates the convolution integral of recombination data
    ~empiricalConvolution(){}

    ///computes the convolution integral of the product of the pulse and empirical data
    double computeNfpMu(const double time, const double timestep, double voltage);


  private:

    ///adds a pulse to the list of pulses that have occured
    void addMu(int count);

    ///calculates the value of mu at a given time
    void calculateMu(double currentTime, double timeStep);

    ///Stores an old value of mu
    void storeOldMu();

    ///restores an old value of mu
    void restoreOldMu();

   //Working Objects
    Teuchos::RCP<muData> muCalculator;

    Teuchos::RCP<charonConvolute> convolution;

    std::vector<double> mu,mutStar;
    std::vector<double> muOld;

    std::vector<double>::size_type pulseCount;
    double runningVoltage;
    double timeOfRecord;
    double NfpOfRecord;
    std::string pulseShape;

    PulseDamage_Spec damageSpec;
    bool pulseIsRate;
  };

//--------------------------------------------------------
/// This class contains general spline support for the empirical model (but is generic
//--------------------------------------------------------
class charonSpline
{

public:

  ///Charon spline constructor
  charonSpline(){}
  ///Charon spline constructor
  charonSpline(std::vector<double> x, std::vector<double> y);

  ///Creates a spline from a set of points.
  bool createSpline(std::vector<double> x, std::vector<double> y);

  ///Calculate y given x on the spline
  double evaluateSpline(double x);
  ///Calculate the slope of the spline
  double evaluateSplineDerivative(double x);

  ///calculate x given y on the spline
  double reverseEvaluateSpline(double y);

  //print out the spline
  void printSpline();

  ///destructs the spline
  ~charonSpline(){};


private:

  ///Performs a binary search
  int binarySearch(double x);
  std::vector<double> a_,b_,c_,d_,x_,y_;


};


  //--------------------------------------------------------
  /// This class contains integral convolution support for the empirical model
  //--------------------------------------------------------
class charonConvolute
{

public:

  ///construct convolution object
  charonConvolute(){}

  ///destruct convolution object
  ~charonConvolute(){};

  ///calculate a convolution integral from a generic set of data
  double convolve(bool pIR, double currentTime, double factor, std::vector<double> vector1,
                  std::vector<double> vector1Time, std::vector<double> vector2);

private:

};

}  //namespace charon

#endif //CHARON_EMPIRICAL_CONVOLUTION
