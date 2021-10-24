
#ifndef CHARON_CLUSTER_INTERPOLATOR
#define CHARON_CLUSTER_INTERPOLATOR

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <mpi.h>
#include <vector>

// Charon
#include "Charon_Vector.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

namespace charon {

  /**
   * \brief Interpolates arbitrarily placed cluster data to an arbitrary point
   */

  //------------------------------------------------------
  //struct containing data for a single cluster
  //------------------------------------------------------
  class dataPointSet
  {

  public:

    double xp,yp,zp;

    std::vector<double> time;

    std::vector<double> sinkLifetime;

    void clear(){time.clear();
      sinkLifetime.clear();}
  };

  //------------------------------------------------------
  //  interpolant Method Class -- A base class for various interp schemes
  //------------------------------------------------------
  class interpolantMethod
  {

  public:


    virtual void interpolateToPoint(std::vector<dataPointSet> dPS_,
                                    double xn, double yn,
                                    double zn, double tn,
                                    double& pn, double RofI) = 0;

    virtual ~interpolantMethod(){}

  protected:

    virtual void findTimeIndexes(double time,
                                 std::vector<double>& timeStepVector,
                                 int indexes[2]);



  };

  //------------------------------------------------------
  //  Read File Class -- A base class for various input file formats
  //------------------------------------------------------

  class readFileType
  {

  public:

    virtual bool readFiles(std::vector<std::string>,
                           std::vector<dataPointSet>& dPS_) = 0;

    virtual ~readFileType(){}

  private:


  };

  //--------------------------------------------------------
  // The "main" class for interpolations
  //--------------------------------------------------------


  class clusterInterpolator
  {


  public:

    clusterInterpolator() :
      influenceRadius(-1.0),
      influenceRadiusSet(false)
    {;}

    ~clusterInterpolator(){;}

    void interpolateToPoint(double xn, double yn,
                            double zn, double tn,
                            double& pn);

    bool readFiles(std::vector<std::string> fileNames);

    bool setMethod(std::string, float);
    bool setFileType(std::string);
    bool setInfluenceRadius(double iR);


    std::vector<dataPointSet> dPS;

    //These vectors are for moving BCs into the cluster manager
    Vector<std::string> clusterNames;
    Vector<double> clusterX, clusterY, clusterZ;
    Vector<double> clusterElectron, clusterHole, clusterAcceptor,clusterDonor;
    Vector<bool> clusterFound;

    void InitializeClusterBCVectors(Vector<std::string> clusterNames_);


  private:

    bool interpolantMethodFactory(std::string, const float&);
    bool readFileTypeFactory(std::string);
    double distance(double x1, double y1, double z1, double x2, double y2, double z2);

    double calculateInfluenceRadius();

    Teuchos::RCP<interpolantMethod> ipMethod;
    Teuchos::RCP<readFileType> rfType;

    size_t num_dim;

    double influenceRadius;
    bool influenceRadiusSet;

  };

  //------------------------------------------------------
  //  Shepard's Method Class
  //------------------------------------------------------
  class shepardsMethod : public interpolantMethod
  {

  public:

    shepardsMethod(size_t d, const float& shepardPneg):
      num_dim(d),
      p(shepardPneg)
    {;}

    virtual ~shepardsMethod() {}

    virtual void interpolateToPoint(std::vector<dataPointSet> dPS_,
                                    double xn, double yn,
                                    double zn, double tn,
                                    double& pn, double RofI);

  private:
    size_t num_dim;
    float p;
  };


  //------------------------------------------------------
  //  OneD Linear Interpolation Method Class
  //------------------------------------------------------
  class oneDLinearInterpolationMethod : public interpolantMethod
  {

  public:

    oneDLinearInterpolationMethod(){}

    virtual ~oneDLinearInterpolationMethod() {}

    virtual void interpolateToPoint(std::vector<dataPointSet> dPS_,
                                    double xn, double yn,
                                    double zn, double tn,
                                    double& pn, double RofI);

  };


  //------------------------------------------------------
  //  read cluster files
  //------------------------------------------------------
  class clusterFiles : public readFileType
  {

  public:

    clusterFiles(size_t d): num_dim(d)
    {;}

    virtual ~clusterFiles() {}

    virtual bool readFiles(std::vector<std::string>,
                           std::vector<dataPointSet>& dPS_);

  private:

  std::ifstream infile;
  std::string line;
  dataPointSet dPS_temp;

  double xn,yn,zn,tn,ps;

  double trap;

  size_t num_dim;
  };

}  //namespace charon

#endif //CHARON_CLUSTER_INTERPOLATOR
