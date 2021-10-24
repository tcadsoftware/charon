
#ifndef Charon_XYCECLUSTER_HPP
#define Charon_XYCECLUSTER_HPP

#include "Xyce_config.h"
#include <N_CIR_Xyce.h>

#include <string>
#include <map>
#include <vector>

#include <boost/mpi.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/lexical_cast.hpp>

namespace Xyce {
namespace Circuit {

class clusterParameters
{

private:

  friend class  boost::serialization::access;

  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & clusterName;
    ar & clusterCanonicalName;
    ar & electronBC;
    ar & holeBC;
    ar & acceptorConcentration;
    ar & donorConcentration;
    ar & x_;
    ar & y_;
    ar & z_;
    ar & nameSet;
    ar & eBCSet;
    ar & hBCSet;
    ar & aCSet;
    ar & dCSet;
    ar & locationSet;
    ar & saveCoefficients;
    ar & onsetTime;
    ar & myProc;
    ar & deviceProc;
    ar & wonkyTime;
    ar & rank;
    ar & file;
    ar & updateMe;
  }


 public:

  clusterParameters() :
    clusterName(""),
    clusterCanonicalName(""),
    electronBC(0.0),
    holeBC(0.0),
    acceptorConcentration(0.0),
    donorConcentration(0.0),
    x_(0.0),
    y_(0.0),
    z_(0.0),
    nameSet(false),
    eBCSet(false),
    hBCSet(false),
    aCSet(false),
    dCSet(false),
    locationSet(false),
    saveCoefficients(false),
    onsetTime(0.0),
    myProc(0),
    deviceProc(0),
    wonkyTime(0.0),
    rank(-1),
    file(-1),
    updateMe(false)
  { }


  clusterParameters(std::string name, double eBC, double hBC, double aC, double dC,
                    double x=0.0, double y=0.0, double z=0.0) :
    clusterName(name),
    clusterCanonicalName(name),
    electronBC(eBC),
    holeBC(hBC),
    acceptorConcentration(aC),
    donorConcentration(dC),
    x_(x),
    y_(y),
    z_(z),
    nameSet(true),
    eBCSet(true),
    hBCSet(true),
    aCSet(true),
    dCSet(true),
    locationSet(false),
    saveCoefficients(false),
    onsetTime(0.0),
    myProc(0),
    deviceProc(0),
    wonkyTime(0.0),
    rank(-1),
    file(-1),
    updateMe(false)
  {}

  ~clusterParameters(){}//delete [] clusterName;}

    friend std::ostream& operator<<(std::ostream& os, clusterParameters& obj)
      {
        if(obj.locationSet)
          {
            double x,y,z;
            obj.getClusterLocation(x,y,z);
            os<<" The particle located at ("<<x<<","<<y<<","<<z<<"):"<<std::endl;
          }
        os<<" The name        is "<<obj.getName()<<std::endl;
        os<<" The electron BC is "<<obj.getElectronBC()<<std::endl;
        os<<" The hole     BC is "<<obj.getHoleBC()<<std::endl;
        os<<" The acceptor    is "<<obj.getAcceptorConc()<<std::endl;
        os<<" The donor       is "<<obj.getDonorConc()<<std::endl;
        os<<" The onset Time  is "<<obj.getOnsetTime()<<std::endl;
        os<<" I'm on proc      # "<<obj.getMyProc()<<std::endl;
        os<<" My rank is         "<<obj.getRank()<<std::endl;
        os<<" My file is         "<<obj.getFile()<<std::endl;
        return os;
      }


  std::string getName() { return clusterName;}
  std::string getCanonicalName() { return clusterCanonicalName;}
  std::string getElectronBC() { return boost::lexical_cast<std::string>(electronBC);}
  std::string getHoleBC() { return boost::lexical_cast<std::string>(holeBC);}
  double getElectronBCNumerical() { return electronBC;}
  double getHoleBCNumerical() { return holeBC;}
  std::string getAcceptorConc() { return boost::lexical_cast<std::string>(acceptorConcentration);}
  std::string getDonorConc() { return boost::lexical_cast<std::string>(donorConcentration);}
  double getAcceptorConcNumerical(){ return acceptorConcentration;}
  double getDonorConcNumerical(){ return donorConcentration;}
  void getClusterLocation(double& x, double& y, double& z){x=x_ ; y=y_ ; z=z_; locationSet=true;}
  bool getSaveCoefficients(){return saveCoefficients;}
  double getOnsetTime(){return onsetTime;}
  int getMyProc(){return myProc;}
  int getDeviceProc(){return deviceProc;}
  double getWonkyTime(){return wonkyTime;}
  int getRank(){return rank;}
  int getFile(){return file;}
  bool getUpdateMe(){return updateMe;}

  void setName(std::string name) {  clusterName = name; nameSet=true;}
  void setCanonicalName(std::string name) {  clusterCanonicalName = name; nameSet=true;}
  void setElectronBC(double eBC) {  electronBC = eBC; eBCSet=true;}
  void setHoleBC(double  hBC) {  holeBC = hBC; hBCSet=true;}
  void setAcceptorConc(double aC) {  acceptorConcentration = aC; aCSet=true;}
  void setDonorConc(double dC) {  donorConcentration = dC; dCSet=true;}
  void setClusterLocation(double x, double y, double z){x_=x ; y_=y ; z_=z;}
  void setSaveCoefficients(bool save){saveCoefficients=save;}
  void setOnsetTime(double oT){onsetTime = oT;}
  void setMyProc(int mP){myProc = mP;}
  void setDeviceProc(int dP){deviceProc = dP;}
  void setWonkyTime(double wT){wonkyTime = wT;}
  void setRank(int rnk){rank = rnk;}
  void setFile(int fle){file = fle;}
  void setUpdateMe(bool uM){updateMe = uM;}


private:

  //----------------------------------------------------------
  // WARNING!  DO NOT MODIFY THE CLASS MEMBERS WITHOUT ALSO
  // MODIFYING THE SERIALIZE METHOD ABOVE.
  // THIS IS NECESSARY TO ENSURE PROPER MPI COMMUNICATIONS
  //----------------------------------------------------------
  std::string clusterName,clusterCanonicalName;
  double electronBC, holeBC;
  double acceptorConcentration,donorConcentration;
  double x_,y_,z_;
  bool nameSet,eBCSet,hBCSet,aCSet,dCSet,locationSet,saveCoefficients;
  double onsetTime;
  int myProc,deviceProc;
  double wonkyTime;
  int rank,file;
  bool updateMe;

  //----------------------------------------------------------
  // WARNING!  DO NOT MODIFY THE CLASS MEMBERS WITHOUT ALSO
  // MODIFYING THE SERIALIZE METHOD ABOVE.
  // THIS IS NECESSARY TO ENSURE PROPER MPI COMMUNICATIONS
  //----------------------------------------------------------
};


  class CharonXyceInterface : public Simulator
  {

  public:

    CharonXyceInterface(MPI_Comm comm); // : Simulator(comm);

    //---------------------------------------------------------------------------
    // Function      : getData
    // Purpose       : Returns cluster recombination data
    // Special Notes :
    // Scope         : public
    // Creator       : Lawrence C Musson
    // Creation Date : 06/22/2015
    //---------------------------------------------------------------------------
    void getData(std::string clusterName, std::vector<double> & t, std::vector<double> & cc);

    //---------------------------------------------------------------------------
    // Function      : createClusterInputFile
    // Purpose       : create cluster input file for Xyce
    // Special Notes :
    // Scope         : public
    // Creator       : Lawrence C Musson
    // Creation Date : 06/22/2015
    //---------------------------------------------------------------------------
    std::string createClusterInputFile(std::string clusterTemplate, clusterParameters cP);


    //---------------------------------------------------------------------------
    // Function      : modifyClusterBCs
    // Purpose       : Modify the carrier BCs on the clusters
    // Special Notes :
    // Scope         : public
    // Creator       : Lawrence C Musson
    // Creation Date : 01/06/2016
    //---------------------------------------------------------------------------
    void modifyClusterBCs(std::string clusterName, double electronConcentration, double holeConcentration);

  };

}
}

#endif // Charon_XYCECLUSTER_HPP
