
/*
This class manages all in-situ clusters, their parallel distribution and the
interface to the xyce library where the clusters are computed.
*/

#ifndef CHARON_CLUSTERMANAGER
#define CHARON_CLUSTERMANAGER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <fstream>
#include <mpi.h>
#include <new>
#include <sstream>
#include <string>
#include <vector>

// Boost
#include <boost/mpi.hpp>
#include <boost/mpi/exception.hpp>
#include <boost/serialization/serialization.hpp>

// Charon
#include "Charon_config.hpp"
#include "Charon_interp.hpp"
#include "Charon_Vector.hpp"
#include "Charon_XyceCluster.hpp"
#include "N_ERH_ErrorMgr.h"
#include "Xyce_config.h"

// Teuchos
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace charon
{

  //This class provides communication paths between the cluster manager and cluster interpolator
  class clusterManagerInterpolatorCommunicator
  {

  private:
    Teuchos::RCP<charon::clusterInterpolator> interpolator_;

  public:

    void registerInterpolator(Teuchos::RCP<charon::clusterInterpolator> interp);
    void fillInterpolatorWithData(std::vector<std::vector<double> > clusterDataTableTime,
                                  std::vector<std::vector<double> > clusterDataTableCoefficients);

  };


  //This class provides a packing for cluster data;
  class clusterDataSuitcase
  {

  private:

    friend class  boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
      ar & clusterName;
      ar & clusterNameList;
      ar & time;
      ar & coefficients;

      ar & clusterX;
      ar & clusterY;
      ar & clusterZ;
      ar & clusterFound;

      //The following are for communicating BCs--may not last
      ar & tempClusterNameList;
      ar & tempElectron;
      ar & tempHole;
      ar & tempAcceptor;
      ar & tempDonor;
      ar & tempdeviceProc;

    }

  public:

  //----------------------------------------------------------
  // WARNING!  DO NOT MODIFY THE CLASS MEMBERS WITHOUT ALSO
  // MODIFYING THE SERIALIZE METHOD ABOVE.
  // THIS IS NECESSARY TO ENSURE PROPER MPI COMMUNICATIONS
  //----------------------------------------------------------
    std::string clusterName;
    std::vector<std::string> clusterNameList;
    std::vector<double> time;
    std::vector<double> coefficients;
    std::vector<double> clusterX, clusterY, clusterZ;
    std::vector<bool> clusterFound;

    std::vector<std::string> tempClusterNameList;
    std::vector<double> tempElectron,tempHole;
    std::vector<double> tempAcceptor,tempDonor;
    std::vector<int> tempdeviceProc;
  //----------------------------------------------------------
  // WARNING!  DO NOT MODIFY THE CLASS MEMBERS WITHOUT ALSO
  // MODIFYING THE SERIALIZE METHOD ABOVE.
  // THIS IS NECESSARY TO ENSURE PROPER MPI COMMUNICATIONS
  //----------------------------------------------------------


    void packClusterNames(std::vector<std::string> names);
    void unpackClusterNames(std::vector<std::string> &names);


    //Add master list arrays for point data for the clusters
    //This is part of a KLUDGE and may get blasted later -- LCM
    void packClusterNames(std::vector<std::string> names, std::vector<double> x_,
                          std::vector<double> y_, std::vector<double> z_, std::vector<bool> cF);
    void unpackClusterNames(std::vector<std::string> &names, std::vector<double> &x_,
                          std::vector<double> &y_, std::vector<double> &z_, std::vector<bool> &cF);
    void packClusterBCs(std::vector<std::string> names_,  std::vector<double> electron_,
                        std::vector<double> hole_, std::vector<double> acceptor_, std::vector<double> donor_,
                        std::vector<int> deviceProc_);
    void unpackClusterBCs(std::vector<std::string> &names_,  std::vector<double> &electron,
                          std::vector<double> &hole, std::vector<double> &acceptor, std::vector<double> &donor,
                          std::vector<int> &deviceProc_);
    void clearTemps();
    //Add master list arrays for point data for the clusters
    //This is part of a KLUDGE and may get blasted later -- LCM

    void packClusterData(std::string name, std::vector<double> time, std::vector<double> coefficients);
    void unpackClusterData(std::string &name, std::vector<double> &time, std::vector<double> &coefficients);

  };


  //This class manages all things about the clusters
  class ClusterManager
  {

  public:

    ClusterManager();
    ClusterManager(Teuchos::ParameterList*);
    ~ClusterManager();

    void buildClusters();
    void obtainClusterParametersAndDistribute();
    void executeClusters(double stopTime);
    void registerInterpolatorClusters(std::vector<charon::dataPointSet> dp);
    void gatherClusterData();
    void setTemplateName(std::string tN){templateName = tN;}
    void transferToInterpolator();
    void registerInterpolator(Teuchos::RCP<charon::clusterInterpolator> interp){interpolator_ = interp;}
    void extractCarriersAndDopants(double time);
    void extractCarriersForUpdate(double time);
    void updateCarrierBCs(double time);
    void constructClusters();
    std::vector<double> extractBCBreakpoints();
    void addOnsetTimes(double oT);
    void broadcastOnsetTimesToMasters();
    void addClustersToTheRankAndFile(double time, int currentBreakPoint);

    Teuchos::RCP<clusterManagerInterpolatorCommunicator> cMIC;

    //This is part of a kludge to find data at cluster locations -- LCM
    void insertClusterLocationsIntoParameterList();


  private:  //member data

    int me, nprocs,clusterVerbosity;

    MPI_Group orig_group, my_group;
    MPI_Comm my_comm;
    int *ranks,sendbuf,recvbuf,my_rank;

    //create a vector of xyce objects
    std::vector<Xyce::Circuit::CharonXyceInterface*> xyce;
    //First the parameter objects;
    std::vector<Xyce::Circuit::clusterParameters> cP,cPMaster;
    std::vector<std::string> inputFiles;
    std::vector<std::string> clusterMasterNameList;
    std::string templateName;
    size_t numClusters,numClustersPerProc,leftoverClusters;
    std::vector<std::vector<double> > clusterDataTableTime,clusterDataTableCoefficients;
    Teuchos::RCP<charon::clusterInterpolator> interpolator_;

    //Some behavioral stuff
    bool saveBCHistory;

    //rank and file (space,time) refers to the array of clusters in the device
    //The clusters in a rank all share a common file (time)
    //The clusters in a file all share a common rank (spatial location)
    bool finitePulseOn;
    int rankSize,fileSize,currentFile;
    std::vector<double> pulseOnsetTimes,pulseMagnitudes;
    std::vector<double> breakPoints;
    bool breakPointsInitialized;
    std::vector<int> pulseBreakPoints;
    std::vector<std::string> deadRanks;
    bool useSingleFileClusters;

    //Add master list arrays for point data for the clusters
    //This is part of a KLUDGE and may get blasted later -- LCM
    Vector<std::string> clusterNames;
    Vector<double> clusterX, clusterY, clusterZ;
    Vector<double> clusterElectron, clusterHole, clusterAcceptor,clusterDonor;
    Vector<bool> clusterFound;
    int mostRecentBreakPoint;

    std::vector<std::ofstream*> convoData;  //DEBUGGING THING

  private:  //member functions

    void obtainClusterParameters();
    void readClusterLocations(std::string clusterLocationsFile);
    void addClusterParameterElement(std::vector<double> lineInputs);
    void distributeClusters();
    void createClusterInputFiles();
    void createClusters();
    void initializeClusters();
    void fillClusterDataTable(std::vector<clusterDataSuitcase> cDS);
    void broadcastAProcsData(const int, double, bool&);
    //Teuchos::RCP<Teuchos::ParameterList> clusterParameterList;
    Teuchos::ParameterList *clusterParameterList;
    std::vector<std::vector<double> > BCHistoryHoles,BCHistoryElectrons;
    std::vector<double> BCHistoryTime;
    void distributeNewClusters(int);
    void createNewClusters(int);
    void createNewClusterInputFiles(int);
    void initializeNewClusters(int);

    void convolveData(int irank, std::vector<double> & timeVector,
                      std::vector<double> & coefficientVector);
    double interpolateCoefficient(int irank, int ifile, int offset, double currentTime,
                                  std::vector<double> time, std::vector<double> coefficient);
    std::vector<double> getConvoTimesVector(double Tstart, double Tend);
    double interpolatePulse(int iwindow, double Tinterp);
    double integrate(int irank, int ifile, int offset, int pulseWindow, double currentTime,
                     std::vector<std::vector <double> > time,
                     std::vector<std::vector<double> > localCoefficient);
    double integrateWindow(int extent, std::vector<double> pulseVector, std::vector<double> responseVector,
                           std::vector<double> time);
    double getExtendedPulseTime();

  } ;  //end class ClusterManager


} //end namespace charon


#endif //CHARON_CLUSTERMANAGER



