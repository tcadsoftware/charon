
///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <inttypes.h>
#include <sstream>
#include <stdlib.h>
#include <ios>

// Boost
#include <boost/mpi/collectives.hpp>
#include <boost/tokenizer.hpp>
#include <boost/serialization/vector.hpp>

// Charon
#include "Charon_ClusterManager.hpp"
#include "Charon_Vector.hpp"

//-----------------------------------------------------------------------------
//Methods for the cluster manager class
//-----------------------------------------------------------------------------

//--------------------------------------------------------
// Name:      ClusterManager::ClusterManager
// Purpose:   Constructor
// Notes:     This constructor sets up the parallel environment
//            for the in-situ clusters
//
// Author: Lawrence C Musson
// Date: 7 June 2015
//--------------------------------------------------------

charon::ClusterManager::ClusterManager(Teuchos::ParameterList* cPL) :
  clusterVerbosity(0),
  saveBCHistory(false),
  finitePulseOn(false),
  fileSize(1),  //There'll always be at least one pulse
  currentFile(-1),   //of a series (in time) of clusters
  breakPointsInitialized(false),  //Need to make sure this is done only once
  useSingleFileClusters(false),
  mostRecentBreakPoint(-10)
{

  //The way clusters are computed in a parallel environment is
  //significantly different than the way devices are computed.  clusters
  //are always confined to a single core.  To prevent Xyce from
  //trying to split things up its own way, charon will control it.
  //First must split the mpi job into groups with new communicators of one-core
  //jobs each equal to the total number of cores in the charon job.  Then, the
  //in situ clusters will be distributed amongst them according to charon rules.

  //Save cluster parameters in the manager

  clusterParameterList = cPL;

  /* Extract the original group handle */

  MPI_Comm_rank(MPI_COMM_WORLD,&me);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  //Sub divide the procs into new groups/communicators

  ranks = new int;

  ranks[0] = me;
  sendbuf = me;

  MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
  MPI_Group_incl(orig_group, 1, ranks, &my_group);

  MPI_Comm_create(MPI_COMM_WORLD, my_group, &my_comm);
  MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, my_comm);
  MPI_Group_rank (my_group, &my_rank);

  //create the hook for the interpolator

  cMIC = Teuchos::rcp(new clusterManagerInterpolatorCommunicator);

  //extract some parameters
  if( clusterParameterList->isParameter("Save Cluster BC History"))
    saveBCHistory = clusterParameterList->get<bool>("Save Cluster BC History");

  if( clusterParameterList->isParameter("Single File Clusters"))
    useSingleFileClusters = clusterParameterList->get<bool>("Single File Clusters");

  //If this is going to be a finite width pulse problem, read in the series of impulses
  if(clusterParameterList->isParameter("Finite Pulse"))
    {
      finitePulseOn = true;
      std::string pulseFilename = clusterParameterList->get<std::string>("Finite Pulse");
      std::ifstream pulsesFile(pulseFilename);
      if(!pulsesFile.is_open())
        {
          std::string msg = "Hold it, bubba.  The Cluster Manager couldn't open the requested filename, ";
          msg += pulseFilename;
          TEUCHOS_TEST_FOR_EXCEPTION(!pulsesFile.is_open(), std::logic_error, msg);
        }

      double time, magnitude;

      while(pulsesFile >> time >> magnitude )
        {
          pulseOnsetTimes.push_back(time);
          pulseMagnitudes.push_back(magnitude);
        }

      pulsesFile.close();

    }


  //Get the verbosity level of the output from the cluster models
  if(clusterParameterList->isParameter("Cluster Verbosity"))
    {
      clusterVerbosity = clusterParameterList->get<int>("Cluster Verbosity");
    }


}


//--------------------------------------------------------
// Name:      ClusterManager::ClusterManager
// Purpose:   Destructor
// Notes:     Self explanatory
//
// Author: Lawrence C Musson
// Date: 7 June 2015
//--------------------------------------------------------

charon::ClusterManager::~ClusterManager()
{

  //Record BC History

  if(saveBCHistory && me==0)
    {
      //Loop over the fileSize (number of pulses) and produce output for each one.
      for(int ifile=0 ; ifile<fileSize ; ++ifile)
        {
          std::stringstream houtputFile;
          std::stringstream eoutputFile;

          houtputFile<<"BCHistoryHoles_"<<ifile<<".dat";
          eoutputFile<<"BCHistoryElectrons_"<<ifile<<".dat";


          std::ofstream ofileHoles(houtputFile.str());
          std::ofstream ofileElectrons(eoutputFile.str());

          ofileHoles<<"Time    ";
          ofileElectrons<<"Time    ";
          for(size_t icluster=0 ; icluster<clusterNames.size() ; ++icluster)
            {
              ofileHoles<<clusterNames[icluster]<<"    ";
              ofileElectrons<<clusterNames[icluster]<<"    ";
            }

          ofileHoles<<std::endl;
          ofileElectrons<<std::endl;

          size_t historyOffset = ifile*rankSize;

	  if(BCHistoryHoles.size() > 0)
          for(size_t bcUpdate=0 ; bcUpdate<BCHistoryHoles[historyOffset].size() ; ++bcUpdate)
            {
              ofileHoles<<BCHistoryTime[bcUpdate]<<"    ";
              ofileElectrons<<BCHistoryTime[bcUpdate]<<"    ";
              //for(size_t icluster=0 ; icluster<BCHistoryHoles.size() ; ++icluster)
              for(int icluster=0 ; icluster<rankSize ; ++icluster)
                {
                  ofileHoles<<BCHistoryHoles[icluster+historyOffset][bcUpdate]<<"    ";
                  ofileElectrons<<BCHistoryElectrons[icluster+historyOffset][bcUpdate]<<"    ";
                }
              ofileHoles<<std::endl;
              ofileElectrons<<std::endl;
            }
          ofileHoles.close();
          ofileElectrons.close();
        }

    }

  delete ranks;
  for ( size_t i=0 ; i<xyce.size() ; ++i)
    delete xyce[i];

}

//--------------------------------------------------------
// Name:      ClusterManager::buildClusters
// Purpose:   create cluster models
// Notes:
//
// Author: Lawrence C Musson
// Date: 7 June 2015
//--------------------------------------------------------

void charon::ClusterManager::buildClusters()
{

  //get cluster parameters by whatever means
  if(me==0)
    obtainClusterParameters();

  //distribute the parameter objects--and thus clusters--across procs
  distributeClusters();

  //create the clusters
  createClusters();

  //create the input netlists for xyce
  createClusterInputFiles();

  //initialize the clusters
  initializeClusters();

  //All done

}

//--------------------------------------------------------
// Name:      ClusterManager::buildClusters
// Purpose:   create cluster models
// Notes:    This is different from buildClusters in that we need
//           a two-part process of construction to couple the
//           clusters to the device.  The companion method to this one
//           is constructClusters
//
// Author: Lawrence C Musson
// Date: 7 June 2015
//--------------------------------------------------------

void charon::ClusterManager::obtainClusterParametersAndDistribute()
{

  //get cluster parameters by whatever means
  if(me==0)
    obtainClusterParameters();

  //distribute the parameter objects--and thus clusters--across procs
  distributeClusters();

}

//--------------------------------------------------------
// Name:      ClusterManager:constructClusters
// Purpose:   create cluster models
// Notes:    This is different from buildClusters in that we need
//           a two-part process of construction to couple the
//           clusters to the device.  The companion method to this one
//           is obtainClusterParametersAndDistribute.
//
// Author: Lawrence C Musson
// Date: 17 Februray 2016
//--------------------------------------------------------

void charon::ClusterManager::constructClusters()
{

  //create the clusters
  createClusters();

  //create the input netlists for xyce
  createClusterInputFiles();

  //initialize the clusters
  initializeClusters();

  //All done

}

//--------------------------------------------------------
// Name:      ClusterManager:addClustersToTheRankAndFile
// Purpose:   add new clusters
// Notes:    This method adds new clusters to the master list.
//           There is a tacit assumption that new clusters stack
//           the old but have different onset times.
//
// Author: Lawrence C Musson
// Date: 17 Februray 2016
//--------------------------------------------------------

void charon::ClusterManager::addClustersToTheRankAndFile(double time, int currentBreakPoint)
{

  int newFile=-1;
  bool startNewFile=false;

  //If useSingleFileClusters is true, the initial file of clusters gets used for the entire pulse
  if(useSingleFileClusters)
    return;

  //This is used in determining extent over which convolution integral is performed
  mostRecentBreakPoint = currentBreakPoint; //This is used in determining extent over which convolution integral is performed

  //This method potentially gets called frequently...most often as a no-op.
  //The first order of business is to determine if this really is a
  //stage where new clusters should be added and to which file they
  //should be added.  If not, bail.

  //N.B. the 0th pulse is what started this whole mess.  Start for the 1st
  for(size_t newCluster=1 ; newCluster < pulseOnsetTimes.size() ; ++newCluster)
    {
      //if(abs((time -  pulseOnsetTimes[newCluster])/pulseOnsetTimes[newCluster]) < 1.0e-12)
      if(currentBreakPoint == pulseBreakPoints[newCluster])
        {
          startNewFile = true;
          newFile=newCluster;
          break;
        }
    }

  if(!startNewFile)
    return;

  //Now some assumptions are made. The new file clusters are stacked on equivalent ranks
  //Push the new clusters onto the master and copy parameters

  for(int newCluster=0 ; newCluster<rankSize ; ++newCluster)
    {
      Xyce::Circuit::clusterParameters localCP=cPMaster[newCluster];

      cPMaster.push_back(localCP);
      cPMaster[rankSize*newFile+newCluster].setFile(newFile);
      cPMaster[rankSize*newFile+newCluster].setOnsetTime(pulseOnsetTimes[newFile]);
      std::string oldName = cPMaster[rankSize*newFile+newCluster].getCanonicalName();
      cPMaster[rankSize*newFile+newCluster].setCanonicalName(oldName);
      std::stringstream newName;
      newName<<oldName<<"_"<<newFile;
      cPMaster[rankSize*newFile+newCluster].setName(newName.str());
    }

  distributeNewClusters(newFile);

  createNewClusterInputFiles(newFile);

  createNewClusters(newFile);

  initializeNewClusters(newFile);

  //The new file size will be equal to newFile + 1
  fileSize = newFile+1;

  //Need to resize the cluster data tables for the larger system of clusters
  clusterDataTableTime.resize(numClusters);
  clusterDataTableCoefficients.resize(numClusters);

}


//--------------------------------------------------------
// Name:      ClusterManager::createNewClusters
// Purpose:   create new cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 28 March 2016
//--------------------------------------------------------

void charon::ClusterManager::createNewClusters(int newFile)
{


  //create the uninitialized clusters
  for(size_t i=0 ; i<cP.size() ; ++i)
    {
      if(cP[i].getFile() == newFile)
        {
          xyce.push_back( new Xyce::Circuit::CharonXyceInterface(my_comm) );
        }
    }

}

//--------------------------------------------------------
// Name:      ClusterManager::createClusterNewInputFiles
// Purpose:   create new cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 28 March 2016
//--------------------------------------------------------

void charon::ClusterManager::createNewClusterInputFiles(int newFile)
{


  //create the cluster input files for each of the clusters
  for(size_t i=0 ; i<cP.size() ; ++i)
    {
      if(cP[i].getFile() == newFile)
        inputFiles.push_back( xyce[i]->createClusterInputFile(templateName,cP[i]));
    }

}

//--------------------------------------------------------
// Name:      ClusterManager::buildClusters
// Purpose:   create cluster models
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::obtainClusterParameters()
{

  //NOTE: Even though there's no code in this method to reflect it, this method
  //is ONLY called on proc 0

  std::string clusterLocationsFile = clusterParameterList->get<std::string>("Cluster Locations File");
  //std::string clusterLocationsFile = clusterParameterList.get<std::string>("Cluster Locations File");

  if(clusterLocationsFile != "")
    {
      readClusterLocations(clusterLocationsFile);
      //oddly, the master will be one cluster too big
      cPMaster.erase(cPMaster.end()-1);
    }

  //Take care of saving coefficients
  bool saveCoefficients;
  if( clusterParameterList->isParameter("Save Cluster Coefficients"))
    {
      saveCoefficients = clusterParameterList->get<bool>("Save Cluster Coefficients");
      //saveCoefficients = clusterParameterList.get<bool>("Save Cluster Coefficients");

      for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
        cPMaster[icluster].setSaveCoefficients(saveCoefficients);
    }

  double wonkyTime;
  if(clusterParameterList->isParameter("Cluster Wonky Time"))
    {
      wonkyTime = clusterParameterList->get<double>("Cluster Wonky Time");
      for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
        cPMaster[icluster].setWonkyTime(wonkyTime);
    }

  //Set the rank, file and onset time in each of these clusters.  This is the first shot, so
  //the file of each cluster will be 0

  double pulseOnset=0;
  if(finitePulseOn)
    pulseOnset = pulseOnsetTimes[0];

  if(finitePulseOn)
    for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
      {
        cPMaster[icluster].setRank(icluster);
        cPMaster[icluster].setFile(0);
        cPMaster[icluster].setOnsetTime(pulseOnset);
      }

  BCHistoryElectrons.resize(cPMaster.size());
  BCHistoryHoles.resize(cPMaster.size());

}


//--------------------------------------------------------
// Name:      ClusterManager::readClusterLocationsFile
// Purpose:   read locations of clusters and potentially dopes and carriers
// Notes:
//
// Author: Lawrence C Musson
// Date: 12 August 2015
//--------------------------------------------------------

void charon::ClusterManager::readClusterLocations(std::string clusterLocationsFile)
{

  //NOTE: Even though there's no code in this method to reflect it, this method
  //is ONLY called on proc 0

  std::ifstream clusterLocations(clusterLocationsFile);

  if(clusterLocations.fail())
    {
      std::string msg = "Listen, buster.  Either give me a cluster locations filename I can open or move on.\n";
      msg += clusterLocationsFile + " either doesn't exist or is not readable.";
      TEUCHOS_TEST_FOR_EXCEPTION(!clusterLocations.is_open(), std::logic_error, msg);
    }

  //The first line determines the number of columns in the file.  This file must be space delimited

  std::string firstLine;

  std::getline(clusterLocations,firstLine);

  boost::char_separator<char> sep(" ");

  boost::tokenizer<boost::char_separator<char> > firstLineTokens(firstLine, sep);

  std::vector<double> lineInputs;

  for (const auto& t : firstLineTokens)
    {
      lineInputs.push_back(strtod(t.c_str(),NULL));
    }

  size_t numberOfColumns = lineInputs.size();

  addClusterParameterElement(lineInputs);
  lineInputs.clear();

  while(!clusterLocations.eof())
    {
      for(size_t iLine = 0 ; iLine<numberOfColumns ; ++iLine)
        {
          double value;
          clusterLocations>>value;
          lineInputs.push_back(value);
        }
      addClusterParameterElement(lineInputs);
      lineInputs.clear();
    }


  //Take care of saving coefficients
  bool saveCoefficients;
  saveCoefficients = clusterParameterList->get<bool>("Save Cluster Coefficients");
  //saveCoefficients = clusterParameterList.get<bool>("Save Cluster Coefficients");

  for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
    cPMaster[icluster].setSaveCoefficients(saveCoefficients);
}


//--------------------------------------------------------
// Name:      ClusterManager::addClusterParameterElement
// Purpose:   create an element of a vector of cluster parameter objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 12 August 2015
//--------------------------------------------------------

void charon::ClusterManager::addClusterParameterElement(std::vector<double> lineInputs)
{

  //NOTE: Even though there's no code in this method to reflect it, this method
  //is ONLY called on proc 0

  //There are only a few possibilities for the length of line inputs.

  //lineInputs.zise() = 2,  2D, positions only
  //lineInputs.zise() = 3,  3D, positions only
  //lineInputs.zise() = 6,  2D, positions and dopes and carriers
  //lineInputs.zise() = 7,  3D, positions and dopes and carriers

  double x=0,y=0,z=0;
  double acceptors=0,donors=0;
  double holes=0,electrons=0;

  switch(lineInputs.size())
    {

    case 2:  //2D positions only


        x=lineInputs[0];
        y=lineInputs[1];
        z=0.0;
        {
          Xyce::Circuit::clusterParameters cPLocal("",electrons,holes,acceptors,donors,x,y,z);
          cPMaster.push_back(cPLocal);
        }
        break;

    case 3:  //3D positions only


      x=lineInputs[0];
      y=lineInputs[1];
      z=lineInputs[2];
      {
        Xyce::Circuit::clusterParameters cPLocal("",electrons,holes,acceptors,donors,x,y,z);
        cPMaster.push_back(cPLocal);
      }
      break;


    case 6:  //2D positions, dopes and carriers

      x=lineInputs[0];
      y=lineInputs[1];
      z=0.0;
      electrons=lineInputs[2];
      holes=lineInputs[3];
      acceptors=lineInputs[4];
      donors=lineInputs[5];
      {
        Xyce::Circuit::clusterParameters cPLocal("",electrons,holes,acceptors,donors,x,y,z);
        cPMaster.push_back(cPLocal);
      }
      break;


    case 7:  //3D positions, dopes and carriers

      x=lineInputs[0];
      y=lineInputs[1];
      z=lineInputs[2];
      electrons=lineInputs[3];
      holes=lineInputs[4];
      acceptors=lineInputs[5];
      donors=lineInputs[6];
      {
        Xyce::Circuit::clusterParameters cPLocal("",electrons,holes,acceptors,donors,x,y,z);
        cPMaster.push_back(cPLocal);
      }
      break;

    default:
      std::string msg = "There was an invalid number of columns in the cluster locations file.  Try again.  2,3,6 or 7.\n";
      TEUCHOS_TEST_FOR_EXCEPTION(false, std::logic_error, msg);

    }


}


//--------------------------------------------------------
// Name:      ClusterManager::distributeClusters
// Purpose:   distribute cluster models across procs
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::distributeClusters()
{

  //boost::mpi::environment env;
  boost::mpi::communicator world;

  //Get all cluster parameters on proc 0, then scatter

  int *recvcounts,*sendcounts;

  sendcounts = new int[nprocs];
  recvcounts = new int[nprocs];

  numClusters=0;
  numClustersPerProc=0;
  leftoverClusters=0;


  if(me==0)
    {
      //Scatter them to all procs

      numClusters=cPMaster.size();
      numClustersPerProc=numClusters/nprocs;

      leftoverClusters = numClusters % nprocs;

      for (int i=0 ; i<nprocs ; ++i)
        sendcounts[i] = numClustersPerProc;

      for(size_t i=0 ; i<leftoverClusters ; ++i)
        sendcounts[i]++;

      //at this point, create machine generated names for each of the clusters

      //Create a unique name for each cluster

      std::string myName;
      int iclusterNumber=0;
      for(size_t icluster=0 ; icluster<numClustersPerProc ; ++icluster)
	{
	  for(int iproc=0 ; iproc<nprocs ; ++iproc)
	    {
	      std::stringstream myNameStream;
	      myNameStream<<iproc<<"Cluster"<<icluster;
	      cPMaster[iclusterNumber].setCanonicalName(myNameStream.str());
	      myNameStream<<"_0";
	      cPMaster[iclusterNumber].setName(myNameStream.str());
	      cPMaster[iclusterNumber].setMyProc(iproc);
	      iclusterNumber++;
	    }
	}


      //And now the leftovers
      for(size_t iproc=0 ; iproc<leftoverClusters ; ++iproc)
        {
          std::stringstream myNameStream;
          myNameStream<<iproc<<"Cluster"<<numClustersPerProc;
          cPMaster[iclusterNumber].setCanonicalName(myNameStream.str());
          myNameStream<<"_0";
          cPMaster[iclusterNumber].setName(myNameStream.str());
          cPMaster[iclusterNumber].setMyProc(iproc);
          iclusterNumber++;
        }
    }


  MPI_Bcast(&numClusters,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numClustersPerProc,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&leftoverClusters,1,MPI_INT,0,MPI_COMM_WORLD);

  //This method is executed only the first time clusters are created //set the ranksize now.
  rankSize = numClusters;

  //Pass a master list of cluster CANONICAL names to all procs.  This is necessary for book-keeping
  clusterMasterNameList.resize(numClusters);

  //It's also turned out necessary to have a master parameter list on each proc
  if(me != 0)
    cPMaster.resize(numClusters);
  for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
      broadcast(world,cPMaster[icluster],0);


  //Add master list arrays for point data for the clusters
  //This is part of a KLUDGE and may get blasted later -- LCM
  clusterNames.resize(numClusters);
  clusterX.resize(numClusters);
  clusterY.resize(numClusters);
  clusterZ.resize(numClusters);
  clusterElectron.resize(numClusters);
  clusterHole.resize(numClusters);
  clusterAcceptor.resize(numClusters);
  clusterDonor.resize(numClusters);
  clusterFound.resize(numClusters);


  if(me==0)
    {
      for(size_t i=0 ; i<numClusters ; ++i)
        clusterMasterNameList[i] = cPMaster[i].getCanonicalName();

      for(size_t i=0 ; i<numClusters ; ++i)
        {
          double cx,cy,cz;
          cPMaster[i].getClusterLocation(cx,cy,cz);
          clusterX[i] = cx;
          clusterY[i] = cy;
          clusterZ[i] = cz;
          clusterFound[i] = false;
        }

      //ALL PART OF THE KLUDGE - LCM
      for(size_t i=0 ; i<numClusters ; ++i)
        clusterNames[i] = cPMaster[i].getName();

    }

  //Above this is a kludge -- LCM

  clusterDataSuitcase cDS;

  if(me == 0)
    cDS.packClusterNames(clusterMasterNameList);

  broadcast(world,cDS,0);

  if(me != 0)
    cDS.unpackClusterNames(clusterMasterNameList);

  //At this point, we can create the clustersize dimension of the double vector that stores cluster data
  clusterDataTableTime.resize(numClusters);
  clusterDataTableCoefficients.resize(numClusters);

  //Let each proc know how many clusters they're getting
  int *mycount;
  mycount = new int;
  MPI_Scatter(sendcounts,1,MPI_INT,mycount,1,MPI_INT,0,MPI_COMM_WORLD);

  cP.resize(*mycount);

  //The serialization brought by boost is really nice and I want to use it.  But boost's
  //collective communication is lacking.  Deal with it now.  First, do scatters through
  //a loop equal to the smallest number of clusters on any one node.  Will deal with the
  //leftovers later.
  for(size_t i=0 ; i<numClustersPerProc ; ++i)
    {
      int offset = i*nprocs;
      scatter(world,&cPMaster[offset],cP[i],0);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  int numReqs;
  if(me==0)
    numReqs = leftoverClusters+1;
  else
    numReqs = 1;


  boost::mpi::request *req;
  req = new boost::mpi::request[numReqs];
  //And now the leftovers with point-to-point comms
  //post the receives
  for (size_t i=0 ; i<leftoverClusters ; ++i)
    if((size_t)me==i)
          req[0] = world.irecv(0,i,cP[*mycount-1]);

  MPI_Barrier(MPI_COMM_WORLD);

  //do the sends
  if(me == 0)
    for (size_t i=0 ; i<leftoverClusters ; ++i)
      req[i+1] = world.isend(i,i,cPMaster[nprocs*numClustersPerProc+i]);

  boost::mpi::wait_all(req, req+numReqs);

  delete [] sendcounts;
  delete [] recvcounts;
  delete [] req;
  delete mycount;

  
}


//--------------------------------------------------------
// Name:      ClusterManager::distributeNewClusters
// Purpose:   distribute the added file of cluster models across procs
// Notes:    Going to try doing this without parallel communication
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::distributeNewClusters(int newFile)
{

  Xyce::Circuit::clusterParameters cPLocal;

  //newFile is is the current index of the pulses vector, i.e. the index of the file
  int filloutClusters=nprocs-leftoverClusters;   //equal to nprocs - leftover from previous fileadd
  if(leftoverClusters == 0)
    filloutClusters = 0;
  numClusters = cPMaster.size();

  //Need to reset leftover clusters
  leftoverClusters = numClusters % nprocs;
  //Need to reset num clusters per proc
  numClustersPerProc = numClusters/nprocs;


  //create local clusters from fillout clusters first

  int clusterCounter=newFile*rankSize;  //This equals the number of clusters that already exist
  //By coincidence, the fillout clusters quantity is equal to the first proc to get the new clusters

  //  std::cout<<me<<" fillouts "<<clusterCounter<<"    "<<filloutClusters<<std::endl;

  if(filloutClusters != 0)
    for(int icluster=filloutClusters ; icluster<nprocs ; ++icluster)  //if the leftovers are zero, this no-ops
      {
        if(me == icluster)
          {
            cPLocal = cPMaster[clusterCounter];
            cP.push_back(cPLocal);
            std::cout<<"  Distributing fillout clusters "<<std::endl;
            std::cout<<me<<"    "<<clusterCounter<<std::endl;
          }
        clusterCounter++;
      }

  int localNumClustersPerProc = (rankSize - filloutClusters)/nprocs;
  int startIndex = clusterCounter;

  for(int icluster=0 ; icluster < localNumClustersPerProc ; ++icluster)
    {
      int masterIndex = startIndex + me + icluster*nprocs;  //This is a little formula to get the
                                                            //index of the CPMaster that belongs
                                                            //on a given proc
      cP.push_back(cPMaster[masterIndex]);
    }

  clusterCounter = numClusters - leftoverClusters;

  //And finish them off

  for(size_t icluster=0 ; icluster<leftoverClusters ; ++icluster)
    {
      if(me == (int)icluster)
        {
          cP.push_back(cPMaster[clusterCounter]);
        }
      ++clusterCounter;
    }

  //Resize the BC history vectors to accommodate new clusters
  if(me==0)
    {
      BCHistoryElectrons.resize(cPMaster.size());
      BCHistoryHoles.resize(cPMaster.size());
    }

}



//--------------------------------------------------------
// Name:      ClusterManager::createClusters
// Purpose:   create cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::createClusters()
{


  //create the uninitialized clusters
  for(size_t i=0 ; i<cP.size() ; ++i)
    xyce.push_back( new Xyce::Circuit::CharonXyceInterface(my_comm) );

}

//--------------------------------------------------------
// Name:      ClusterManager::createClusterInputFiles
// Purpose:   create cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::createClusterInputFiles()
{


  //create the cluster input files for each of the clusters
  for(size_t i=0 ; i<cP.size() ; ++i)
    inputFiles.push_back( xyce[i]->createClusterInputFile(templateName,cP[i]));
}


//--------------------------------------------------------
// Name:      ClusterManager::initializeClusters
// Purpose:   initialize cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::initializeClusters()
{

  //Initialize each of the clusters
  for(size_t i=0 ; i<xyce.size() ; ++i)
    {
      int myargc = 4;
      char** myargv = new char*[myargc];

      for(int j=0 ; j<myargc ; ++j)
        myargv[j] = new char[100];

      strncpy(myargv[0],"a.out",100);
      strncpy(myargv[1],inputFiles[i].c_str(),100);

      strncpy(myargv[2],"-l",100);
      std::string options = inputFiles[i]+".log";
      strncpy(myargv[3],options.c_str(),100);

      bool bsuccess = xyce[i]->initialize(myargc, myargv);

      if(!bsuccess)
        std::cout<<" WARNING!!!  "<<cP[i].getName()<<" failed to initialize properly!"<<std::endl;

      for (int j=0; j < myargc; ++j)
        delete[] myargv[j];

      delete[] myargv;

    }
}


//--------------------------------------------------------
// Name:      ClusterManager::initializeNewClusters
// Purpose:   initialize new cluster model objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::initializeNewClusters(int newFile)
{

  //Initialize each of the clusters
  for(size_t i=0 ; i<xyce.size() ; ++i)
    {

      //Only initialize new clusters
      if(cP[i].getFile() != newFile)
        continue;

      int myargc = 4;
      char** myargv;
      myargv = new char*[myargc];

      for(int j=0 ; j<myargc ; ++j)
        myargv[j] = new char[100];

      strncpy(myargv[0],"a.out",100);
      strncpy(myargv[1],inputFiles[i].c_str(),100);

      strncpy(myargv[2],"-l",100);
      std::string options = inputFiles[i]+".log";
      strncpy(myargv[3],options.c_str(),100);

      bool bsuccess = xyce[i]->initialize(myargc, myargv);

      if(!bsuccess)
        std::cout<<" WARNING!!!  "<<cP[i].getName()<<" failed to initialize properly!"<<std::endl;

      delete [] myargv;

    }
}


//--------------------------------------------------------
// Name:      ClusterManager::executeClusters
// Purpose:   execute clusters
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::ClusterManager::executeClusters(double stopTime)
{

  boost::mpi::communicator world;

  bool success;
  double completedTime;

  //N.B.  When this method is called, the time will be in device time, not cluster time.
  //      We need to offset the stopTime by the cluster onsetTime in order to keep things in sync

  std::cout<<std::setprecision(8);

  if(me==0 && clusterVerbosity > 0)
    {
      std::cout<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<std::endl;
      std::cout<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<std::endl;
      std::cout<<"||       /======== //         //      //  /======  ============   /=======  /=====/      /======    ||"<<std::endl;
      std::cout<<"||     //         //         //      // //              //      //        //     //    //           ||"<<std::endl;
      std::cout<<"||    //         //         //      // //              //      //        //     //    //            ||"<<std::endl;
      std::cout<<"||   //         //         //      //  /=====/        //      //=====// // ====/      /=====/       ||"<<std::endl;
      std::cout<<"||  //         //         //      //        //       //      //        //    \\\\            //       ||"<<std::endl;
      std::cout<<"|| //         //         //      //        //       //      //        //      \\\\          //        ||"<<std::endl;
      std::cout<<"|| /========  /========  /======/   ======/        ==       /======= //        \\\\  ======/          ||"<<std::endl;
      std::cout<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<std::endl;
      std::cout<<"||"<<std::endl;
      std::cout<<"|| Executing all clusters to device time "<<stopTime<<"s"<<std::endl;
      std::cout<<"||"<<std::endl;

    }

  //Find the longest cluster name for formatted output
  size_t nameLength=0;
  for(size_t iname=0 ; iname<cPMaster.size() ; ++iname)
    if(cPMaster[iname].getName().length() > nameLength)
      nameLength = cPMaster[iname].getName().length();

  std::vector<std::string> deadRanksLocal;
  deadRanksLocal.resize(0);
  for(size_t i=0 ; i<xyce.size() ; ++i)
    {
      bool skipThis=false;
      for(size_t icluster=0 ; icluster<deadRanks.size() ; ++icluster)
        if(cP[i].getCanonicalName() == deadRanks[icluster])
          {
            if(clusterVerbosity > 1)
              std::cout<<"||  *** Skipping Rank "<<cP[i].getName()<<" because rank "<<deadRanks[icluster]<<" is dead. ***"<<std::endl;
            skipThis = true;
            break;
          }

      if(skipThis) continue;

      success = xyce[i]->simulateUntil(stopTime-cP[i].getOnsetTime(),completedTime);

      if(!success && clusterVerbosity > 1)
        {
          std::cout<<
            " WARNING!! Cluster "<<cP[i].getName()<<" has failed to converge.\n Adding it to the dead soldierslist."
                                 <<std::endl;
          std::cout<<cP[i]<<std::endl;
          deadRanksLocal.push_back(cP[i].getCanonicalName());
        }
      else
        if(clusterVerbosity > 1)
	  {
	    std::string message = " has completed execution to cluster time ";
	    std::cout<<"||  "<<std::setw(nameLength)<<std::left<<cP[i].getName()<<std::setw(message.length())<<message<<completedTime<<"s on processor number "<<me<<std::endl;
	  }
    }

  //Need to wait for all processes to be done with cluster obligations before moving on.
  MPI_Barrier(MPI_COMM_WORLD);

  //Broadcast dead ranks to the whole world

  //Make some dummy arrays just to take advantage of the packClusterNames serializer

  std::vector<double> dummyX,dummyY,dummyZ;
  std::vector<bool> dummyFC;

  charon::clusterDataSuitcase deadcdP;

  for(int iproc=0 ; iproc<nprocs ; ++iproc)
    {
      if(iproc == me)
        deadcdP.packClusterNames(deadRanksLocal,dummyX,dummyY,dummyZ,dummyFC);

      broadcast(world,deadcdP,iproc);

      std::vector<std::string> r_deadRanksLocal;

      deadcdP.unpackClusterNames(r_deadRanksLocal,dummyX,dummyY,dummyZ,dummyFC);

      for(size_t icluster=0 ; icluster<r_deadRanksLocal.size() ; ++icluster)
        deadRanks.push_back(r_deadRanksLocal[icluster]);

    }


  if(me==0 && clusterVerbosity > 0)
    {
      std::cout<<"||"<<std::endl;
      std::cout<<"||  Returning to device integration. "<<std::endl;
      std::cout<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<std::endl;
      std::cout<<"-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"<<std::endl;
    }
}

//--------------------------------------------------------
// Name:      ClusterManager::gatherClusterData
// Purpose:   gather data from the clusters and store in a table
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 July 2015
//--------------------------------------------------------

void charon::ClusterManager::gatherClusterData()
{

  boost::mpi::communicator world;

  charon::clusterDataSuitcase cDS;
  std::vector<charon::clusterDataSuitcase> cDSGlobal;

  cDSGlobal.resize(nprocs);

  //This starts with a loop over the number of local clusters

  for(size_t icluster=0 ; icluster<numClustersPerProc ; ++icluster)
    {

      std::vector<double> time,coefficients;
      xyce[icluster]->getData(cP[icluster].getName(),time,coefficients);

      cDS.packClusterData(cP[icluster].getName(),time,coefficients);

      all_gather(world,cDS,cDSGlobal);

      //Add this set of data to the tables
      fillClusterDataTable(cDSGlobal);
    }

  cDSGlobal.resize(1);

  //Now need to take care of the leftovers
  int offset=numClustersPerProc;

  for(size_t icluster=0 ; icluster<leftoverClusters ; ++icluster)
    {
      if((size_t)me==icluster)
        {
          std::vector<double> time,coefficients;
          xyce[icluster]->getData(cP[offset].getName(),time,coefficients);
          cDSGlobal[0].packClusterData(cP[offset].getName(),time,coefficients);
        }
      broadcast(world,cDSGlobal,icluster);
      fillClusterDataTable(cDSGlobal);
    }
}

//--------------------------------------------------------
// Name:      ClusterManager::fillClusterDataTable
// Purpose:   gather data from the clusters and store in a table
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 July 2015
//--------------------------------------------------------

void charon::ClusterManager::fillClusterDataTable(std::vector<charon::clusterDataSuitcase> cDS)
{

  for(size_t icluster=0 ; icluster<cDS.size() ; ++icluster)
    {
      std::string name=cDS[icluster].clusterName;
      int columnNumber=-1;
      bool columnFound = false;
      //Loop over all clusters here.

      for(size_t column=0 ; column<cPMaster.size() ; ++column)
        {
          if(name == cPMaster[column].getName())
            {
              columnNumber = column;
              columnFound = true;
              break;
            }
        }
      if(!columnFound)
        std::cout<<" ERROR!! ERROR!! ERROR!!  Missing matching cluster names when trying to extract data. "<<name<<std::endl;
      if(columnNumber >= (int)clusterDataTableTime.size())
        std::cout<<" ERROR!! ERROR!! ERROR!!  The cluster data tables have not been sized correctly "<<name<<std::endl;

      clusterDataTableTime[columnNumber].clear();
      clusterDataTableCoefficients[columnNumber].clear();
      clusterDataTableTime[columnNumber] = cDS[icluster].time;
      clusterDataTableCoefficients[columnNumber] = cDS[icluster].coefficients;

      //Adjust for wonky time.  I.e. often the beginning timesteps of a cluster are messy because
      //of the impulsive nature of the model.  Ignore them by setting all coefficients prior to
      //"wonky time" (model user input).

      double wonkyTime = clusterParameterList->get<double>("Cluster Wonky Time");
      double wonkyValue=0;
      for (size_t i=0 ; i<clusterDataTableCoefficients[columnNumber].size() ; ++i)
        {
          if(clusterDataTableTime[columnNumber][i] > wonkyTime)
            {
              wonkyValue = clusterDataTableCoefficients[columnNumber][i];
              break;
            }
        }

      for (size_t i=0 ; i<clusterDataTableCoefficients[columnNumber].size() ; ++i)
        {
          if(clusterDataTableTime[columnNumber][i] < wonkyTime)
            clusterDataTableCoefficients[columnNumber][i] = wonkyValue;
          else
            break;
        }

    }

}

//--------------------------------------------------------
// Name:      ClusterManager::transferToInterpolator
// Purpose:   move cluster data into interpolator
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 July 2015
//--------------------------------------------------------

void charon::ClusterManager::transferToInterpolator()
{

  interpolator_->dPS.clear();

  //interpolator_->dPS.resize(rankSize);
  //interpolator_->dPS.resize(rankSize - deadRanks.size());   //accommodate dead ranks. //Handle later


  if(deadRanks.size())
    {
      std::cout<<"  Removed "<<deadRanks.size()<<" from the ranks on "<<me<<std::endl;

      std::cout<<"Dead ranks are: "<<std::endl;
      for(size_t i=0 ; i<deadRanks.size() ; ++i)
        std::cout<<me<<"    "<<i<<"    "<<deadRanks[i]<<std::endl;
    }

  //Compile a vector of all ranks that wil be sent to the interpolator

  std::vector<bool> amIAlive;

  for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
    {
      amIAlive.push_back(true);
      for(size_t idead=0 ; idead<deadRanks.size() ; ++idead)
        {
          if(cPMaster[icluster].getCanonicalName() == deadRanks[idead])
            {
              amIAlive[icluster] = false;
              continue;
            }
        }
    }


  //interpolator_->dPS.resize(cPMaster.size());

  //for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
  size_t dPSSize=0;
  for(int icluster=0 ; icluster<rankSize ; ++icluster)
    {
      if(!amIAlive[icluster]) continue;
      ++dPSSize;
      interpolator_->dPS.resize(dPSSize);
      interpolator_->dPS[dPSSize-1].xp = clusterX[icluster];
      interpolator_->dPS[dPSSize-1].yp = clusterY[icluster];
      interpolator_->dPS[dPSSize-1].zp = clusterZ[icluster];
    }

  std::vector<double> timeVector, coefficientVector;

 //for(size_t i=0 ; i<clusterDataTableTime.size() ; ++i)
  dPSSize = 0;

  if((int)interpolator_->dPS.size() != rankSize)
    std::cout<<" Something is wrong.  dPS is the wrong size in the interpolator. "<<
      rankSize<<"   "<<interpolator_->dPS.size()<<std::endl;

  for(int irank=0 ; irank<rankSize ; ++irank)
    {
      //If the rank is dead, move on
      if(!amIAlive[irank]) continue;

      //compute convolution integral over rank irank
      convolveData(irank,timeVector,coefficientVector);

      //interpolator_->dPS[i].time = clusterDataTableTime[i];
      //interpolator_->dPS[i].sinkLifetime = clusterDataTableCoefficients[i];

      //interpolator_->dPS[irank].time = timeVector;
      //interpolator_->dPS[irank].sinkLifetime = coefficientVector;

      ++dPSSize;
      interpolator_->dPS[dPSSize-1].time = timeVector;
      interpolator_->dPS[dPSSize-1].sinkLifetime = coefficientVector;

      if(me==0)
        {
          if(timeVector.size() != coefficientVector.size())
            std::cout<<" There is something wrong.  Sizes of time and coefficient don't match "<<timeVector.size()<<"    "<<coefficientVector.size()<<std::endl;
	  //DEBUG STUFF

	  //if(me==0)
	  //{	      
	  //for(size_t icluster=0 ; icluster<cPMaster.size() ; ++icluster)
	  //{
	  //	  std::string name = cPMaster[icluster].getCanonicalName();
	  	  std::string name = cPMaster[irank].getCanonicalName();
		  std::string filename = "ConvolvedData."+name+".dat";
		  convoData.push_back(new std::ofstream(filename));
		  //}
		  //}
          for(size_t ivec=0 ; ivec<timeVector.size() ; ++ivec)
            {
              *convoData[irank]<<std::setprecision(18)<<timeVector[ivec]<<"    "<<std::setprecision(5)<<coefficientVector[ivec]<<std::endl;
            }
	  convoData[irank]->close();
        }


    }

  //cMIC->fillInterpolatorWithData(clusterDataTableTime,clusterDataTableCoefficients);

  convoData.clear();


}


//--------------------------------------------------------
// Name:      ClusterManager::convolveData
// Purpose:   perform convolution integral over clusters
//            of identical rank
// Notes:
//
// Author: Lawrence C Musson
// Date: 29 March 2016
//--------------------------------------------------------

void charon::ClusterManager::convolveData(int irank, std::vector<double> & timeVector, std::vector<double> & coefficientVector)
{

  //The offset across the data table is rankSize
  int offset = rankSize;
  //The number of things to convolve is equal to the file size
  //int cSize = fileSize;

  //The end times of each cluster in the file are recored in cluster time
  //When adjusted for onset time, they all end at identical device times.
  double currentTime = clusterDataTableTime[irank][clusterDataTableTime[irank].size()-1] + cPMaster[irank].getOnsetTime();
  //double wonkyTime = cPMaster[irank].getWonkyTime(); //flag for removal
  //double onsetTime = cPMaster[irank].getOnsetTime(); //flag for removal
  double extendedTime = getExtendedPulseTime();
  double window=extendedTime - pulseOnsetTimes[pulseOnsetTimes.size()-1];
  bool pulseCompleted = false;
  if(currentTime-extendedTime > 0 &&
     currentTime-extendedTime > 0.05*window)  //The pulse has completed
    //if(currentTime > pulseOnsetTimes[pulseOnsetTimes.size()-1])
    pulseCompleted = true;

  std::vector<std::vector<double> > timeVectorExtract, coefficientVectorExtract;

  //extract the raw data.  This just makes it look a little nicer for smoothing.
  timeVectorExtract.resize(fileSize);
  coefficientVectorExtract.resize(fileSize);
  for(int ifile=0 ; ifile<fileSize ; ++ifile)
    {
      timeVectorExtract[ifile] = clusterDataTableTime[irank + ifile*offset];  //This is in cluster time
      coefficientVectorExtract[ifile] = clusterDataTableCoefficients[irank + ifile*offset];
    }

  //Clear the arrays before pushing back
  timeVector.clear();
  coefficientVector.clear();

  //Create a vector of times at which to evaluate the convolution integral
  //run through and count up the number of pulse windows we're working with'
  int pulseWindows = 0;

  if(currentTime > pulseOnsetTimes[pulseOnsetTimes.size()-1])
    {
      pulseWindows = pulseOnsetTimes.size();
      if(pulseCompleted)
        pulseWindows = pulseOnsetTimes.size()+1;
    }
  else
    for(size_t indexes=0 ; indexes < pulseOnsetTimes.size() ; ++indexes)
      {
        if(pulseOnsetTimes[indexes] >= currentTime)
          {
            //This is a little arbitrary, but if time is less than 5% of next window, ignore it.
            if(indexes == 0)pulseWindows=1;
            if(indexes>=1)
              {
                double window=pulseOnsetTimes[indexes] - pulseOnsetTimes[indexes-1];
                if(currentTime-pulseOnsetTimes[indexes] < 0.05*window)
                  pulseWindows = indexes;
                else
                  pulseWindows = indexes+1;
              }
            else
              pulseWindows=1;
            break;
          }
        pulseWindows=1;
      }

  std::vector<std::vector<double> > convoTimes;

  convoTimes.resize(pulseWindows);

 //Do ten increments of each decade and in each pulse window
  //Take this opportunity to identify a partial window--Do it later--maybe unnecessary
  for(int iwindow = 0 ; iwindow<pulseWindows ; ++iwindow)
    {
      //First time will always be 0 + wonkyTime--Working in cluster time
      //double firstTime =  pulseOnsetTimes[iwindow] - onsetTime + wonkyTime;
      //double firstTime =  wonkyTime;
      double firstTime=0;
      if(iwindow == (int)pulseOnsetTimes.size())
        firstTime = extendedTime;
      else
        firstTime = pulseOnsetTimes[iwindow];  //Do this in device time
      //Make a quick check--should do something more sophisticated later
      if(firstTime == 0)
        firstTime = 1e-10;
      //Last time should also be in cluster time---so it's the start
      //of the next pulse minus the start of the current pulse
      double lastTime = 0;
      if(iwindow < (int)pulseOnsetTimes.size()-1)
        lastTime = pulseOnsetTimes[iwindow+1];  //width of the pulse window
        //lastTime = pulseOnsetTimes[iwindow+1] - pulseOnsetTimes[iwindow];
        //lastTime = pulseOnsetTimes[iwindow+1] - onsetTime;
      else if(pulseCompleted && iwindow==pulseWindows-1) //This captures the time past even the extended pulse window
        lastTime = clusterDataTableTime[irank][clusterDataTableTime[irank].size()-1];  //Out to the end of the simulation
        //lastTime = currentTime - pulseOnsetTimes[pulseOnsetTimes.size()-1];
      else
        lastTime = extendedTime; //This captures the extended time past the final pulse
        //lastTime = extendedTime - pulseOnsetTimes[iwindow];
        //lastTime = extendedTime - onsetTime;

      int startExp = floor(log10(firstTime));
      int endExp = floor(log10(lastTime));
      int numberOfDecades = abs(endExp - startExp);
      if(numberOfDecades == 0)numberOfDecades=1;
      int numberOfTimes= 20*numberOfDecades; //equal incrememnts in log time
      double width = log10(lastTime)-log10(firstTime);
      double logTimeIncrement = width/(double)(numberOfTimes-1);

      double currentLogTime = log10(firstTime);
 
      for(int itime=0 ; itime<numberOfTimes ; ++itime)
        {
          if(itime==numberOfTimes-1)
            convoTimes[iwindow].push_back(lastTime);  //Doing this avoids the painful heartache of roundoff in doing logs and powers
          else
            convoTimes[iwindow].push_back(pow(10,currentLogTime));
          currentLogTime += logTimeIncrement;
        }
    }

  std::vector<std::vector<double> > localCoefficient,pulseVector;
  localCoefficient.resize(pulseWindows);
  pulseVector.resize(pulseWindows);

  for(int iwindow=0 ; iwindow<pulseWindows ; ++iwindow)
    {
      for(size_t itime=0 ; itime<convoTimes[iwindow].size() ; ++itime)
        {
          double integratedCoefficient=0;
          int jwindowExtent = iwindow+1;
          if(pulseCompleted && iwindow == pulseWindows-1)
            jwindowExtent = iwindow;
          for(int jwindow=0 ; jwindow<jwindowExtent ; ++jwindow)
            {
              //Figure out where we are in the file
              int ifile=jwindow;
              if(useSingleFileClusters)
                ifile = 0;

              double startTimeOffset = 0;
              if(jwindow<iwindow && iwindow < (int)pulseOnsetTimes.size())
                startTimeOffset=pulseOnsetTimes[jwindow+1];

              if(iwindow == (int)pulseOnsetTimes.size())
                {
                  if(jwindow < (int)pulseOnsetTimes.size()-1)
                    startTimeOffset = pulseOnsetTimes[jwindow+1];
                  else
                    startTimeOffset = extendedTime;
                }

              if(jwindow==iwindow)
                startTimeOffset = pulseOnsetTimes[jwindow];


              //Calculate the time in the pulse where integration begins
              double TcStart = convoTimes[iwindow][itime] - startTimeOffset;
              if(jwindow==iwindow)
                TcStart = 0;

              /*
              if(iwindow == pulseOnsetTimes.size())
                TcStart = convoTimes[iwindow][itime] - extendedTime;
              */

              //double TcEnd = TcStart + convoTimes[iwindow][itime] - startTimeOffset; // when i==j
              //double TcEnd = TcStart + extendedTime-pulseOnsetTimes[pulseOnsetTimes.size()-1];  //This is actually the final one
              double TcEnd = 0;
              if(jwindow==iwindow)
                TcEnd = convoTimes[iwindow][itime] - pulseOnsetTimes[iwindow];  //TcStart should be 0 in this case
              if(jwindow < iwindow)
                TcEnd = convoTimes[iwindow][itime] - pulseOnsetTimes[jwindow];

              if(TcEnd == 0 || TcEnd-TcStart == 0)
                continue;

              std::vector<double> intConvoTimes;
              intConvoTimes = getConvoTimesVector(TcStart,TcEnd);

              //Interpolate the response function across the window at convoTimes
              localCoefficient[jwindow].clear();
              for(size_t iTinterp=0 ; iTinterp<intConvoTimes.size() ; ++iTinterp)
                {
                  //double Tinterp = TcStart + intConvoTimes[iTinterp];
                  double Tinterp = intConvoTimes[iTinterp];
                  double myCoefficient = interpolateCoefficient(irank,ifile,offset,Tinterp,
                                                                timeVectorExtract[ifile],coefficientVectorExtract[ifile]);
                  localCoefficient[jwindow].push_back(myCoefficient);
                }

              //Interpolate the pulse function across the window at convoTimes
              pulseVector[jwindow].clear();
              for(size_t iPinterp=0 ; iPinterp<intConvoTimes.size() ; ++iPinterp)
                {
                  //Note that the interpolation of pulses is done in device time and the spacing has to be reversed from the response
                  double tPStart = pulseOnsetTimes[jwindow];
                  double Tinterp = tPStart + intConvoTimes[intConvoTimes.size()-1] - intConvoTimes[intConvoTimes.size()-1-iPinterp];
                  //double Tinterp = pulseOnsetTimes[jwindow] + intConvoTimes[iPinterp] - intConvoTimes[0];
                  double myPulse = interpolatePulse(jwindow,Tinterp);
                  pulseVector[jwindow].push_back(myPulse);
                }


              double ICThisWindow = integrateWindow(intConvoTimes.size()-1,pulseVector[jwindow],
                                                    localCoefficient[jwindow],intConvoTimes);
              integratedCoefficient += ICThisWindow;

 
            }

          //Calculate the device time and push it into the time vector along with the coefficient we're calculating
          double integralTime = convoTimes[iwindow][itime];
          //double deviceTime = pulseOnsetTimes[iwindow] + integralTime;  //flag for removal
          timeVector.push_back(integralTime);
          coefficientVector.push_back(integratedCoefficient);
        }
    }
}

//--------------------------------------------------------
// Name:      clusterManager::getExtendedPulseTime
// Purpose:
// Notes:
//
// Author: Lawrence C Musson
// Date: 20 July 2016
//--------------------------------------------------------

double charon::ClusterManager::getExtendedPulseTime()
{

  double extendedPulseTime = pulseOnsetTimes[pulseOnsetTimes.size()-1] +
    0.1*(pulseOnsetTimes[pulseOnsetTimes.size()-1] - pulseOnsetTimes[pulseOnsetTimes.size()-2]);

  return extendedPulseTime;

}




//--------------------------------------------------------
// Name:      clusterManager::integrate
// Purpose:
// Notes:
//
// Author: Lawrence C Musson
// Date: 29 March 2016
//--------------------------------------------------------


double charon::ClusterManager::integrateWindow(int extent, std::vector<double> pulseVector,
                                               std::vector<double> responseVector,
                                               std::vector<double> time)
{

  double integral=0;
  for(int convo=0 ; convo<extent ; ++convo)
    {
      int convoReverse = extent - convo;

      integral += 0.5*(responseVector[convo]*pulseVector[convoReverse]
                       + responseVector[convo+1]*pulseVector[convoReverse-1])
        *(time[convo+1] - time[convo]);
    }

  return integral;

}

std::vector<double> charon::ClusterManager::getConvoTimesVector(double Tstart, double Tend)
{

  std::vector<double> time;

  //First time will always be 0 + wonkyTime--Working in cluster time
  //double firstTime =  pulseOnsetTimes[iwindow] - onsetTime + wonkyTime;
  double firstTime =  Tstart;
  //Make a quick check--should do something more sophisticated later
  if(firstTime < 1e-10)
    firstTime = 1e-10;
  //Last time should also be in cluster time---so it's the start
      //of the next pulse minus the start of the current pulse
  double lastTime = Tend;

  int startExp = floor(log10(firstTime));
  int endExp = floor(log10(lastTime));
  int numberOfDecades = abs(endExp - startExp);
  if(numberOfDecades == 0)numberOfDecades=1;
  int numberOfTimes= 20*numberOfDecades; //equal incrememnts in log time
  double width = log10(lastTime)-log10(firstTime);
  double logTimeIncrement = width/(double)(numberOfTimes-1);

  double currentLogTime = log10(firstTime);


  if(Tstart < 1e-10)
    time.push_back(0);
  else
    time.push_back(Tstart);
  for(int itime=1 ; itime<numberOfTimes-1 ; ++itime)
    {
      currentLogTime += logTimeIncrement;
      time.push_back(pow(10,currentLogTime));
    }
  time.push_back(Tend);

  return time;

}

//--------------------------------------------------------
// Name:      clusterManager::integrate
// Purpose:
// Notes:
//
// Author: Lawrence C Musson
// Date: 29 March 2016
//--------------------------------------------------------


double charon::ClusterManager::integrate(int irank, int ifile, int offset, int pulseWindow, double currentTime,
                                         std::vector<std::vector <double> > time, std::vector<std::vector<double> > localCoefficient)
{

  //First and foremost, if there is only onecoefficient, return 0
  if(localCoefficient[pulseWindow].size() == 1) return 0.0;

  //First, construct an array of pulse values that coincide with the reverse time of the response coefficients
  //double onsetTime = cPMaster[irank + ifile*offset].getOnsetTime();  flag for removal

  //double timeMin = pulseOnsetTimes[pulseWindow]-onsetTime;  //shift to cluster time
  double timeMin = time[pulseWindow][0]; //This is already in cluster time.
  double timeMax = currentTime;
  double extendedPulseTime = getExtendedPulseTime()-pulseOnsetTimes[pulseWindow];  //shift to cluster time
  double lateTime;

  if(pulseWindow == (int)pulseOnsetTimes.size())
    lateTime = extendedPulseTime;
  else
    lateTime = pulseOnsetTimes[pulseWindow+1]-pulseOnsetTimes[pulseWindow];

  if(currentTime > lateTime)  //This should probably never happen
    timeMax = lateTime;

  //Calculate the extents of the pulse magnitudes
  double pulseMin = pulseMagnitudes[pulseWindow];
  //If it's past the end of the total pulse, max is zero
  double pulseAbsMax = 0, pulseMax=0;
  //If it's in the pulse, maximum could be at the extent of the pulse windows
  if(pulseWindow < (int)pulseOnsetTimes.size()-1)
    pulseAbsMax = pulseMagnitudes[pulseWindow+1];

  //If it's in the middle of a pulse window, need to interpolate
  if(currentTime <= lateTime)
    {
      double fraction = (currentTime - timeMin)/(lateTime - timeMin);
      pulseMax = pulseMin + fraction*(pulseAbsMax - pulseMin);
    }
  else
    pulseMax = pulseAbsMax;


  std::vector<double> pulseVector;
  pulseVector.resize(localCoefficient[pulseWindow].size());

  //for(size_t i=0 ; i<time[pulseWindow].size() ; ++i)
  for(size_t i=0 ; i<localCoefficient[pulseWindow].size() ; ++i)
    {
      double interpolationTime = timeMax - (time[pulseWindow][i]-timeMin);
      double fraction = (interpolationTime - timeMin)/(timeMax - timeMin);
      pulseVector[localCoefficient[pulseWindow].size() - i - 1] = pulseMin + fraction*(pulseMax - pulseMin);
      //pulseVector[i] = pulseMax + fraction*(pulseMin - pulseMax);
    }


  double convolution=0;

  for(size_t convo=0 ; convo<localCoefficient[pulseWindow].size()-1 ; ++convo)
    {
      int convoReverse = localCoefficient[pulseWindow].size() - convo - 1;

      convolution += 0.5*(localCoefficient[pulseWindow][convo]*pulseVector[convoReverse]
                          + localCoefficient[pulseWindow][convo+1]*pulseVector[convoReverse-1])
        *(time[pulseWindow][convo+1] - time[pulseWindow][convo]);

    }


  return convolution;

  /* OLD STUFF

  double convolution = 0.0;
  size_t totalPoints=0;
  double extendedPulseTime = getExtendedPulseTime();
  bool includePartial = false;

  //run through and count up the number of vectors we're working with'
  if(currentTime > pulseOnsetTimes[pulseOnsetTimes.size()-1])
    totalPoints = pulseOnsetTimes.size();
  else
    for(size_t indexes=0 ; indexes < pulseOnsetTimes.size() ; ++indexes)
      if(pulseOnsetTimes[indexes] > currentTime)
        {
          totalPoints = indexes;
          includePartial = true;
          break;
        }

  //Check for the extended time

  if(totalPoints == 0) return 0.0;

  //local coefficient index--normally this would be identical to the pulse we're integrating.
  //However, if single file clusters are being used, the index will always be zero.
  int lCIndex=0;

  //This section is for responses that have completed their fraction of the total pulse.
  for(size_t convInt=0 ; convInt<totalPoints ; ++convInt)
    {
      if(!useSingleFileClusters)
        lCIndex = convInt;

      //Late time magnitudes and times
      double lateTime=0,lateMagnitude=0;
      if(convInt == pulseOnsetTimes.size()-1)  //When the last pulse is in the extended
        {
          if(currentTime < extendedPulseTime)
            {
              lateTime = currentTime;
              double pulseHi=0.0;
              double pulseLow = pulseMagnitudes[convInt];
              lateMagnitude = pulseMagnitudes[convInt] +
                (pulseHi - pulseLow) *
                (currentTime - pulseOnsetTimes[convInt]) /
                (pulseOnsetTimes[convInt+1] - pulseOnsetTimes[convInt]);

            }
          else
            {
              lateTime = extendedPulseTime;
              lateMagnitude = 0.0;
            }
        }
      else if(convInt < pulseOnsetTimes.size()-1)  //When the last pulse is incomplete, but not in the extended range
        {
          if(currentTime < pulseOnsetTimes[convInt+1])
            {
              lateTime = currentTime;
              double pulseHi  = pulseMagnitudes[convInt+1];
              double pulseLow = pulseMagnitudes[convInt];
              lateMagnitude = pulseMagnitudes[convInt] +
                (pulseHi - pulseLow) *
                (currentTime - pulseOnsetTimes[convInt]) /
                (pulseOnsetTimes[convInt+1] - pulseOnsetTimes[convInt]);

            }
          else
            {
              lateTime = pulseOnsetTimes[convInt+1];
              lateMagnitude = pulseMagnitudes[convInt+1];
            }
        }

      convolution += 0.5*localCoefficient[lCIndex]*
        (pulseMagnitudes[convInt] + lateMagnitude) *
        (lateTime - pulseOnsetTimes[convInt]);

    }

  return convolution;

  */

}

//--------------------------------------------------------
// Name:      clusterManager::interpolatePulse
// Purpose:
// Notes:
//
// Author: Lawrence C Musson
// Date: 1 September 2017
//--------------------------------------------------------

double charon::ClusterManager::interpolatePulse(int jwindow, double Tinterp)
{

  double firstPulse = pulseMagnitudes[jwindow];
  double lastPulse = pulseMagnitudes[0];
  if(lastPulse > firstPulse)
    lastPulse = 0.1*pulseMagnitudes[jwindow];
  if(jwindow < (int)pulseMagnitudes.size()-1)
    lastPulse = pulseMagnitudes[jwindow+1];

  double firstTime = pulseOnsetTimes[jwindow];
  double lastTime=0;
  if(jwindow < (int)pulseMagnitudes.size()-1)
    lastTime = pulseOnsetTimes[jwindow+1];
  else
    lastTime = getExtendedPulseTime();

  double fraction = (Tinterp - firstTime)/(lastTime - firstTime);
  double pulse = firstPulse + fraction*(lastPulse - firstPulse);

  return pulse;

}

//--------------------------------------------------------
// Name:      clusterManager::interpolateCoefficient
// Purpose:
// Notes:
//
// Author: Lawrence C Musson
// Date: 29 March 2016
//--------------------------------------------------------

double charon::ClusterManager::interpolateCoefficient(int irank, int ifile, int offset, double currentTime, std::vector<double> time, std::vector<double> coefficient)
{

  //First check to see if the cluster in question has even started at the time requested.
  //double localOnsetTime = cPMaster[irank + ifile*offset].getOnsetTime();
  //double extendedPulseTime = getExtendedPulseTime();

  if(currentTime<0)  //This should never happen
    return 0;

  double coeff = 0;

  if(currentTime < time[0])
    return coefficient[0];

  for(size_t i=1 ; i<time.size() ; ++i)
    {
      if(currentTime < time[i] &&
         currentTime > time[i-1])
        {
          double fraction = (currentTime - time[i-1])/(time[i] - time[i-1]);
          coeff = coefficient[i-1] + fraction*(coefficient[i] - coefficient[i-1]);
          break;
        }
    }


  /*
  //Hoke in an estimate
  double tlo = time[0];
  double thi = time[time.size()-1];
  double clo = coefficient[0];
  double chi = coefficient[coefficient.size()-1];

  coeff = clo + (currentTime-tlo)/(thi-tlo)*(chi-clo);
  */

  return coeff;


}


//--------------------------------------------------------
// Name:      clusterManagerInterpolatorCommunicator::insertClusterLocationsIntoParameterList()
// Purpose:   Add cluster locations into parameter list.
//            This is part of a kludge to find data at cluster locations -- LCM
// Notes:
//
// Author: Lawrence C Musson
// Date: 24 November 2015
//--------------------------------------------------------


void charon::ClusterManager::insertClusterLocationsIntoParameterList()
{
  using charon::Vector;

  //boost::mpi::environment env;
  boost::mpi::communicator world;

  clusterDataSuitcase cDS;

  //if(me == 0)
  //  cDS.packClusterNames(clusterMasterNameList);

  //This is part of the kludge...should figure out what to do with this
  //Will probably become unnecessary at some point.  If so, use code commented out above.
  if(me==0)
    cDS.packClusterNames(clusterMasterNameList,static_cast<std::vector<double> >(clusterX),
                         static_cast<std::vector<double> >(clusterY),
                         static_cast<std::vector<double> >(clusterZ),
                         static_cast<std::vector<bool> >(clusterFound));

  broadcast(world,cDS,0);

  //if(me != 0)
  //  cDS.unpackClusterNames(clusterMasterNameList);

  //This is part of the kludge...should figure out what to do with this
  //Will probably become unnecessary at some point.  If so, use code commented out above.
  if(me!=0)
    cDS.unpackClusterNames(clusterMasterNameList, clusterX, clusterY, clusterZ, clusterFound);
  
  //Stuff these things into the Teuchos parameter list.

  //This is so damn stupid
  Vector<std::string> localNames;
  for(size_t iname=0 ; iname<clusterMasterNameList.size() ; ++iname)
	localNames.push_back(clusterMasterNameList[iname]);
      
  clusterParameterList->set<Vector<std::string> >("Cluster Names",localNames);
  clusterParameterList->set<Vector<double> >("Cluster X Location",clusterX);
  clusterParameterList->set<Vector<double> >("Cluster Y Location",clusterY);
  clusterParameterList->set<Vector<double> >("Cluster Z Location",clusterZ);
  clusterParameterList->set<Vector<double> >("Cluster Electron Concentration",clusterElectron);
  clusterParameterList->set<Vector<double> >("Cluster Hole Concentration",clusterHole);
  clusterParameterList->set<Vector<double> >("Cluster Acceptor Concentration",clusterAcceptor);
  clusterParameterList->set<Vector<double> >("Cluster Donor Concentration",clusterDonor);
  clusterParameterList->set<Vector<bool> >("Cluster Found",clusterFound);

}

//--------------------------------------------------------
// Name:      clusterManager::extractCarriersAndDopants
// Purpose:   get carrier and doping concentrations from the cluster
//            parameter list and pass around to respective clusters
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


void charon::ClusterManager::extractCarriersAndDopants(double onsetTime)
{

  //This is probably super inefficient, but it should work.
  //the data to be broadcast is fairly small in size and overhead may not be noticable.
  bool setTime=false;

  for(int iproc=0 ; iproc<nprocs ; ++iproc)
    broadcastAProcsData(iproc,onsetTime,setTime);

}

//--------------------------------------------------------
// Name:      clusterManager::brodcastAProcsData
// Purpose:   This method broadcasts the local conditions data for a cluster
//            to all procs.  A check is made on each proc to see if it owns
//            the cluster.  If so, it stores the local conditions in the
//            appropriate cluster parameter object.
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


void charon::ClusterManager::broadcastAProcsData(const int iproc, double onsetTime, bool& setTime)
{


  //boost::mpi::environment env;
  boost::mpi::communicator world;

  clusterDataSuitcase cDS;
  std::vector<std::string> names_;
  std::vector<double> electron_,hole_,acceptor_,donor_;
  std::vector<int> dP;
  std::vector<std::string> rnames_;
  std::vector<double> relectron_,rhole_,racceptor_,rdonor_;
  std::vector<int> rdP;

  if(iproc == me)
    {
      for(size_t icluster=0 ; icluster<interpolator_->clusterNames.size() ; ++icluster)
        {
          if(interpolator_->clusterFound[icluster])
            {
              names_.push_back(interpolator_->clusterNames[icluster]);
              electron_.push_back(interpolator_->clusterElectron[icluster]);
              hole_.push_back(interpolator_->clusterHole[icluster]);
              acceptor_.push_back(interpolator_->clusterAcceptor[icluster]);
              donor_.push_back(interpolator_->clusterDonor[icluster]);
              dP.push_back(me);
            }
        }
      cDS.packClusterBCs(names_,electron_,hole_,acceptor_,donor_,dP);
    }

  //broadcast data
  broadcast(world,cDS,iproc);

  cDS.unpackClusterBCs(rnames_,relectron_,rhole_,racceptor_,rdonor_,rdP);
  //Does it belong here?
  for(size_t icluster=0 ; icluster<rnames_.size() ; ++icluster)
    {
      //loop over local clusters to check if cluster is on this proc
      for(size_t jcluster=0 ; jcluster<cP.size() ; ++jcluster)
        {
          if(rnames_[icluster] == cP[jcluster].getCanonicalName())
            { //store it
              cP[jcluster].setElectronBC(relectron_[icluster]);
              cP[jcluster].setHoleBC(rhole_[icluster]);
              cP[jcluster].setAcceptorConc(racceptor_[icluster]);
              cP[jcluster].setDonorConc(rdonor_[icluster]);
              cP[jcluster].setDeviceProc(rdP[icluster]);
              //cP[jcluster].setOnsetTime(onsetTime); // Don't need to do this here
            }
        }
      //loop over master cluster vector to keep master up to date
      for(size_t jcluster=0 ; jcluster<cPMaster.size() ; ++jcluster)
        {
          if(rnames_[icluster] == cPMaster[jcluster].getCanonicalName())
            { //store it
              cPMaster[jcluster].setElectronBC(relectron_[icluster]);
              cPMaster[jcluster].setHoleBC(rhole_[icluster]);
              cPMaster[jcluster].setAcceptorConc(racceptor_[icluster]);
              cPMaster[jcluster].setDonorConc(rdonor_[icluster]);
              cPMaster[jcluster].setDeviceProc(rdP[icluster]);
              //cP[jcluster].setOnsetTime(onsetTime); // Don't need to do this here
            }
        }
    }
  if( me==0) //update the BC history
    {
      for(size_t icluster=0 ; icluster<rnames_.size() ; ++icluster)
        {
          for(size_t jcluster=0 ; jcluster<cPMaster.size() ; ++jcluster)
            {
              if(rnames_[icluster] == cPMaster[jcluster].getName())
                {
                  if(!setTime)
                    {
                      BCHistoryTime.push_back(onsetTime);
                      setTime = true;
                    }
                  BCHistoryElectrons[jcluster].push_back(relectron_[icluster]);
                  BCHistoryHoles[jcluster].push_back(rhole_[icluster]);
                }
            }
        }
    }


}

//--------------------------------------------------------
// Name:      clusterManager::updateCarrierBCs
// Purpose:   update cluster BCs
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------

void charon::ClusterManager::updateCarrierBCs(double time)
{


  extractCarriersForUpdate(time);

  for(size_t icluster=0 ; icluster < cP.size() ; ++icluster)
    {
      //Kill off negative values
      double eConc=cP[icluster].getElectronBCNumerical();
      double hConc=cP[icluster].getHoleBCNumerical();

      //Don't update if concentrations are negative
      if(eConc < 0 || hConc <0)continue;

      if(eConc<0)
        eConc = 1e3;
      if(hConc<0)
        hConc = 1e3;


      /*  Suppress this for now.
      if(cP[icluster].getUpdateMe())
        {
          xyce[icluster]->modifyClusterBCs(cP[icluster].getName(),eConc,hConc);
          cP[icluster].setUpdateMe(false);
        }
      */
    }
}

//--------------------------------------------------------
// Name:      clusterManager::extractCarriersForUpdate
// Purpose:   get carrier concentrations to update clusters
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


void charon::ClusterManager::extractCarriersForUpdate(double time)
{
  //boost::mpi::environment env;
  boost::mpi::communicator world;

  clusterDataSuitcase cDS;
  std::vector<std::string> names_;
  std::vector<double> electron_,hole_,acceptor_,donor_;
  std::vector<int> dP;
  std::vector<std::string> rnames_;
  std::vector<double> relectron_,rhole_,racceptor_,rdonor_;
  std::vector<int> rdP;

  bool setTime=false;

  names_.clear();
  electron_.clear();
  hole_.clear();
  acceptor_.clear();
  donor_.clear();
  dP.clear();

  //This is probably super inefficient, but it should work.
  //the data to be broadcast is fairly small in size and overhead may not be noticable.

  for(int iproc=0 ; iproc<nprocs ; ++iproc)
    {
      if(iproc == me)
        {
          //clear out transfer vectors so we dont screw up the next pass through iproc
          //Note that this process works on the canonical name.  When conditions are extracted
          //in the charon evaluator, it is stored by canonical name.  This is valid because
          //clusters of the same rank and different file share a physical locatation in the
          //device.
          names_.clear();
          electron_.clear();
          hole_.clear();
          acceptor_.clear();
          donor_.clear();
          dP.clear();
          for(size_t icluster=0 ; icluster<interpolator_->clusterNames.size() ; ++icluster)
            {
              if(interpolator_->clusterFound[icluster])
                {
                  names_.push_back(interpolator_->clusterNames[icluster]);
                  electron_.push_back(interpolator_->clusterElectron[icluster]);
                  hole_.push_back(interpolator_->clusterHole[icluster]);
                  acceptor_.push_back(0);
                  donor_.push_back(0);
                  dP.push_back(0);
                }
            }
          cDS.clearTemps();
          cDS.packClusterBCs(names_,electron_,hole_,acceptor_,donor_,dP);
        }

      //broadcast data
      broadcast(world,cDS,iproc);

      cDS.unpackClusterBCs(rnames_,relectron_,rhole_,racceptor_,rdonor_,rdP);
      cDS.clearTemps();
      //Does it belong here?
      for(size_t icluster=0 ; icluster<rnames_.size() ; ++icluster)
        {
          //loop over local clusters to check if cluster is on this proc
          for(size_t jcluster=0 ; jcluster<cP.size() ; ++jcluster)
            {
              if(rnames_[icluster] == cP[jcluster].getCanonicalName())
                { //store it
                  cP[jcluster].setElectronBC(relectron_[icluster]);
                  cP[jcluster].setHoleBC(rhole_[icluster]);
                }
            }
          //loop over master clusters and set.  Masters are on all procs
          for(size_t jcluster=0 ; jcluster<cPMaster.size() ; ++jcluster)
            {
              if(rnames_[icluster] == cPMaster[jcluster].getCanonicalName())
                { //store it
                  cPMaster[jcluster].setElectronBC(relectron_[icluster]);
                  cPMaster[jcluster].setHoleBC(rhole_[icluster]);
                }
            }
        }
      if( me==0) //update the BC history
        {
          for(size_t icluster=0 ; icluster<rnames_.size() ; ++icluster)
            {
              for(size_t jcluster=0 ; jcluster<cPMaster.size() ; ++jcluster)
                {
                  if(rnames_[icluster] == cPMaster[jcluster].getName())
                    {
                      if(!setTime)
                        {
                          BCHistoryTime.push_back(time);
                          setTime = true;
                        }
                      BCHistoryElectrons[jcluster].push_back(relectron_[icluster]);
                      BCHistoryHoles[jcluster].push_back(rhole_[icluster]);
                    }
                }
            }
        }
    }
}


//--------------------------------------------------------
// Name:      clusterManager::extractBCBreakpoints
// Purpose:   get BC update breakpoints from parameter list
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


std::vector<double> charon::ClusterManager::extractBCBreakpoints()
{

  //The first line determines the number of columns in the file.  This file must be space delimited

  std::string bPLine = clusterParameterList->get<std::string>("Update Cluster BCs Break Points");

  boost::char_separator<char> sep(",");

  boost::tokenizer<boost::char_separator<char> > bPLineTokens(bPLine, sep);

  std::vector<double> localBreakPoints;

  for (const auto& t : bPLineTokens)
    {
      localBreakPoints.push_back(strtod(t.c_str(),NULL));
    }

  //Now extract the pulse break points
  std::vector<int> localPulseBreakPoints;

  bPLine = clusterParameterList->get<std::string>("Pulse Break Points");

  boost::tokenizer<boost::char_separator<char> > pBPLineTokens(bPLine, sep);

  for (const auto& t : pBPLineTokens)
    {
      localPulseBreakPoints.push_back(std::stoi(t.c_str()));
    }


  if(!breakPointsInitialized)
    {
      breakPoints = localBreakPoints;
      pulseBreakPoints = localPulseBreakPoints;
      breakPointsInitialized = true;
    }


  return localBreakPoints;

}

//--------------------------------------------------------
// Name:      clusterManager::addOnsetTImes
// Purpose:   add onset times to the cluster parameter objects
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


void charon::ClusterManager::addOnsetTimes(double oT)
{

  for(size_t icluster=0 ; icluster<cP.size() ; ++icluster)
    {
      cP[icluster].setOnsetTime(oT);
    }

  broadcastOnsetTimesToMasters();

}

//--------------------------------------------------------
// Name:      clusterManager::brodcastOnsetTimesToMasters
// Purpose:   This method broadcasts the onset times to cluster master list
//            on every proc.
// Notes:
//
// Author: Lawrence C Musson
// Date: 15 February 2016
//--------------------------------------------------------


void charon::ClusterManager::broadcastOnsetTimesToMasters()
{
  //boost::mpi::environment env;
  boost::mpi::communicator world;

  std::vector<Xyce::Circuit::clusterParameters> cPLocal;

  for(int iproc=0 ; iproc<nprocs ; ++iproc)
    {

      if(iproc == me)
        cPLocal = cP;

      broadcast(world,cPLocal,iproc);

      for(size_t icluster=0 ; icluster<cPLocal.size() ; ++icluster)
        {
          for(size_t jcluster=0 ; jcluster<cPMaster.size() ; ++jcluster)
            {
              if(cPLocal[icluster].getName() == cPMaster[jcluster].getName())
                {
                  cPMaster[jcluster].setOnsetTime(cPLocal[icluster].getOnsetTime());
                  break;
                }
            }
        }
      cPLocal.clear();
    }

}


//--------------------------------------------------------
// Name:      clusterManagerInterpolatorCommunicator::registerInterpolatorClusters
// Purpose:   register Interpolator clusters with cluster manager
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterManagerInterpolatorCommunicator::
registerInterpolator(Teuchos::RCP<charon::clusterInterpolator> interp)
{

  interpolator_ = interp;

}


//--------------------------------------------------------
// Name:      clusterManagerInterpolatorCommunicator::registerInterpolatorClusters
// Purpose:   register Interpolator clusters with cluster manager
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterManagerInterpolatorCommunicator::
fillInterpolatorWithData( std::vector<std::vector<double> > clusterDataTableTime,
                          std::vector<std::vector<double> > clusterDataTableCoefficients)
{

  //NOTA BENE:  THis method is obsolete and targeted for demolition.

  for(size_t i=0 ; i<clusterDataTableTime.size() ; ++i)
    {
      interpolator_->dPS[i].time = clusterDataTableTime[i];
      interpolator_->dPS[i].sinkLifetime = clusterDataTableCoefficients[i];
    }


}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::packClusterNames
// Purpose:   pack up cluster master name list for broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
packClusterNames(std::vector<std::string> names_)
{

  clusterNameList = names_;

}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::unpackClusterNames
// Purpose:   unpack cluster master name list from broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
unpackClusterNames(std::vector<std::string> &names_)
{

  names_ = clusterNameList;

}


//--------------------------------------------------------
// Name:      clusterDataSuitcase::packClusterNames
// Purpose:   pack up cluster master name list and locations for broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
packClusterNames(std::vector<std::string> names_,  std::vector<double> x_,
                 std::vector<double> y_, std::vector<double> z_, std::vector<bool> cF_)
{

  clusterNameList = names_;
  clusterX = x_;
  clusterY = y_;
  clusterZ = z_;
  clusterFound = cF_;

}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::unpackClusterNames
// Purpose:   unpack cluster master name list and locations from broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
unpackClusterNames(std::vector<std::string> &names_, std::vector<double> &x_,
                   std::vector<double> &y_, std::vector<double> &z_, std::vector<bool> &cF_)
{

  names_ = clusterNameList;
  x_ = clusterX;
  y_ = clusterY;
  z_ = clusterZ;
  cF_ = clusterFound;

}


//--------------------------------------------------------
// Name:      clusterDataSuitcase::clearTemps()
// Purpose:   pack up cluster master name list and BCs
// Notes:
//
// Author: Lawrence C Musson
// Date: 16 February 2016
//--------------------------------------------------------
void charon::clusterDataSuitcase::
clearTemps()
{
  tempClusterNameList.clear();
  tempElectron.clear();
  tempHole.clear();
  tempAcceptor.clear();
  tempDonor.clear();
  tempdeviceProc.clear();

}


//--------------------------------------------------------
// Name:      clusterDataSuitcase::packClusterBCs
// Purpose:   pack up cluster master name list and BCs
// Notes:
//
// Author: Lawrence C Musson
// Date: 16 February 2016
//--------------------------------------------------------

void charon::clusterDataSuitcase::
packClusterBCs(std::vector<std::string> names_,  std::vector<double> electron_,
               std::vector<double> hole_, std::vector<double> acceptor_, std::vector<double> donor_,
               std::vector<int> deviceProc_)
{

  tempClusterNameList = names_;
  tempElectron = electron_;
  tempHole = hole_;
  tempAcceptor = acceptor_;
  tempDonor = donor_;
  tempdeviceProc = deviceProc_;
}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::unpackClusterNames
// Purpose:   unpack cluster master name list and BCs
// Notes:
//
// Author: Lawrence C Musson
// Date: 16 February 2016
//--------------------------------------------------------

void charon::clusterDataSuitcase::
unpackClusterBCs(std::vector<std::string> &names_,  std::vector<double> &electron_,
               std::vector<double> &hole_, std::vector<double> &acceptor_, std::vector<double> &donor_,
               std::vector<int> &deviceProc_)
{

  names_ = tempClusterNameList;
  electron_ = tempElectron;
  hole_ = tempHole;
  acceptor_ = tempAcceptor;
  donor_ = tempDonor;
  deviceProc_ = tempdeviceProc;
}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::packClusterData
// Purpose:   pack up cluster data for broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
packClusterData(std::string names_,std::vector<double> time_, std::vector<double> coefficients_)
{

  clusterName = names_;
  time = time_;
  coefficients = coefficients_;
}

//--------------------------------------------------------
// Name:      clusterDataSuitcase::unpackClusterData
// Purpose:   unpack cluster data from broadcast
// Notes:
//
// Author: Lawrence C Musson
// Date: 14 July 2015
//--------------------------------------------------------

void charon::clusterDataSuitcase::
unpackClusterData(std::string &names_,std::vector<double> &time_, std::vector<double> &coefficients_)
{

  names_ = clusterName;
  time_ = time;
  coefficients_ = coefficients;

}
//-----------------------------------------------------------------------











// ClusterManager
