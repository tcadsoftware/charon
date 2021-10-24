#include <Xyce_config.h>
#include <Charon_config.hpp>
#include <Charon_XyceCluster.hpp>


// ---------- Standard Includes ----------

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <stdexcept>
#include <ctime>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


#include <N_DEV_1D_QASPR.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <streambuf>

namespace Xyce {
namespace Circuit {


//-----------------------------------------------------------------------------
// Function      : Constructor
// Purpose       : Construct
// Special Notes :
// Scope         : public
// Creator       : LawrenceCMusson
// Creation Date :
//-----------------------------------------------------------------------------
CharonXyceInterface::CharonXyceInterface(MPI_Comm comm) : Simulator(comm)
  {
    Xyce::Device::One_D_QASPR::clusterData& cD = Xyce::Device::One_D_QASPR::clusterData::instance();
    cD.setClusterCoupledBool(true);
  }


  //-----------------------------------------------------------------------------
// Function      : getData
// Purpose       : return cluster recombination data
// Special Notes :
// Scope         : public
// Creator       : LawrenceCMusson
// Creation Date :
//-----------------------------------------------------------------------------
void CharonXyceInterface::getData(std::string clusterName, std::vector<double>& time, std::vector<double>& coefficients)
  {
    Xyce::Device::One_D_QASPR::clusterData& cD = Xyce::Device::One_D_QASPR::clusterData::instance();
    cD.getClusterData(clusterName,time,coefficients);
  }


//-----------------------------------------------------------------------------
// Function      : createClusterInputFile
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Lawrence C Musson
// Creation Date : 06/22/2015
//-----------------------------------------------------------------------------
std::string
CharonXyceInterface::createClusterInputFile(std::string clusterTemplate, clusterParameters cP)
{

  std::ifstream cTemp(clusterTemplate.c_str());
  std::string inputFile;

  cTemp.seekg(0, std::ios::end);
  inputFile.reserve(cTemp.tellg());
  cTemp.seekg(0, std::ios::beg);
  inputFile.assign((std::istreambuf_iterator<char>(cTemp)),
             std::istreambuf_iterator<char>());

  size_t f = inputFile.find("NAME");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for NAME.  "<<std::endl;
  inputFile.replace(f, std::string("NAME").length(), cP.getName());

  f = inputFile.find("ELECTRONBC");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for ELECTRONBC.  "<<std::endl;
  inputFile.replace(f, std::string("ELECTRONBC").length(), cP.getElectronBC());

  f = inputFile.find("HOLEBC");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for HOLEBC.  "<<std::endl;
  inputFile.replace(f, std::string("HOLEBC").length(), cP.getHoleBC());

  f = inputFile.find("ACCEPTORCONCENTRATION");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for ACCEPTORCONCENTRATION.  "<<std::endl;
  inputFile.replace(f, std::string("ACCEPTORCONCENTRATION").length(), cP.getAcceptorConc());

  f = inputFile.find("DONORCONCENTRATION");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for DONORCONCENTRATION.  "<<std::endl;
  inputFile.replace(f, std::string("DONORCONCENTRATION").length(), cP.getDonorConc());

  bool saveClusters = cP.getSaveCoefficients();
  if(saveClusters)
    {
      f = inputFile.find("SAVECOEFFICIENTS");
      inputFile.replace(f,std::string("SAVECOEFFICIENTS").length(),"true");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for SAVECOEFFICIENTS.  "<<std::endl;
      f = inputFile.find("CLUSTERCOEFFICIENTFILENAME");
      std::stringstream coefficientFilename;
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for CLUSTERCOEFFICIENTFILENAME.  "<<std::endl;
      coefficientFilename<<cP.getName()<<".coefficients";
      inputFile.replace(f,std::string("CLUSTERCOEFFICIENTFILENAME").length(),coefficientFilename.str());
    }
  else
    {
      f = inputFile.find("SAVECOEFFICIENTS");
  if(f>inputFile.length())std::cout<<inputFile.length()<<"  BOOOOO ERRRORR:  Check your cluster template file for SAVECOEFFICIENTS.  "<<std::endl;
      inputFile.replace(f,std::string("SAVECOEFFICIENTS").length(),"false");
    }

  cTemp.close();

  std::string inputFileDirectory = "clusterInputFiles";

  //Does this directory exit? if not, create it.
  if ( !boost::filesystem::is_directory(inputFileDirectory))
   {
     boost::filesystem::create_directories(inputFileDirectory);
   }

  std::string inputFileName = inputFileDirectory + "/" + cP.getName()+".cir";

  std::ofstream ofile(inputFileName.c_str());

  ofile<<inputFile;

  return inputFileName;

}


//-----------------------------------------------------------------------------
// Function      : modifyClusterBCs
// Purpose       : change the cluster carrier boundary condition
// Special Notes :
// Scope         : public
// Creator       : LawrenceCMusson
// Creation Date : 01/06/2016
//-----------------------------------------------------------------------------
void CharonXyceInterface::modifyClusterBCs(std::string clusterName, double electronConcentration, double holeConcentration)
  {
    Xyce::Device::One_D_QASPR::clusterData& cD = Xyce::Device::One_D_QASPR::clusterData::instance();
    cD.setClusterBoundaryConditions(clusterName, electronConcentration, holeConcentration);
    cD.setBCsModified(clusterName, true);
  }

}
}
