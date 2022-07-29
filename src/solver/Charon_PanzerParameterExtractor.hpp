
#ifndef CHARON_PANZERPARAMETEREXTRACTOR_HPP
#define CHARON_PANZERPARAMETEREXTRACTOR_HPP

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Comm.hpp"
#include "Thyra_ModelEvaluator.hpp"


#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Panzer_GlobalData.hpp"
#include <vector>
#include <string>


/********************************************************************************************************
  Create a class that can grab and return patameters from the panzer parameter library.
*********************************************************************************************************/

namespace charon {

  class panzerParameterExtractor {
    
  private:

    Teuchos::RCP<panzer::ParamLib> parameterLibrary_;
    int me, nprocs;
    std::vector<std::pair<std::string,double> > tempMap;
    std::pair<std::string,double> endPair;
    int numParameters,strLen;


  public:

    panzerParameterExtractor(Teuchos::RCP<panzer::ParamLib> const& parameterLibrary)
      : parameterLibrary_(parameterLibrary),
	numParameters(0)
    {endPair.first = "AllDone"; endPair.second = 0.0;}
    

    /***********************************************************************************
    / get the names of the contacts out of the parameter library
    ************************************************************************************/
     std::map<std::string,double>  get_VoltageParameters()
    {
      MPI_Comm_rank(MPI_COMM_WORLD,&me);
      MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

						
      std::map<std::string,double> workingMap;
      panzer::ParamLib::iterator pLit = parameterLibrary_->begin();
      panzer::ParamLib::iterator pLend = parameterLibrary_->end();
      std::string voltageTag = "_Voltage";
      tempMap.clear();

      while (pLit != pLend)
	{
	  const std::string name = (*pLit).first;
	  std::pair<std::string,double> temptempMap;
	  size_t pos = name.find(voltageTag);
	  int length = voltageTag.length();
	  if (pos != std::string::npos and
	      pos + length == name.length())
	    {
	      temptempMap.first = name;
	      temptempMap.second = parameterLibrary_->getValue<panzer::Traits::Residual>(name);
	      tempMap.push_back(temptempMap);
	    }
	  ++pLit;
	}

      //If on proc 0, put the tempMap into the workingMap
      if (me == 0)
	for (size_t imap=0 ; imap<tempMap.size() ; ++imap)
	  if (tempMap[imap].first != "AllDone")
	    workingMap[tempMap[imap].first] = tempMap[imap].second;

      std::vector<std::string> bufSendString;
      std::string bufRecString;
      MPI_Status status;
      std::string bufString;
      std::string recString;
      std::vector<int> expected(nprocs,0);

      //First order of business is to tell proc 0 how many items to expect.  Minimum is one
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  if (me==0)
	    {
	      int recNum;
	      MPI_Recv(&recNum,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);
	      expected[iproc] = recNum;
	    }
	  else if (me == iproc)
	    {
	      int sendNum = tempMap.size();
	      expected[me] = sendNum;
	      MPI_Send(&sendNum,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
	    }
	}


      //pass all the necessary stuff
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  for (int ipair=0 ; ipair<expected[iproc] ; ++ipair)
	    {
	      //Set up the receive
	      if (me == 0)
	      {
		strLen = 0;
		MPI_Recv(&strLen,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);

		std::vector<char> receivedString(strLen + 1,0);
		MPI_Recv(&receivedString[0], strLen, MPI_CHAR, iproc, iproc, MPI_COMM_WORLD, &status);
		recString = &receivedString[0];

		double gotDouble = -1;
		MPI_Recv(&gotDouble,1,MPI_DOUBLE,iproc,  iproc,MPI_COMM_WORLD,&status);

		workingMap[&receivedString[0]] = gotDouble;

	      }
	    else if (iproc == me)
	      //set up the send
	      {
		
		strLen = tempMap[ipair].first.size();
		MPI_Send(&strLen,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
		
		bufString = tempMap[ipair].first;
		MPI_Send(bufString.c_str(), bufString.size(), MPI_CHAR, 0, iproc, MPI_COMM_WORLD);

		double sendDouble = tempMap[ipair].second;
		MPI_Send(&sendDouble,1,MPI_DOUBLE,0,iproc,MPI_COMM_WORLD);
		
	      }
	    }
	}

      //This workingMap must now be broadcast to every proc so that exodus is happy.  What a pain.
      int expectThisMany=workingMap.size();
      MPI_Bcast(&expectThisMany,1,MPI_INT,0,MPI_COMM_WORLD);

      //pass all the necessary stuff.  this feels like the dumbest way to do this, but it's this or screw around
      //with data serialization for a tiny amount of data.

      tempMap.clear();
      if (me == 0)
	{
	  std::map<std::string,double>::iterator wmit = workingMap.begin();
	  while( wmit != workingMap.end())
	    {
	      tempMap.push_back(*wmit);
	      ++wmit;
	    }
	}

      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  for (int ipair=0 ; ipair<expectThisMany; ++ipair)
	    {
	      //Set up the receive
	      if (me == iproc)
	      {
		strLen = 0;
		MPI_Recv(&strLen,1,MPI_INT,0,  iproc,MPI_COMM_WORLD,&status);

		std::vector<char> receivedString(strLen + 1,0);
		MPI_Recv(&receivedString[0], strLen, MPI_CHAR, 0, iproc, MPI_COMM_WORLD, &status);
		recString = &receivedString[0];

		double gotDouble = -1;
		MPI_Recv(&gotDouble,1,MPI_DOUBLE, 0,  iproc,MPI_COMM_WORLD,&status);

		workingMap[&receivedString[0]] = gotDouble;

	      }
	    else if (me == 0)
	      //set up the send
	      {
		
		strLen = tempMap[ipair].first.size();
		MPI_Send(&strLen,1,MPI_INT,iproc,iproc,MPI_COMM_WORLD);
		
		bufString = tempMap[ipair].first;
		MPI_Send(bufString.c_str(), bufString.size(), MPI_CHAR, iproc, iproc, MPI_COMM_WORLD);

		double sendDouble = tempMap[ipair].second;
		MPI_Send(&sendDouble,1,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD);
		
	      }
	    }
	}

      return workingMap;
    }


    /***********************************************************************************
    / If the continuation parameter name is stored in the parameter library, retrieve it now.
    ************************************************************************************/
    std::string  getContinuationParameterName()
    {

      MPI_Comm_rank(MPI_COMM_WORLD,&me);
      MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

      panzer::ParamLib::iterator pLit = parameterLibrary_->begin();
      panzer::ParamLib::iterator pLend = parameterLibrary_->end();

      std::string noName =  "NoName";
      std::string returnName = noName;
      std::string locaTag = "IS LOCA";
      while (pLit != pLend)
	{
	  const std::string name = (*pLit).first;
	  size_t pos = name.find(locaTag);
	  if (pos != std::string::npos)
	    {
	      //Pull off the loca tag
	      returnName = name;
	      returnName.erase(pos, locaTag.length());
	    }
	  ++pLit;
	}

      std::vector<std::string> bufSendString;
      std::string bufRecString;
      MPI_Status status;
      std::string bufString;
      std::vector<int> expected(nprocs,0);

      //First order of business is to tell proc 0 how many items to expect.  Minimum is one
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  if (me==0)
	    {
	      int recNum;
	      MPI_Recv(&recNum,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);
	      expected[iproc] = recNum;
	    }
	  else if (me == iproc)
	    {
	      int sendNum = 0;
	      if (returnName != noName)
		sendNum = 1;
	      expected[me] = sendNum;
	      MPI_Send(&sendNum,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
	    }
	}


      //pass all the necessary stuff
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  for (int ipair=0 ; ipair<expected[iproc] ; ++ipair)
	    {
	      //Set up the receive
	      if (me == 0)
	      {
		strLen = 0;
		MPI_Recv(&strLen,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);

		std::vector<char> receivedString(strLen + 1,0);
		MPI_Recv(&receivedString[0], strLen, MPI_CHAR, iproc, iproc, MPI_COMM_WORLD, &status);
		returnName = &receivedString[0];
	      }
	    else if (iproc == me)
	      //set up the send
	      {
		
		strLen = tempMap[ipair].first.size();
		MPI_Send(&strLen,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
		
		bufString = returnName;
		MPI_Send(bufString.c_str(), bufString.size(), MPI_CHAR, 0, iproc, MPI_COMM_WORLD);
	      }
	    }
	}



      return returnName;
    }

    /***********************************************************************************
    / Pull out the values of the contact voltages and send the back to the writers.
    ************************************************************************************/
   double getVoltage(std::string name)
    {
      return parameterLibrary_->getValue<panzer::Traits::Residual>(name);
    }


	/***********************************************************************************
    / get the names of the norms out of the parameter library
    ************************************************************************************/
     std::map<std::string,double>  get_NormParameters()
    {
      MPI_Comm_rank(MPI_COMM_WORLD,&me);
      MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

						
      std::map<std::string,double> workingMap;
      panzer::ParamLib::iterator pLit = parameterLibrary_->begin();
      panzer::ParamLib::iterator pLend = parameterLibrary_->end();
      //std::string normTag = "_L2_Norm_Error";
	  std::string normTag = "_Norm";
      tempMap.clear();

      while (pLit != pLend)
	{
	  const std::string name = (*pLit).first;
	  std::pair<std::string,double> temptempMap;
	  size_t pos = name.find(normTag);
	  int length = normTag.length();
	  if (pos != std::string::npos and
	      pos + length == name.length())
	    {
	      temptempMap.first = name;
	      temptempMap.second = parameterLibrary_->getValue<panzer::Traits::Residual>(name);
	      tempMap.push_back(temptempMap);
	    }
	  ++pLit;
	}

      //If on proc 0, put the tempMap into the workingMap
      if (me == 0)
	for (size_t imap=0 ; imap<tempMap.size() ; ++imap)
	  if (tempMap[imap].first != "AllDone")
	    workingMap[tempMap[imap].first] = tempMap[imap].second;

      std::vector<std::string> bufSendString;
      std::string bufRecString;
      MPI_Status status;
      std::string bufString;
      std::string recString;
      std::vector<int> expected(nprocs,0);

      //First order of business is to tell proc 0 how many items to expect.  Minimum is one
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  if (me==0)
	    {
	      int recNum;
	      MPI_Recv(&recNum,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);
	      expected[iproc] = recNum;
	    }
	  else if (me == iproc)
	    {
	      int sendNum = tempMap.size();
	      expected[me] = sendNum;
	      MPI_Send(&sendNum,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
	    }
	}


      //pass all the necessary stuff
      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  for (int ipair=0 ; ipair<expected[iproc] ; ++ipair)
	    {
	      //Set up the receive
	      if (me == 0)
	      {
		strLen = 0;
		MPI_Recv(&strLen,1,MPI_INT,iproc,  iproc,MPI_COMM_WORLD,&status);

		std::vector<char> receivedString(strLen + 1,0);
		MPI_Recv(&receivedString[0], strLen, MPI_CHAR, iproc, iproc, MPI_COMM_WORLD, &status);
		recString = &receivedString[0];

		double gotDouble = -1;
		MPI_Recv(&gotDouble,1,MPI_DOUBLE,iproc,  iproc,MPI_COMM_WORLD,&status);

		workingMap[&receivedString[0]] = gotDouble;

	      }
	    else if (iproc == me)
	      //set up the send
	      {
		
		strLen = tempMap[ipair].first.size();
		MPI_Send(&strLen,1,MPI_INT,0,iproc,MPI_COMM_WORLD);
		
		bufString = tempMap[ipair].first;
		MPI_Send(bufString.c_str(), bufString.size(), MPI_CHAR, 0, iproc, MPI_COMM_WORLD);

		double sendDouble = tempMap[ipair].second;
		MPI_Send(&sendDouble,1,MPI_DOUBLE,0,iproc,MPI_COMM_WORLD);
		
	      }
	    }
	}

      //This workingMap must now be broadcast to every proc so that exodus is happy.  What a pain.
      int expectThisMany=workingMap.size();
      MPI_Bcast(&expectThisMany,1,MPI_INT,0,MPI_COMM_WORLD);

      //pass all the necessary stuff.  this feels like the dumbest way to do this, but it's this or screw around
      //with data serialization for a tiny amount of data.

      tempMap.clear();
      if (me == 0)
	{
	  std::map<std::string,double>::iterator wmit = workingMap.begin();
	  while( wmit != workingMap.end())
	    {
	      tempMap.push_back(*wmit);
	      ++wmit;
	    }
	}

      for (int iproc = 1 ; iproc<nprocs ; ++iproc)
	{
	  for (int ipair=0 ; ipair<expectThisMany; ++ipair)
	    {
	      //Set up the receive
	      if (me == iproc)
	      {
		strLen = 0;
		MPI_Recv(&strLen,1,MPI_INT,0,  iproc,MPI_COMM_WORLD,&status);

		std::vector<char> receivedString(strLen + 1,0);
		MPI_Recv(&receivedString[0], strLen, MPI_CHAR, 0, iproc, MPI_COMM_WORLD, &status);
		recString = &receivedString[0];

		double gotDouble = -1;
		MPI_Recv(&gotDouble,1,MPI_DOUBLE, 0,  iproc,MPI_COMM_WORLD,&status);

		workingMap[&receivedString[0]] = gotDouble;

	      }
	    else if (me == 0)
	      //set up the send
	      {
		
		strLen = tempMap[ipair].first.size();
		MPI_Send(&strLen,1,MPI_INT,iproc,iproc,MPI_COMM_WORLD);
		
		bufString = tempMap[ipair].first;
		MPI_Send(bufString.c_str(), bufString.size(), MPI_CHAR, iproc, iproc, MPI_COMM_WORLD);

		double sendDouble = tempMap[ipair].second;
		MPI_Send(&sendDouble,1,MPI_DOUBLE,iproc,iproc,MPI_COMM_WORLD);
		
	      }
	    }
	}

      return workingMap;
    }

  };  //end class panzerParameterExtractor

}  //end namespace charon


#endif
