// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Charon_CoupledModelEvaluator_impl.hpp
//
///////////////////////////////////////////////////////////////////////////////

#ifndef __Charon_CoupledModelEvaluator_impl_hpp__
#define __Charon_CoupledModelEvaluator_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

#include "Charon_config.hpp"
#include "Charon_PanzerParameterExtractor.hpp"

#include "Thyra_DefaultFiniteDifferenceModelEvaluator_decl.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_describeLinearOp_decl.hpp"
#include "Thyra_VectorStdOps_def.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// include Tpetra stuff to modify Thyra vectors
#include "Teuchos_dyn_cast.hpp"
#include "Epetra_CrsMatrix.h"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Panzer stuff for DOF indexing
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

// Panzer stuff for hooking in the LOCA parameters
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

#include "NOX_TpetraTypedefs.hpp" // For underlying types for LOCA Bordering
#include "Thyra_DetachedVectorView.hpp"
//#include "Thyra_TpetraThyraWrappers.hpp" // Convert Thyra to Tpetra
//#include "Thyra_DefaultSpmdVector_decl.hpp" // Convert Thyra to Tpetra

#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
#include <N_CIR_SecondLevelSimulator.h>
#include <N_TIA_TwoLevelError.h>
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER

namespace charon {

// Constructors/initializers/accessors/utilities
template<class Scalar>
CoupledModelEvaluator<Scalar>::CoupledModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
    MPI_Comm                                           comm,
    const Teuchos::RCP<Teuchos::ParameterList>&        parameters,
    charon::CurrentConstraintList&                     constraints,
    const Teuchos::RCP<panzer_stk::STK_Interface>&     mesh)
  :
  Thyra::ModelEvaluatorDelegatorBase<Scalar>(physics),
  physics_(physics),
  comm_(Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(
    Teuchos::opaqueWrapper(comm)))),
  parameters_(parameters),
  dimension_(mesh->getDimension()),
  mesh_(mesh),
  constraints_(constraints)
{
  // Ensure we were actually given a physics ModelEvaluator.
  TEUCHOS_ASSERT(not physics_.is_null())

  // Ensure we have a global indexer
  TEUCHOS_ASSERT(parameters_->isParameter("Unique Global Indexer"))

  globalIndexer_ = parameters->get<Teuchos::RCP<panzer::GlobalIndexer > >("Unique Global Indexer");
  Teuchos::RCP<panzer::GlobalData> globalData =  parameters->get<Teuchos::RCP<panzer::GlobalData > >("Global Data");
  paramLib_ = globalData->pl;

  xyceCouplingPL_ = Teuchos::rcp(new Teuchos::ParameterList("Xyce Coupling Params"));
  xyceCouplingPL_->set<int>("Coupling Step Number", 0); // for MMvC
  xyceCouplingPL_->set<double>("Initial Voltage", 1.0); // for MMvV
  xyceCouplingPL_->set<int>("Initialize Coupling Step", parameters_->get<int>("Initial Xyce Coupling on Step Number",0));
  xyceCouplingPL_->set<std::string>("Netlist", parameters_->get<std::string>("Netlist","No_netlist_for_MMvV"));
  xyceCouplingPL_->set<double>("Xyce DC Voltage", 0.0);
  mixedModeLOCA = paramLib_->isParameterForType<panzer::Traits::Residual>("Xyce DC Voltage");

#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
  if(xyceCouplingPL_->get<std::string>("Netlist") != "No_netlist_for_MMvV"){
    simulator_ = new Xyce::Circuit::SecondLevelSimulator();
  
    std::cout << "Charon-Xyce coupling will use this netlist: " \
              << xyceCouplingPL_->get<std::string>("Netlist") << std::endl;
  
    int argc = 2;
    char *argv[3];
    argv[0] = strdup("Xyce");
    if (printDebug) std::cout << "initializing Xyce... ";
    argv[1] = strdup((xyceCouplingPL_->get<std::string>("Netlist")).c_str());
    argv[2] = 0;
    simulator_->initialize(argc, argv);
    if (printDebug) std::cout << "done!" << std::endl;
    delete[] argv[0];
    delete[] argv[1];
  
    simulator_->startupSolvers();
  } // for MMvC
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER

}

// Public function overridden from Teuchos::Describable
template<class Scalar>
std::string CoupledModelEvaluator<Scalar>::description() const
{
  return this->getUnderlyingModel()->description();
}

// Private function overridden from ModelEvaulatorDefaultBase
template<class Scalar>
void CoupledModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  using Thyra::get_ele;

  typename ::Thyra::ModelEvaluatorBase::InArgs<Scalar> physicsInArgs(inArgs);
  typename ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> physicsOutArgs(outArgs);

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > thyraModel = 
    Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar> >(this->getUnderlyingModel());

#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
  // begin MMvC block
  if(xyceCouplingPL_->get<std::string>("Netlist") != "No_netlist_for_MMvV"){
    Xyce::TimeIntg::TwoLevelError tlError;  // data structure related to global error control
    bool initJctFlag(true);
  
    // get the contact_coupling_map (map between g parameter index and contact name)
    Teuchos::RCP<charon::panzerParameterExtractor> pPE = 
      Teuchos::rcp(new charon::panzerParameterExtractor(paramLib_));
    std::map<std::string,double> paramMap =  pPE->get_VoltageParameters();
    Teuchos::RCP<Teuchos::ParameterList> contact_coupling_map =
      parameters_->get<Teuchos::RCP<Teuchos::ParameterList> >("Contact Coupling Map");
    if(printDebug){
      std::cout << "This is the contact coupling map: " << std::endl 
  	      << *contact_coupling_map << std::endl;
    }
    // set the voltageInputMap by iterating through the constraints_, which requires figuring out the index of the responses g
    // in this loop, also set up the Xyce coupling sensitivities (dI/dV)
    std::map<std::string,double> voltageInputMap;
    for (int i=0; i < constraints_.size(); ++i){
      int pIdx = constraints_[i]->parameterIndex(); // response index = pidx
      std::string ssid = constraints_[i]->sidesetId();
      std::string xyceNodeName = contact_coupling_map->get<std::string>(ssid);
      voltageInputMap[xyceNodeName] = paramMap[ssid+"_Voltage"]; // using the panzerParameterExtractor is MPI safe
      // note: can't do paramLib_->getValue<panzer::Traits::Residual>(ssid+"_Voltage") here since it won't be available on all procs
      if(printDebug){
        std::cout << "The contact voltages of the Charon device are:" << std::endl
                  << "g " << pIdx << " / " << ssid << " which uses node named " 
                  << contact_coupling_map->get<std::string>(ssid) << std::endl;
      }
      std::string couplingSensitivity("Xyce Coupling Sensitivity"+std::to_string(i));
      panzer::registerScalarParameter(couplingSensitivity, *paramLib_, 0.0);
    }
    // get the numContacts from contact_coupling_map length
    int numContacts = constraints_.size() + (mixedModeLOCA ? 1 : 0); // NOTE: this assumes all constraints are xyce contacts!
    // initialize outputVector to store Xyce contact currents
    std::vector<double> outputVector(numContacts, 0.0);
  
    // print the contact currents and voltages
    if(printDebug){
      std::cout << "Xyce coupling voltages are: " << std::endl;
      for(int i = 0 ; i < constraints_.size() ; i++){
        int pIdx = constraints_[i]->parameterIndex();
        std::string ssid = constraints_[i]->sidesetId();
        std::string xyceNodeName = contact_coupling_map->get<std::string>(ssid);
        std::cout << "constraint number " << pIdx << ": " <<  ssid << "/" << xyceNodeName << ": " << voltageInputMap[xyceNodeName] << std::endl;
      }
    }
  
    // make a square Jacobian stamp for the coupling terms
    std::vector< std::vector<double> > innerJacobianStamp(numContacts, std::vector<double>(numContacts, 0.0));
  
    // run Xyce with the new Charon device's contact voltages
    if(mixedModeLOCA){
      std::string  xyceSweptNodeName = "vswept";
      voltageInputMap[xyceSweptNodeName] = xyceCouplingPL_->get<double>("Xyce DC Voltage",0.0);
      if(printDebug || printVerbose) std::cout << xyceSweptNodeName << " = " << voltageInputMap[xyceSweptNodeName] << std::endl;
    }
    bool xyceStepTaken = simulator_->simulateStep(initJctFlag, voltageInputMap, outputVector, innerJacobianStamp, tlError);
    if(printDebug){
      std::cout << "Did Xyce take a step? " << xyceStepTaken << std::endl;
      std::cout << "Here are the results: ";
      for (auto outputParam : outputVector)
        std::cout << outputParam << " ";
      std::cout << std::endl;
    }
  
    // the xyce inner Jacobian stamp data is exaactly what we need to update dg/dp entries
    // the diagonal entries are what we need to push to the vector of Xyce coupling sensitivities
    if(xyceCouplingPL_->get<int>("Coupling Step Number") >= xyceCouplingPL_->get<int>("Initialize Coupling Step")){
    //if(xyceCouplingPL_->get<int>("Coupling Step Number") % xyceCouplingPL_->get<int>("Initialize Coupling Step") < 1){
      if(printDebug){
        std::cout << "Mixed mode Step Number: " << xyceCouplingPL_->get<int>("Coupling Step Number") << std::endl; 
        std::cout << "Updating the current constraint derivative information." << std::endl;
      }
      for(int i = 0; i < constraints_.size() ; i++)
        paramLib_->setValue<panzer::Traits::Residual>("Xyce Coupling Sensitivity"+std::to_string(i), innerJacobianStamp[i][i]);
    }
    else
      std::cout << "Don't run Xyce until step " << xyceCouplingPL_->get<int>("Initialize Coupling Step") << std::endl; 
  
    if(printDebug){
      std::cout << "Here is the inner Jacobian stamp " << std::endl;
      for(int i = 0; i < constraints_.size() ; i++)
        std::cout << "dI_i/dV_j[" << i << "," << i << "] = " << innerJacobianStamp[i][i] << std::endl;
    }
  
    // use the Xyce coupled current (outputVector) to
    // update the current constraint values used by the Householder solver
    Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_ = physics_->getNominalValues();
    for (int i=0; i < constraints_.size(); ++i){    
      int pIdx = constraints_[i]->parameterIndex();
      double scaleFactor = constraints_[i]->contactLength() / constraints_[i]->contactArea() * 1.0e4;
      if(printDebug)
        std::cout << "The scale factor on contact " << i << " is: " << scaleFactor << std::endl;
      if(xyceCouplingPL_->get<int>("Coupling Step Number") >= xyceCouplingPL_->get<int>("Initialize Coupling Step")){
      //if(xyceCouplingPL_->get<int>("Coupling Step Number") % xyceCouplingPL_->get<int>("Initialize Coupling Step") < 1){
        constraints_[i]->currentValueUpdate(-1.0*outputVector[i]*scaleFactor); 
        // check: the -1 factor is explained by: what *goes out of Xyce* will *go in to Charon*
      }
      else
        std::cout << "Don't run Xyce on the initial coupling step." << std::endl; 
      if(printDebug){
        std::cout << "This is the constraint index: " << i << std::endl
                  << "This is the constraint response index: " << constraints_[i]->responseIndex() << std::endl
                  << "This is the constraint parameter index: " << pIdx << " " << std::endl
                  << "This is the constraint parameter after update:" <<  constraints_[i]->currentValue() << std::endl;
      }
    }
  } // end MMvC block

  if(xyceCouplingPL_->get<std::string>("Netlist") == "No_netlist_for_MMvV"){
    // Allocate vectors for each of the responses, using a new non const outArgs
    Thyra::ModelEvaluatorBase::OutArgs<Scalar> nonConstOutArgs = outArgs;
    for(int i=0; i < outArgs.Ng(); ++i)
      nonConstOutArgs.set_g(i,Thyra::createMember(*(thyraModel->get_g_space(i))));
  
    std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>> responseValues(nonConstOutArgs.Ng());
    //std::vector<Teuchos::RCP<Thyra::VectorBase<Scalar>>> responseValues(outArgs.Ng());
    /* output for debugging
    // Extract the response information that is required on output.
    std::cout << "After evalModel, there are " << responseValues.size() << " response values for nonConst outArgs: " << std::endl;
    for (int ri(0); ri < nonConstOutArgs.Ng(); ++ri){
      responseValues[ri] = nonConstOutArgs.get_g(ri);
      if(nonnull(responseValues[ri])){
        //std::cout << *responseValues[ri] << std::endl;
        //std::cout << get_ele<Scalar>(*responseValues[ri], 0);
        for (Thyra::Ordinal k(0); k < responseValues[ri]->space()->dim(); ++k)
          std::cout << " " << get_ele<Scalar>(*responseValues[ri], k);
        std::cout << std::endl;
      }
    }
    */
  
    //////////////////////////////////////////////////////////////////////
    //
    //  crudely couple a 1k ohm resistor to a Charon diode
    //
    //       o ---- resistor ---- o --- anode| diode |cathode --- o
    //
    //      3V                 Vcouple                            1V
    //
    //  Solve for Vcouple iteratively
    //  Initialize: Vold = 1.5V
    //  Iterate:
    //    Charon solve: apply Vold/0V at anode/cathode -> get I from Charon
    //    After Newton step, calculate Vnew = Vleft - Iold*(1k ohms)
    //    update Vold to take the value of Vnew
    //  When converged, Vcouple = Vold
    //
    /////////////////////////////////////////////////////////////////////
  
    std::string paramName("Xyce Coupled Voltage");
    double paramValue = paramLib_->getValue<panzer::Traits::Residual>(paramName);
    std::cout << "The current coupling value was: " << paramValue << std::endl;
    // why does the following need to be set to the initial value? if we dont update it here, it won't converge
    paramLib_->setValue<panzer::Traits::Residual>(paramName, xyceCouplingPL_->get<double>("Initial Voltage"));
    thyraModel->evalModel(inArgs, outArgs);
    // update that value
    double Vleft = 3.0;
    double R = 1000.0;
    double Iold = 0;
    //for (int ri(0); ri < nonConstOutArgs.Ng(); ++ri){
    for (int ri(0); ri < outArgs.Ng(); ++ri){
      //responseValues[ri] = nonConstOutArgs.get_g(ri);
      responseValues[ri] = outArgs.get_g(ri);
      if(Teuchos::nonnull(responseValues[ri])){
        for (Thyra::Ordinal k(0); k < responseValues[ri]->space()->dim(); ++k)
  	if(get_ele<Scalar>(*responseValues[ri], k) > 0.0)
            Iold = get_ele<Scalar>(*responseValues[ri], k);
      }
    }
    double Vnew = Vleft - Iold*R;
    paramLib_->setValue<panzer::Traits::Residual>(paramName, Vnew);
  
  } // MMvV block
  
  // print the parameter library for debugging
  if(printDebug)
    paramLib_->print(std::cout);

#else
  std::cout << "Error: Can't couple Charon to Xyce without either "
            << "ENABLE_XYCE_CLUSTER or ENABLE_MIXED_MODE bool set to ON "
	    << "in opts file for build configuration." << std::endl;
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER

  // finally, call this evalModel layer with the modified p in the physicsInArgs
  thyraModel->evalModel(physicsInArgs, physicsOutArgs);
  xyceCouplingPL_->set<int>("Coupling Step Number",  xyceCouplingPL_->get<int>("Coupling Step Number") + 1); 
  if(printDebug){
    std::cout << "Mixed Mode step number: " << xyceCouplingPL_->get<int>("Coupling Step Number") << std::endl;
    std::cout << "Xyce initial step number" <<  xyceCouplingPL_->get<int>("Initialize Coupling Step") << std::endl;
  }

  if(xyceCouplingPL_->get<std::string>("Netlist") != "No_netlist_for_MMvV"){
    // if performing a LOCA sweep, print the modelevaluator parameters for debugging
    std::string locaParamName("Xyce DC Voltage");
    for (int i=0; i < physics_->Np() / 2 ; ++i) { // dont need their tangents
      auto p_names = physics_->get_p_names(i);
      for (int j=0; j < p_names->size(); ++j){
        Teuchos::ArrayRCP<const double> p_data;
        dynamic_cast<const Thyra::SpmdVectorBase<double> &>(*(physicsInArgs.get_p(i))).getLocalData(Teuchos::ptrFromRef(p_data));
        if(printDebug || printVerbose)
          std::cout << "ME Parameter(dim=" << i << ", index=" << j << ") called " << (*p_names)[j]  << " = " << p_data[0] << std::endl;
        if((*p_names)[j] == locaParamName){
  	std::cout << "Let's use this value for Xyce! " << p_data[0] << std::endl;
          xyceCouplingPL_->set<double>("Xyce DC Voltage",  p_data[0]);
        }
      }
    }
  } // end MMvC LOCA debug print block
  
}

} // namespace charon

#endif // __Charon_CoupledModelEvaluator_impl_hpp__
