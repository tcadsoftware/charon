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

#include "Thyra_DefaultFiniteDifferenceModelEvaluator_decl.hpp"
#include "Thyra_ScaledLinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_describeLinearOp_decl.hpp"

// include Tpetra stuff to modify Thyra vectors
#include "Teuchos_dyn_cast.hpp"
#include "Epetra_CrsMatrix.h"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Panzer stuff for DOF indexing
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

namespace charon {


// Constructors/initializers/accessors/utilities


template<class Scalar>
CoupledModelEvaluator<Scalar>::CoupledModelEvaluator(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
    MPI_Comm                                           comm,
    const Teuchos::RCP<Teuchos::ParameterList>&        parameters,
    const Teuchos::RCP<panzer_stk::STK_Interface>&     mesh)
  :
  Thyra::ModelEvaluatorDelegatorBase<Scalar>(physics),
  physics_(physics),
  //comm_(Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(
  //                                                        Teuchos::opaqueWrapper(comm)))),
  parameters_(parameters),
  dimension_(mesh->getDimension()),
  mesh_(mesh)
{
  // Ensure we were actually given a physics ModelEvaluator.
  TEUCHOS_ASSERT(not physics_.is_null())

  // Ensure we have a global indexer
  TEUCHOS_ASSERT(parameters_->isParameter("Unique Global Indexer"))

  globalIndexer_ = parameters->get<Teuchos::RCP<panzer::GlobalIndexer > >("Unique Global Indexer");

  // from Panzer EquationSet Default Impl...
  //Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer;
  //if(lof!=Teuchos::null)
  //  globalIndexer = lof->getRangeGlobalIndexer();
}

// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string CoupledModelEvaluator<Scalar>::description() const
{
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::CoupledModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

// Private functions overridden from ModelEvaulatorDefaultBase


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

  std::cout << "Got in the CoupledModelEvaluator!" << std::endl;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::CoupledModelEvaluator",inArgs,outArgs
    );

  thyraModel->evalModel(inArgs, outArgs);

  // grab the residual vector and the Jacobian
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f();
  const Teuchos::RCP<Thyra::LinearOpBase<Scalar> > W_op = outArgs.get_W_op();

  //std::cout << "Does ELECTRIC_POTENTIAL exist in block named silicon?: " 
  //          << (globalIndexer_->fieldInBlock("ELECTRIC_POTENTIAL", "silicon")? "yes" : "no") << std::endl;

  // modify the residual vector
  if (nonnull(f)){
    double preValue = 0.0;
    double couplingValue = 0.0;
    for(int index = 0 ; index < (f->space()->dim()) ; ++index){
      // contribute to the index^th entry of the residual
      std::cout << "Bef: f(" << index << ") = " << Thyra::get_ele(*f,index) << std::endl;
      preValue = Thyra::get_ele(*f,index);
      couplingValue = 0.0;
      Thyra::set_ele(index,preValue + couplingValue,f.ptr());
      std::cout << "Aft: f(" << index << ") = " << Thyra::get_ele(*f,index) << std::endl;
    }
  }

  // modify the Jacobian
  if(!(W_op.is_null())){
    Teuchos::RCP<Thyra::EpetraLinearOpBase > W_op_base = 
      rcp_dynamic_cast<Thyra::EpetraLinearOpBase>(W_op);
    const Teuchos::RCP<Epetra_Operator> W_epetra_op = 
      Thyra::get_Epetra_Operator(*W_op);
    Epetra_CrsMatrix &W_crs = Teuchos::dyn_cast<Epetra_CrsMatrix>(*W_epetra_op);
    W_crs.Print(std::cout);
    int numvals = 1;
    double * rowvals = new double [numvals];
    int    * rowinds = new int    [numvals];
    // get these entry values
    //W_crs.ExtractGlobalRowCopy(0, numvals, numvals, rowvals, rowinds); // get the (0,0) entry
    // edit these values
    //rowvals[0] = 1.0;
    // set these values
    //W_crs.ReplaceGlobalValues(0 , // int GlobalRow - (In) Row number (in global coordinates) to put elements.
    //                          numvals , // int NumEntries - (In) Number of entries.
    //                          rowvals , // const double* Values - (In) Values to enter. 
    //                          rowinds); // const int* Indices - (In) Global column indices corresponding to values
    //xW_crs.Print(std::cout);
    delete [] rowvals;
    delete [] rowinds;
  }


  typedef int LocalOrdinal;

  // begin mesh debug info
  {
  std::string name = "cathode"; // provide a sideset id name (e.g., "anode" and "cathode" are the sidesets of the pndiode.dd.equ test

  std::vector<stk::mesh::Entity> side_entities;

  std::vector<std::string> block_names;
  mesh_->getElementBlockNames(block_names);
  std::vector<LocalOrdinal> v_elems;
  std::vector<int> v_nodes;
  std::vector<int> v_nodes_global;

  std::size_t num_nodes = 0;
  for (std::size_t i=0;i<block_names.size(); ++i) {
    auto my_block_name = block_names[i];
    mesh_->getAllSides(name,my_block_name,side_entities); // get all sides of the block specified
    std::vector<stk::mesh::Entity> elements;
    std::vector<std::size_t> local_subcell_ids, subcell_dim;
    panzer_stk::workset_utils::getSideElementCascade(*mesh_, my_block_name,
                                                     side_entities,subcell_dim,local_subcell_ids,elements);

    // count number of nodes on sideset
    for(std::size_t elem = 0; elem < elements.size(); elem++){
      if(subcell_dim[elem] == 0){
        v_elems.push_back(mesh_->elementLocalId(elements[elem]));
        v_nodes.push_back(local_subcell_ids[elem]);
        v_nodes_global.push_back(local_subcell_ids[elem]);
        num_nodes++;
      }
    }
    std::cout << "There are " << num_nodes << " (local) nodes on sideset named: " << name << std::endl;
  }
  Kokkos::View<LocalOrdinal*> elems = Kokkos::View<LocalOrdinal*>("Elems for sideset "+name, num_nodes);
  Kokkos::View<int*> nodes = Kokkos::View<int*>("Nodes for sideset "+name, num_nodes);
  for (std::size_t i=0;i<num_nodes; ++i ){
    elems(i) = v_elems[i];
    nodes(i) = v_nodes[i];
    std::cout << "The " << i << "th " << name << " node belongs to (global id) element " << elems(i) << " at local node number " << nodes(i) << std::endl; 
    //std::cout << name << " element(" << i << ") = " << elems(i) << " (this is actually the global element id as seen in ParaView!) " << std::endl;
    //std::cout << name << " nodes(" << i << ") = " << nodes(i) << " (it seems like these are the local node numbers of this element which lie on the sideset) " << std::endl;
  }

  // print out the jacobian and residual at these sideset nodes


  } // end mesh debug info




  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
}


} // namespace charon

#endif // __Charon_CoupledModelEvaluator_impl_hpp__
