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

#ifndef __Charon_CoupledModelEvaluator_hpp__
#define __Charon_CoupledModelEvaluator_hpp__

#include "Charon_config.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

// Panzer
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "Charon_CurrentConstraintList.hpp"

#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
//#include <N_CIR_Xyce.h>
#include <N_CIR_SecondLevelSimulator.h>
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER

namespace charon {


/** \brief This class decorates a ModelEvaluator and returns scaled
 * residual and Jacobian values.
 *
 * Given a scaling vector <tt>s</tt>, this object is treated as a diagonal
 * scaling matrix and applied to <tt>x -> Sf(x)</tt> and <tt>x -> sW</tt>.
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class CoupledModelEvaluator : 
    virtual public Thyra::ModelEvaluatorDelegatorBase<Scalar>
{
public:
  
  // constructor
  CoupledModelEvaluator(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& physics,
        MPI_Comm                                            comm,
        const Teuchos::RCP<Teuchos::ParameterList>&         parameters,
        charon::CurrentConstraintList&                      constraints,
        const Teuchos::RCP<panzer_stk::STK_Interface>&      mesh);

  // returns description of the model evaluator,
  // used  when composing mol evaluators and checekingfor compatibilityy
  std::string description() const;

  // destructor
  ~CoupledModelEvaluator()
  {
#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
    delete simulator_;
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER
  }

private:

  // the core of this model evaluator
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> physics_;
  Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal>> comm_;

  Teuchos::RCP<Teuchos::ParameterList> parameters_;
  const int dimension_;
  const Teuchos::RCP<panzer_stk::STK_Interface> mesh_;
  Teuchos::RCP<panzer::GlobalIndexer> globalIndexer_;

  /// The list of all the current constraints. 
  charon::CurrentConstraintList constraints_;
  Teuchos::RCP<panzer::ParamLib> paramLib_;
  Teuchos::RCP<Teuchos::ParameterList> xyceCouplingPL_;
  bool mixedModeLOCA = false;
  bool printDebug = false;
  bool printVerbose = true;

#if defined(ENABLE_MIXED_MODE) or defined(ENABLE_XYCE_CLUSTER)
  Xyce::Circuit::SecondLevelSimulator * simulator_;
#endif // ENABLE_MIXED_MODE || ENABLE_XYCE_CLUSTER

};


/** \brief Nonmember constructor. */
template<class Scalar>
Teuchos::RCP<CoupledModelEvaluator<Scalar> >
createCoupledModelEvaluator(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar > > &model)
{
  Teuchos::RCP<CoupledModelEvaluator<Scalar> > coupledME(new CoupledModelEvaluator<Scalar>);
  coupledME->initialize(model);
  return coupledME;
}


} // namespace charon


#endif // CHARON_COUPLEDMODELEVALUATOR_HPP
