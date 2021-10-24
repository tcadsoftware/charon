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

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

// Panzer
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_GlobalIndexer.hpp"

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
  
  /** \brief Constructs to uninitialized */
  CoupledModelEvaluator(
        const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
        MPI_Comm                                           comm,
        const Teuchos::RCP<Teuchos::ParameterList>&        parameters,
        const Teuchos::RCP<panzer_stk::STK_Interface>&     mesh);
  
  /** \brief . */
  std::string description() const;

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}
  
private:

   /*  \brief The physics `ModelEvaluator`.
    */
   Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> physics_;

  /**                                                                                                                       
   *  \brief The MPI communicator to use with this `ModelEvaluator`.
   */
  Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal>> comm_;

  Teuchos::RCP<Teuchos::ParameterList> parameters_;
  const int dimension_;
  const Teuchos::RCP<panzer_stk::STK_Interface> mesh_;

  Teuchos::RCP<panzer::GlobalIndexer> globalIndexer_;

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
