
#ifndef __Charon_ResponseEvaluatorFactory_DispCurrent_hpp__
#define __Charon_ResponseEvaluatorFactory_DispCurrent_hpp__

#include <string>


#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

namespace charon {

/** This class defines a response based on a functional. */

template <typename EvalT,typename LO,typename GO>
class ResponseEvaluatorFactory_DispCurrent : public panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>
{
public:

   ResponseEvaluatorFactory_DispCurrent(bool isTransient, bool isSingleFreq, MPI_Comm comm, int cubatureDegree,
                                        const Teuchos::RCP<const charon::Scaling_Parameters> & scaleParams,
                                        const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory=Teuchos::null,
                                        std::string fd_suffix = "", bool isFreqDom = false, double omega = 0.0)
     : panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>(comm,cubatureDegree,false,"",linearObjFactory)
     , scaleParams_(scaleParams)
     , fd_suffix_(fd_suffix)
     , isFreqDom_(isFreqDom)
   {
     isSingleFreq_ = isSingleFreq; 
     omega_ = omega; // radian frequency in physical unit
     isTransient_ = isTransient; 
   }

   virtual ~ResponseEvaluatorFactory_DispCurrent() {}

   /** Build and register evaluators for a response on a particular physics
     * block.
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be built with the evaluators.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     */
   virtual void buildAndRegisterEvaluators(const std::string & responseName,
                                           PHX::FieldManager<panzer::Traits> & fm,
                                           const panzer::PhysicsBlock & physicsBlock,
                                           const Teuchos::ParameterList & user_data) const;

   /** Is this evaluation type supported by the factory. This is used to determine cases
     * where a response may support a particular evaluation type, however at runtime the user
     * decides not to enable the (say) Jacobian evaluation of this response.
     *
     * Note that use of this mechanism is complementary to having the builder return
     * <code>Teuchos::null</code> for a particular evaluation type.
     */
   virtual bool typeSupported() const;

private:
  // std::vector<std::string> currents_;

  Teuchos::RCP<const charon::Scaling_Parameters> scaleParams_;

  std::string fd_suffix_;
  bool isFreqDom_;
  
  bool isSingleFreq_;
  double omega_; 
  
  bool isTransient_;

  Teuchos::RCP<charon::Names> names;
};

template <typename LO,typename GO>
struct DispCurrentResponse_Builder 
{
  MPI_Comm comm;
  bool isTransient;
  bool isSingleFreq;
  bool isFreqDom;
  double omega; 
  int cubatureDegree;
  std::string fd_suffix;
 
  Teuchos::RCP<const charon::Scaling_Parameters> scaleParams;
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer;

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_DispCurrent<T,LO,GO>    (isTransient,isSingleFreq,comm,cubatureDegree,scaleParams,linearObjFactory,fd_suffix,isFreqDom,omega)); }
};

}

#include "Charon_ResponseEvaluatorFactory_DispCurrent_impl.hpp"

#endif
