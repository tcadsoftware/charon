
#ifndef __Charon_ResponseEvaluatorFactory_HOCurrent_hpp__
#define __Charon_ResponseEvaluatorFactory_HOCurrent_hpp__

#include <string>


#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"

#include "Charon_Names.hpp"

namespace charon {

/** This class defines a response based on a functional. */

template <typename EvalT,typename LO,typename GO>
class ResponseEvaluatorFactory_HOCurrent : public panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>
{
public:

   ResponseEvaluatorFactory_HOCurrent(MPI_Comm comm, int cubatureDegree,
                                      const Teuchos::RCP<const charon::Scaling_Parameters> & scaleParams,
                                      const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory=Teuchos::null,
                                      std::string fd_suffix = "",bool isFreqDom = false)
     : panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>(comm,cubatureDegree,false,"",linearObjFactory)
     , scaleParams_(scaleParams)
     , fd_suffix_(fd_suffix)
     , isFreqDom_(isFreqDom)
   {
     TEUCHOS_ASSERT(scaleParams_!=Teuchos::null);

     names = Teuchos::rcp(new charon::Names(1,"","","",fd_suffix_));  // Prefix must be "" in Physics Blocks
   }

   virtual ~ResponseEvaluatorFactory_HOCurrent() {}

   /** Build and register evaluators for a response on a particular physics
     * block.
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
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
  Teuchos::RCP<const charon::Scaling_Parameters> scaleParams_;

  std::string fd_suffix_;

  bool isFreqDom_;

  Teuchos::RCP<charon::Names> names;
};

template <typename LO,typename GO>
struct HOCurrentResponse_Builder {
  MPI_Comm comm;
  int cubatureDegree;
  std::string fd_suffix = "";
  bool isFreqDom = false;
  Teuchos::RCP<const charon::Scaling_Parameters> scaleParams;
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer;

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_HOCurrent<T,LO,GO>(comm,cubatureDegree,scaleParams,linearObjFactory,fd_suffix,isFreqDom)); }
};

}

#include "Charon_ResponseEvaluatorFactory_HOCurrent_impl.hpp"

#endif
