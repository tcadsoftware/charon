
#ifndef __Charon_ResponseEvaluatorFactory_Current_hpp__
#define __Charon_ResponseEvaluatorFactory_Current_hpp__

#include <string>


#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"

#include "Charon_Names.hpp"

namespace charon {

/** This class defines a response based on a functional. */

template <typename EvalT,typename LO,typename GO>
class ResponseEvaluatorFactory_Current : public panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>
{
public:

   ResponseEvaluatorFactory_Current(bool enableHoles, bool enableElectrons, MPI_Comm comm, int cubatureDegree,
                                    const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory=Teuchos::null,
                                    std::string fd_suffix = "", bool isFreqDom = false)
     : panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>(comm,cubatureDegree,false,"",linearObjFactory)
   {
     names = Teuchos::rcp(new charon::Names(1,"","","",fd_suffix));  // Prefix must be "" in Physics Blocks

   //  j0Name_ = names->field.J0;
   //  x0Name_ = names->field.X0;

     enableElectrons_ = enableElectrons;
     enableHoles_ = enableHoles;

     isFreqDom_ = isFreqDom;

     if(enableHoles)
       currents_.push_back(names->field.hole_contact_currdens);
     if(enableElectrons)
       currents_.push_back(names->field.elec_contact_currdens);

     TEUCHOS_TEST_FOR_EXCEPTION(currents_.size()==0,std::logic_error,
                     "Charon Error!: Current calculation requires either holds or electrons to be specified.");
   }

   virtual ~ResponseEvaluatorFactory_Current() {}

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
  std::vector<std::string> currents_;

//  std::string j0Name_;
//  std::string x0Name_;

  bool enableElectrons_;
  bool enableHoles_;

  Teuchos::RCP<charon::Names> names;

  bool isFreqDom_;

};

template <typename LO,typename GO>
struct CurrentResponse_Builder {
  MPI_Comm comm;
  bool enableHoles;
  bool enableElectrons;
  int cubatureDegree;
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer;
  std::string fd_suffix;
  bool isFreqDom;

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_Current<T,LO,GO>(enableHoles,enableElectrons,comm,cubatureDegree,linearObjFactory,fd_suffix,isFreqDom)); }
};

}

#include "Charon_ResponseEvaluatorFactory_Current_impl.hpp"

#endif
