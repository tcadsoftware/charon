

#ifndef CHARON_BCSTRATEGY_NEUMANN_DYNAMICTRAPS_DECL_HPP
#define CHARON_BCSTRATEGY_NEUMANN_DYNAMICTRAPS_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Neumann_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"


namespace charon {

  template <typename EvalT>
  class BCStrategy_Neumann_DynamicTraps : public panzer::BCStrategy_Neumann_DefaultImpl<EvalT> {

  public:

    BCStrategy_Neumann_DynamicTraps(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& side_pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                   const panzer::PhysicsBlock& side_pb,
                                                   const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                   const Teuchos::ParameterList& user_data) const;

    virtual void postRegistrationSetup(typename panzer::Traits::SetupData d,
                                       PHX::FieldManager<panzer::Traits>& vm);

    virtual void evaluateFields(typename panzer::Traits::EvalData d);

  private:
    
    void initDynamicTrapsParams(Teuchos::RCP<const Teuchos::ParameterList>);
    
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
    
    Teuchos::RCP<Teuchos::ParameterList> dynTrapsPList;

    std::string fluxDynTrapsCharge, eFluxDynTrapsRecomb, hFluxDynTrapsRecomb;
    bool withField;
  };

}

#endif
