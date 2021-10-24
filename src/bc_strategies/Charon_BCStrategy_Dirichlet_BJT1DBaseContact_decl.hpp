
#ifndef CHARON_BCSTRATEGY_DIRICHLET_BJT1DBASECONTACT_DECL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_BJT1DBASECONTACT_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Charon_Names.hpp"

namespace charon {

  template <typename EvalT>
  class BCStrategy_Dirichlet_BJT1DBaseContact : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
    
  public:    
    
    //BCStrategy_Dirichlet_BJT1DBaseContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    BCStrategy_Dirichlet_BJT1DBaseContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data, Teuchos::RCP<Teuchos::ParameterList> input_pl = Teuchos::rcp(new Teuchos::ParameterList()));

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    Teuchos::RCP<charon::Names> m_names;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;
    
    bool isFreqDom = true;

  };

}

#endif
