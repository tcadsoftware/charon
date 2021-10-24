
#ifndef CHARON_BCSTRATEGY_DIRICHLET_SCHOTTKYCONTACT_DECL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_SCHOTTKYCONTACT_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Charon_Names.hpp"

namespace charon {

  template <typename EvalT>
  class BCStrategy_Dirichlet_SchottkyContact : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>{

  public:

    BCStrategy_Dirichlet_SchottkyContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& side_pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    Teuchos::RCP<charon::Names> m_names;
    
    //std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;

    private:
    
  };


  

}

#endif
