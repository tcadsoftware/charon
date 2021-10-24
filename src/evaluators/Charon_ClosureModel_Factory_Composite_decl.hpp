#ifndef CHARON_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP
#define CHARON_CLOSURE_MODEL_FACTORY_COMPOSITE_DECL_HPP

#include "Panzer_ClosureModel_Factory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"

namespace charon {

  template<typename EvalT>
  class ClosureModelFactoryComposite : public panzer::ClosureModelFactory<EvalT> {

  public:

    ClosureModelFactoryComposite(const std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >& factories);

    Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
    buildClosureModels(const std::string& model_id,
		       const Teuchos::ParameterList& models,
		       const panzer::FieldLayoutLibrary& fl,
		       const Teuchos::RCP<panzer::IntegrationRule>& ir,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& fm) const;

  private:

    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > > m_factories;

  };

}

#endif
