#ifndef CHARON_CLOSUREMODEL_FACTORY_COMPOSITE_TEMPLATE_BUILDER_HPP
#define CHARON_CLOSUREMODEL_FACTORY_COMPOSITE_TEMPLATE_BUILDER_HPP

#include <string>
#include "Sacado_mpl_apply.hpp"
#include "Teuchos_RCP.hpp"
#include "Charon_ClosureModel_Factory_Composite.hpp"

namespace charon {

  class ClosureModelFactoryComposite_TemplateBuilder {

    std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > > m_factories;

  public:

    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
			   (new charon::ClosureModelFactoryComposite<EvalT>(m_factories)) );
    }
    
    void addFactory(const Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> >& factory)
    {
      m_factories.push_back(factory);
    }

  };
  
}

#endif 
