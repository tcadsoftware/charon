
#ifndef CHARON_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP
#define CHARON_CLOSURE_MODEL_FACTORY_TEMPLATE_BUILDER_HPP

#include <string>
#include "Teuchos_RCP.hpp"
#include "Charon_ClosureModel_Factory.hpp"
#include "Charon_Scaling_Parameters.hpp"

namespace charon {

  class ClosureModelFactory_TemplateBuilder {
    Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;
    bool m_throw_if_model_not_found;
    std::string m_model_name_prefix;
    std::string m_fd_suffix;

  public:
    ClosureModelFactory_TemplateBuilder(const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams,
                                        bool throwIfModelNotFound=true, const std::string & modelNamePrefix="Charon Parameters->", const std::string& fd_suffix = "")
       : m_scale_params(scaleParams), m_throw_if_model_not_found(throwIfModelNotFound), m_model_name_prefix(modelNamePrefix), m_fd_suffix(fd_suffix)  {}

    template <typename EvalT>
    Teuchos::RCP<panzer::ClosureModelFactoryBase> build() const {
      return Teuchos::rcp( static_cast<panzer::ClosureModelFactoryBase*>
			   (new charon::ClosureModelFactory<EvalT>(m_scale_params, m_throw_if_model_not_found, m_model_name_prefix, m_fd_suffix)) );
    }

  };

}

#endif
