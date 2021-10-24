
#ifndef CHARON_EQUATIONSET_NLPOISSON_DECL_HPP
#define CHARON_EQUATIONSET_NLPOISSON_DECL_HPP

#include <vector>
#include <string>

#include "Charon_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Charon_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Charon_Names.hpp"

namespace charon {

  template <typename EvalT>
  class EquationSet_NLPoisson : public charon::EquationSet_DefaultImpl<EvalT> {

  public:

    EquationSet_NLPoisson(const Teuchos::RCP<Teuchos::ParameterList>& params,
                          const int& default_integration_order,
                          const panzer::CellData& cell_data,
                          const Teuchos::RCP<panzer::GlobalData>& global_data,
                          const bool build_transient_support);

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::FieldLibrary& field_library,
                                               const Teuchos::ParameterList& user_data) const;

  protected:

    Teuchos::RCP<charon::Names>  m_names;

  };

}

#endif
