
#ifndef CHARON_EQUATIONSET_DRIFTDIFFUSION_DECL_HPP
#define CHARON_EQUATIONSET_DRIFTDIFFUSION_DECL_HPP

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
  class EquationSet_DriftDiffusion : public charon::EquationSet_DefaultImpl<EvalT> {

  public:

    EquationSet_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
                               const int& default_integration_order,
                               const panzer::CellData& cell_data,
                               const Teuchos::RCP<panzer::GlobalData>& global_data,
                               const bool build_transient_support);

    using charon::EquationSet_DefaultImpl<EvalT>::EquationSet_DefaultImpl;

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::FieldLibrary& field_library,
                                               const Teuchos::ParameterList& user_data) const;


  protected:

    // Options
    Teuchos::RCP<charon::Names> m_names;

    std::string solveElectron;
    std::string solveHole;

    std::string supg_stab;
    std::string tau_e_type;
    std::string tau_h_type;
    std::string ls_type;

    bool haveSource;
    bool add_source_stab;
    bool addTrapCharge; 

  };

}

#endif
