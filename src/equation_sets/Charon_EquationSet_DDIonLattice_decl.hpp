
#ifndef CHARON_EQUATIONSET_DDIONLATTICE_DECL_HPP
#define CHARON_EQUATIONSET_DDIONLATTICE_DECL_HPP

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
  class EquationSet_DDIonLattice : public charon::EquationSet_DefaultImpl<EvalT> {

  public:

    EquationSet_DDIonLattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
            const int& default_integration_order,
            const panzer::CellData& cell_data,
            const Teuchos::RCP<panzer::GlobalData>& global_data,
            const bool build_transient_support);

    EquationSet_DDIonLattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
            const int& default_integration_order,
            const panzer::CellData& cell_data,
            const Teuchos::RCP<panzer::GlobalData>& global_data,
            const bool build_transient_support,
            const std::string& dummy);

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
            const panzer::FieldLibrary& field_library,
            const Teuchos::ParameterList& user_data) const;

    void computeSUPGStabResidual(PHX::FieldManager<panzer::Traits>& fm,
            const panzer::FieldLibrary& fl,
            const Teuchos::ParameterList& user_data,
            const Teuchos::RCP<charon::Names>& m_names,
            const Teuchos::RCP<panzer::IntegrationRule>& ir,
            const Teuchos::RCP<panzer::BasisIRLayout>& basis,
            const std::string& ls_type, const std::string& tau_type,
            const bool& add_source_stab, const std::string& carr_type,
            const int& ion_charge, const bool& includeSoret) const;

  private:

    Teuchos::RCP<charon::Names>  m_names;

    // Options
    std::string solveElectron;
    std::string solveHole;
    std::string solveIon;

    std::string stab_scheme;
    std::string tau_e_type;
    std::string tau_h_type;
    std::string tau_ion_type;
    std::string ls_type;

    bool haveSource;
    bool add_source_stab;
    bool haveBGN;

    std::string field_model;

    int ion_charge;

  };

}

#endif
