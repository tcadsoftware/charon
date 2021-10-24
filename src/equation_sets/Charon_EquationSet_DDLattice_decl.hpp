
#ifndef CHARON_EQUATIONSET_DDLATTICE_DECL_HPP
#define CHARON_EQUATIONSET_DDLATTICE_DECL_HPP

#include <vector>
#include <string>

#include "Charon_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Charon_EquationSet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Charon_Names.hpp"

// Include the DDIonLattice file to use the computeSUPGStabResidual(.) function
#include "Charon_EquationSet_DDIonLattice.hpp"


namespace charon {

  template <typename EvalT>
  class EquationSet_DDLattice : public charon::EquationSet_DefaultImpl<EvalT> {

  public:

    EquationSet_DDLattice(const Teuchos::RCP<Teuchos::ParameterList>& params,
                          const int& default_integration_order,
                          const panzer::CellData& cell_data,
                          const Teuchos::RCP<panzer::GlobalData>& global_data,
                          const bool build_transient_support);

    void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::FieldLibrary& field_library,
                                               const Teuchos::ParameterList& user_data) const;

  protected:

    Teuchos::RCP<charon::Names>  m_names;

    // Options
    std::string solveElectron;
    std::string solveHole;
    std::string field_model;

    std::string supg_stab;
    std::string tau_e_type;
    std::string tau_h_type;
    std::string ls_type;

    bool haveSource;
    bool add_source_stab;
    bool haveBGN;

    // discretization method = FEM-SUPG, FEM-EFFPG, or CVFEM-SG
    std::string discMethod;

    // Define an EquationSet_DDIonLattice RCP to call its computeSUPGStabResidual(.) member function
    Teuchos::RCP<charon::EquationSet_DDIonLattice<EvalT> >  eqnSetDDIonLatt;

  };

}

#endif
