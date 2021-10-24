
#ifndef CHARON_BCSTRATEGY_FREQDOM_DECL_HPP
#define CHARON_BCSTRATEGY_FREQDOM_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_Traits.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Charon_Names.hpp"
#include "Charon_FreqDom_Parameters.hpp"

namespace charon {

  template <typename EvalT>
  class BCStrategy_FreqDom : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {
    
  public:    
    
    BCStrategy_FreqDom(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);
    
    void setup(const panzer::PhysicsBlock& side_pb,
	       const Teuchos::ParameterList& user_data);
    
    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				    const panzer::PhysicsBlock& pb,
				    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				    const Teuchos::ParameterList& models,
				    const Teuchos::ParameterList& user_data) const;

    std::vector<Teuchos::RCP<charon::Names> > m_fd_names;

    Teuchos::RCP<FreqDomParameters> freqDomParamsRCP;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;

    bool isLatTDof; 
    bool isIonDof; 
    bool isFermiPin; 

    Teuchos::RCP<std::vector<std::string> > fd_phi_target_cos_names;
    Teuchos::RCP<std::vector<std::string> > fd_phi_target_sin_names;
    Teuchos::RCP<std::vector<std::string> > fd_elec_target_cos_names;
    Teuchos::RCP<std::vector<std::string> > fd_elec_target_sin_names;
    Teuchos::RCP<std::vector<std::string> > fd_hole_target_cos_names;
    Teuchos::RCP<std::vector<std::string> > fd_hole_target_sin_names;

    Teuchos::RCP<std::vector<std::string> > td_phi_target_names;
    Teuchos::RCP<std::vector<std::string> > td_elec_target_names;
    Teuchos::RCP<std::vector<std::string> > td_hole_target_names;

    // parameters to be passed to the time domain BC Strategy constructors
    Teuchos::RCP<const panzer::BC> td_bc;
    Teuchos::RCP<panzer::GlobalData> td_global_data;
    
  };

}

#endif
