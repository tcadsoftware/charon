

#ifndef CHARON_BCSTRATEGY_NEUMANN_SCHOTTKYCONTACT_DECL_HPP
#define CHARON_BCSTRATEGY_NEUMANN_SCHOTTKYCONTACT_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Neumann_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace charon {

  template <typename EvalT>
  class BCStrategy_Neumann_SchottkyContact : public panzer::BCStrategy_Neumann_DefaultImpl<EvalT> {

  public:

    BCStrategy_Neumann_SchottkyContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& side_pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                   const panzer::PhysicsBlock& side_pb,
                                                   const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                   const Teuchos::ParameterList& user_data) const;

    virtual void postRegistrationSetup(typename panzer::Traits::SetupData d,
                                       PHX::FieldManager<panzer::Traits>& vm);

    virtual void evaluateFields(typename panzer::Traits::EvalData d);

  private:
    int contact_type; // -1 n-type, 1 n-type
    bool withBL;
    bool withTunneling;
    double An;
    double Ap;
    double Wf;
    double alpha;
    double beta;
    double gamma;
    double m_tun;
  };

}

#endif
