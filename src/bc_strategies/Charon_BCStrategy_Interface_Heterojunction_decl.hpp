

#ifndef CHARON_BCSTRATEGY_INTERFACE_HETEROJUNCTION_DECL_HPP
#define CHARON_BCSTRATEGY_INTERFACE_HETEROJUNCTION_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Interface_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"


namespace charon {

/**
 * @brief This bcstrategy code implements the thermionic emission (TE) and local tunneling (LT)
 * boundary conditions at a heterojunction (HJ).

 * Panzer requires a HJ to be defined as a one-sided sideset in Cubit. Discontinuous
 * DOFs on the two sides of the HJ (including HJ) are defined by a DOF name and
 * a suffix number, e.g., ELECTRON_DENSITY1 and ELECTRON_DENSITY2. Due to the limited
 * HJ infrastructure, the implemented TE and LT models require that the left side of
 * a HJ always corresponds to side 1 (with suffix being 1), while the right side
 * always corresponds to side 2 (with suffix being 2).

 * The TE and LT model requires the normal current density at a HJ is given by
 * J_{HJ} = J_{TE} (1+delta), where J_{TE} is determined by the TE model implemented
 * in charon::Heterojunction_CurrentDensity, while (1+delta) is determined by the LT
 * model implemented in charon::Heterojunction_LocalTunneling. The J_{HJ} flux adds
 * contribution to the residual equations of the DOFs at both sides of the HJ.
 * The models are implemented for both FEM-SUPG and CVFEM-SG discretization schemes.
*/

  template <typename EvalT>
  class BCStrategy_Interface_Heterojunction : public panzer::BCStrategy_Interface_DefaultImpl<EvalT> {

  public:

    BCStrategy_Interface_Heterojunction(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

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
    double richConst;
    double bandOffset;
    double fdDensity;

    std::string other_dof_name;
    std::string electric_potential_name;
    std::string indexer_electric_potential_name;
    std::string discMethod;

    const std::string femsupg = "FEM-SUPG";
    const std::string cvfemsg = "CVFEM-SG";

    bool localTunnel;
    double tunnelMass;
    bool field_spy_;

    double surfCharge;
  };

}

#endif
