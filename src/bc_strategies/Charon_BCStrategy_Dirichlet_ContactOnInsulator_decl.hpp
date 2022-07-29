
#ifndef CHARON_BCSTRATEGY_DIRICHLET_CONTACTONINSULATOR_DECL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_CONTACTONINSULATOR_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Charon_Names.hpp"

/**
 * @brief The "Contact On Insulator" boundary condition has the following settings :
 *
 *       <ParameterList>
 *           <Parameter name="Type" type="string" value="Dirichlet"/> 
 *           <Parameter name="Sideset ID" type="string" value="gate"/> 
 *           <Parameter name="Element Block ID" type="string" value="oxide"/> 
 *           <Parameter name="Equation Set Name" type="string" value="ELECTRIC_POTENTIAL"/> 
 *           <Parameter name="Strategy" type="string" value="Contact On Insulator"/>
 *           <ParameterList name="Data">
 *               <Parameter name="Work Function" type="double" value="4.5"/>
 *
 * A time-dependent linear ramp voltage source 
 *               <ParameterList name="Linear Ramp">
 *                   <Parameter name="Initial Time" type="double" value="0.0"/>
 *                   <Parameter name="Initial Voltage" type="double" value="0.0"/>
 *                   <Parameter name="Final Time" type="double" value="3.0"/>
 *                   <Parameter name="Final Voltage" type="double" value="-3.0"/>
 *               </ParameterList>
 * Or a single DC voltage
 *               <Parameter name="Voltage" type="double" value="0.1"/>
 * Or a varying voltage source
 *               <Parameter name="Varying Voltage" type="string" value="Parameter"/>
 *               <Parameter name="Initial Voltage" type="double" value="-5.0" />
 *           </ParameterList>
 *       </ParameterList>
 */

namespace charon {

  template <typename EvalT>
  class BCStrategy_Dirichlet_ContactOnInsulator : public panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT> {

  public:

    BCStrategy_Dirichlet_ContactOnInsulator(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data, Teuchos::RCP<Teuchos::ParameterList> input_pl = Teuchos::rcp(new Teuchos::ParameterList()));

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

  private:

    Teuchos::RCP<charon::Names> m_names;

    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;

    bool isFreqDom = true;
    double small_signal_perturbation;
    bool bLinRamp; 
    bool bTrapezoid;

    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
    Teuchos::RCP<Teuchos::ParameterList> linPList; 
    Teuchos::RCP<Teuchos::ParameterList> trapzPList; 

  };

}

#endif
