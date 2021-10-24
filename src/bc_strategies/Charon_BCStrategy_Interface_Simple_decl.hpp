

#ifndef CHARON_BCSTRATEGY_INTERFACE_SIMPLE_DECL_HPP
#define CHARON_BCSTRATEGY_INTERFACE_SIMPLE_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Interface_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"

namespace charon {

  /** Impose the interface condition
    *     flux = a f_me + b f_other,  (*)
    * where f is a field possibly discontinuous across the interface (2 DOF
    * per interface node), and a,b are scalar coefficients.
    *   ParameterList parameters are as follows:
    *     "Type": "Interface"
    *     "Strategy": "Robin Interface"
    *     "Sideset ID": Name of the interface sideset.
    *     "Element Block ID":  Name of the primary element block, often referred
    *                          to in the code as 'me' or by a similar pronoun.
    *     "Element Block ID2": Name of the element block on the other side of
    *                          the interface, often referred to as 'other'.
    *     "Equation Set Name":  Equation set associated for my element block,
    *                           f_me in (*).
    *     "Equation Set Name2": Equation set for the other's element block,
    *                           f_other in (*).
    *     "Data": The following ParameterList specific to this interface condition:
    *        "a", "b", "c", "d", The coefficients in (*).
    */
  template <typename EvalT>
  class BCStrategy_Interface_Simple : public panzer::BCStrategy_Interface_DefaultImpl<EvalT> {
  public:
    BCStrategy_Interface_Simple(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    virtual void buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                                const panzer::PhysicsBlock& side_pb,
                                                                const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                                const Teuchos::ParameterList& user_data) const;

    virtual void postRegistrationSetup(typename panzer::Traits::SetupData d,
                                       PHX::FieldManager<panzer::Traits>& vm);

    virtual void evaluateFields(typename panzer::Traits::EvalData d);

  private:
    std::string dof_name_, other_dof_name_, coupling_dof_name_, my_coupling_dof_name_;
    bool field_spy_;
    double coeffs_[4];

    static void setCombineValues(Teuchos::ParameterList& p,
                                 const std::string value_name1, const double scalar1,
                                 const std::string value_name2, const double scalar2,
                                 const std::string value_name3 = "", const double scalar3 = 0,
                                 const std::string value_name4 = "", const double scalar4 = 0);
  };

}

#endif
