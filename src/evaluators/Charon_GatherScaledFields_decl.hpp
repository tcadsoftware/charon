
#ifndef _CHARON_GatherScaledFields_DECL_HPP_
#define _CHARON_GatherScaledFields_DECL_HPP_

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_Interface.hpp"


namespace charon
{

/**
 * A class that will extract values from the STK data structures and put
 * them in the local, Panzer, data structures after applying
 * scaling. This is used, for example, when reading an initial guess
 * from an Exodus file and scaling it prioer to using it.
 *
 * Note that this is almost an exact replica of the corresponding Panzer
 * class, <code>panzer_stk::GatherFields</code>.
 */
template<typename EvalT, typename Traits>
class GatherScaledFields :
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>
{

public:

  GatherScaledFields(Teuchos::RCP<panzer_stk::STK_Interface const> const& mesh,
                     Teuchos::ParameterList const& p);

  void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& vm);

  void evaluateFields(typename Traits::EvalData d);

private:
  typedef typename EvalT::ScalarT ScalarT;
  typedef panzer_stk::STK_Interface::SolutionFieldType VariableField;

  std::vector< PHX::MDField<ScalarT,panzer::Cell,panzer::NODE> > gatherFields_;
  std::vector<VariableField*> stkFields_;

  Teuchos::RCP<panzer_stk::STK_Interface const> mesh_;

  // Disabled default ctor for external classes
  GatherScaledFields();

  //! A map of scale factor by variable name.
  Teuchos::RCP<std::map<std::string,double> > scaleFactors_;

};

} // namespace charon

#endif // _CHARON_GatherScaledFields_DECL_HPP_
