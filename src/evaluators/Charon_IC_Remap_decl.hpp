
#ifndef CHARON_IC_REMAP_DECL_HPP
#define CHARON_IC_REMAP_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;

// Set a field to equal the value of another field
namespace charon {

template<typename EvalT, typename Traits>
class IC_Remap
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IC_Remap(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,BASIS> output_field;

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> input_field;

  std::string dof_name;
  std::string dof_name_input;

  Teuchos::RCP<const charon::Names>  m_names;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class IC_Remap


}

#endif
