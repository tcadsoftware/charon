
#ifndef CHARON_BC_DIRICHLETSCHOTTKYCONTACT_DECL_HPP
#define CHARON_BC_DIRICHLETSCHOTTKYCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::BASIS;

namespace charon {

template<typename EvalT, typename Traits>
class BC_DirichletSchottkyContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_DirichletSchottkyContact(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential;   // scaled, no unit

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;   // [eV]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]

  size_type num_basis;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  double cnt_wf; // eV
 
  Teuchos::RCP<const charon::Names> m_names;
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BC_DirichletSchottkyContact


}

#endif
