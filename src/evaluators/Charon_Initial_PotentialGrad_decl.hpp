
#ifndef CHARON_SAVE_GRADPOTENTIAL_DECL_HPP
#define CHARON_SAVE_GRADPOTENTIAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Doping_Params.hpp"
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Charon_MMS_AnalyticFunctions.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include <ostream>

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using panzer::Dim; 

namespace charon {

//! obtain a uniform doping
template<typename EvalT, typename Traits>
class Initial_PotentialGrad
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Initial_PotentialGrad(const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  void preEvaluate(typename Traits::PreEvalData d);

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> potential; 
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi; 

  // output
  PHX::MDField<ScalarT,Cell,BASIS> initial_phi; 
  PHX::MDField<ScalarT,Cell,IP,Dim> initial_grad_phi; 

  // Fields for each workset to store values
  std::vector<PHX::MDField<ScalarT,Cell,BASIS> > initial_phi_wkst; 
  std::vector<PHX::MDField<ScalarT,Cell,IP,Dim> > initial_grad_phi_wkst;
  std::vector<bool> bSaveField_wkst; 

  // workset id
  int worksetId;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ip;
  int num_dim;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;
  int num_basis; 

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Initial_PotentialGrad


}

#endif
