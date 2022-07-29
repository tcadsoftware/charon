#ifndef CHARON_RELATIVE_PERMITTIVITY_DECL_HPP
#define CHARON_RELATIVE_PERMITTIVITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Material_Properties.hpp"


using panzer::Cell;
using panzer::Point;

namespace charon {

template<typename EvalT, typename Traits>
class Relative_Permittivity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  using ScalarT = typename EvalT::ScalarT;

public:
  Relative_Permittivity(const Teuchos::ParameterList& p);

  void evaluateFields(typename Traits::EvalData d);

private:
  // output
  PHX::MDField<ScalarT,Cell,Point> rel_perm; 

  // input
  PHX::MDField<const ScalarT,Cell,Point> xMoleFrac;
  PHX::MDField<const ScalarT,Cell,Point> yMoleFrac;

  int num_points;

  std::string materialName;
  bool withMoleFrac;
  Teuchos::RCP<CompoundMaterial> comp_mat;
  
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
PHX_EVALUATOR_CLASS_END

}

#endif
