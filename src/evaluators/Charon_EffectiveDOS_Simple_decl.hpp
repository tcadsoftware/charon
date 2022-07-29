
#ifndef CHARON_EFFECTIVEDOS_SIMPLE_DECL_HPP
#define CHARON_EFFECTIVEDOS_SIMPLE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_Material_Properties.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

//! calculate the effective dos using Nc\v = Nc\v(300)*(T/300)^Nc\v_F

template<typename EvalT, typename Traits>
class EffectiveDOS_Simple
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    EffectiveDOS_Simple(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<ScalarT,Cell,Point> hole_effdos;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> xMoleFrac;
  PHX::MDField<const ScalarT,Cell,Point> yMoleFrac;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;
  double C0;

  int num_points;

  // material parameters
  double Nc300, Nv300, Nc_F, Nv_F;

  std::string materialName;
  bool withMoleFrac;
  Teuchos::RCP<CompoundMaterial> comp_mat;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class EffectiveDOS_Simple


}

#endif
