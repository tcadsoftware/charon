
#ifndef CHARON_NEUMANNBC_THERMALCONTACT_DECL_HPP
#define CHARON_NEUMANNBC_THERMALCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class NeumannBC_ThermalContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    NeumannBC_ThermalContact(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> heat_flux;   // scaled, no unit

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // Temperature scaling [K]
  double H0;
  double X0;

  int num_points;

  double value;
  double temp;

  std::string paramName;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class NeumannBC_ThermalContact


}

#endif
