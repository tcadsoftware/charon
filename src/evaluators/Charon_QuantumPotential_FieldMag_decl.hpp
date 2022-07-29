
#ifndef CHARON_QUANTUMPOTENTIAL_FIELDMAG_DECL_HPP
#define CHARON_QUANTUMPOTENTIAL_FIELDMAG_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
//using panzer::Point;
using panzer::IP;
using panzer::Dim;

namespace charon {

template<typename EvalT, typename Traits>
class QuantumPotentialFieldMag
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    QuantumPotentialFieldMag(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

    // output
    PHX::MDField<ScalarT,Cell,IP> qp_fieldmag;
    // input
    PHX::MDField<const ScalarT,Cell,IP,Dim> grad_phi;
    PHX::MDField<const ScalarT,Cell,IP,Dim> grad_qp;
    PHX::MDField<const ScalarT,Cell,IP> latt_temp;
  
    double fitParam;
  
    //std::size_t num_points;
    std::size_t num_ip;
    std::size_t num_dim;
  
  private:
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class QuantumPotentialFieldMag


}

#endif
