
#ifndef CHARON_HETEROJUNCTION_SURFACECHARGE_DECL_HPP
#define CHARON_HETEROJUNCTION_SURFACECHARGE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class Heterojunction_SurfaceCharge
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Heterojunction_SurfaceCharge(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

  private:
    using ScalarT = typename EvalT::ScalarT;

    // output fields
    PHX::MDField<ScalarT,Cell,Point> surface_charge;   // scaled, no unit
    
    // input fields

    double C0, X0;  // scaling parameters

    Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > fixedCharge; // in unit of cm^{-2}

    int num_ips, num_nodes; 
    std::string basis_name; 
    std::size_t basis_index;
 
    std::string fluxSurfCharge; 
 
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Heterojunction_SurfaceCharge


}

#endif
