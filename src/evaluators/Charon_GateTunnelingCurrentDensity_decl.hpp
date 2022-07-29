
#ifndef CHARON_GATETUNNELINGCURRENTDENSITY_DECL_HPP
#define CHARON_GATETUNNELINGCURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
//#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {


template<typename EvalT, typename Traits>
class GateTunnelingCurrentDensity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    GateTunnelingCurrentDensity(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

  private:
    using ScalarT = typename EvalT::ScalarT;

   
    // output fields
    PHX::MDField<ScalarT,Cell,Point> tunnel_current;   // scaled, no unit

    // input fields
    PHX::MDField<const ScalarT,Cell,Point> carr_dens;
    PHX::MDField<const ScalarT,Cell,Point> pot;
    PHX::MDField<const ScalarT,Cell,Point> latt_temp;

    // scaling parameters
    Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
    // double T0; // temperature scaling [K]
    double C0; // concentration scaling [cm^-3]
    //double J0; // current density scaling [A/cm^2]

    int num_ips, num_nodes; 
    std::string basis_name; 
    std::size_t basis_index;
 
    std::string sidesetID;
    std::string gate_sidesetID;
    double gate_dist;
    std::string blockID;
    

    std::string tunnel_current_name; 

  
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

 
    
 
}; // end of class GateTunnelingCurrentDensity


}

#endif
