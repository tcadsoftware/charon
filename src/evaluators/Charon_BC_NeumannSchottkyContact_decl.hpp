
#ifndef CHARON_BC_NEUMANNSCHOTTKYCONTACT_DECL_HPP
#define CHARON_BC_NEUMANNSCHOTTKYCONTACT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class BC_NeumannSchottkyContact
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BC_NeumannSchottkyContact(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  // output
  PHX::MDField<ScalarT,Cell,Point> eSurfCurrent;  // scaled
  PHX::MDField<ScalarT,Cell,Point> hSurfCurrent;  // scaled
  
  // input
  PHX::MDField<const ScalarT,Cell,Point> edensity;    // scaled
  PHX::MDField<const ScalarT,Cell,Point> hdensity;    // scaled
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap; // eV
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos; // scaled
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos; // scaled
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;   // scaled
  PHX::MDField<const ScalarT,Cell,Point> effChi;         // eV
  PHX::MDField<const ScalarT,Cell,Point> rel_perm;
  PHX::MDField<const ScalarT,Cell,Point> EdotNorm; 

  double C0, T0, J0, E0;  // scaling parameters

  int num_ips, num_nodes; 
  std::string basis_name; 
  std::size_t basis_index;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  
  int cnt_type; // -1 n-type, 1 n-type
  double An;
  double Ap;
  double Wf;
  bool withBL;
  double BL_alpha;
  double BL_beta;
  double BL_gamma;
  bool withTunneling;
  double tun_m;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class BC_NeumannSchottkyContact


}

#endif
