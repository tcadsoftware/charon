
#ifndef CHARON_NEUMANNBC_DYNAMICTRAPS_DECL_HPP
#define CHARON_NEUMANNBC_DYNAMICTRAPS_DECL_HPP

#include "Charon_Scaling_Parameters.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_RecombRate_DynamicTraps.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class NeumannBC_DynamicTraps
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
public:
    NeumannBC_DynamicTraps(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

private:
  using ScalarT = typename EvalT::ScalarT;

  double initial_time;
  double prev_time;

  // output
  PHX::MDField<ScalarT,Cell,Point> fluxCharge;
  PHX::MDField<ScalarT,Cell,Point> eFluxRecomb;
  PHX::MDField<ScalarT,Cell,Point> hFluxRecomb;
  
  // input 
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> e_gamma;
  PHX::MDField<const ScalarT,Cell,Point> h_gamma;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;
  PHX::MDField<const ScalarT,Cell,Point> elec_field;

  // traps related
  Kokkos::DynRankView<ScalarT,PHX::Device> edens;
  Kokkos::DynRankView<ScalarT,PHX::Device> hdens;
  Kokkos::DynRankView<ScalarT,PHX::Device> lT;
  Kokkos::DynRankView<ScalarT,PHX::Device> egamma;
  Kokkos::DynRankView<ScalarT,PHX::Device> hgamma;
  Kokkos::DynRankView<ScalarT,PHX::Device> e_effdos;
  Kokkos::DynRankView<ScalarT,PHX::Device> h_effdos;
  Kokkos::DynRankView<ScalarT,PHX::Device> eff_bg;
  Kokkos::DynRankView<ScalarT,PHX::Device> field;

  // traps
  Teuchos::RCP<DynamicTraps<EvalT>> traps;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]
  double T0; // temperature scaling, [K]
  double X0; // length scaling, [cm]
  double C0; // concentration scaling, [cm^(-3)]
  double R0; // recomb./gen. scaling,[#/(cm^3.s)]
  double E0;

  int num_ips;
  int num_nodes;
  int num_dims;
  int int_rule_degree;
  std::size_t int_rule_index;
  std::size_t basis_index;
  std::string basis_name; 

  bool withField;

  std::string fluxDynTrapsCharge, eFluxDynTrapsRecomb, hFluxDynTrapsRecomb; 
     
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  void initDynamicTrapsParams(const std::string& matName,
			      const Teuchos::ParameterList& trapsPL);
 
}; // end of class NeumannBC_DynamicTraps


}

#endif
