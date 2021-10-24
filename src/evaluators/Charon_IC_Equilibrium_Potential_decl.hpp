
#ifndef CHARON_IC_EQUILIBRIUM_POTENTIAL_DECL_HPP
#define CHARON_IC_EQUILIBRIUM_POTENTIAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Panzer_STK_GatherFields.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::BASIS;
using std::vector;
using std::string;

namespace charon {

/**
 * @brief Evaluate the equilibrium carrier potential
 * by imposing charge neutrality and equilibrium
 * (n - p = C = Na - Nd  &  Efn = Efp)
 */

template<typename EvalT, typename Traits>
class IC_Equilibrium_Potential
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IC_Equilibrium_Potential(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,BASIS> potential; // scaled

  // input
  PHX::MDField<const ScalarT,Cell,BASIS> elec_effdos;  // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> hole_effdos;  // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> acceptor;     // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> donor;        // scaled
  PHX::MDField<const ScalarT,Cell,BASIS> eff_affinity; // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> eff_bandgap;  // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> ref_energy;   // [eV]
  PHX::MDField<const ScalarT,Cell,BASIS> latt_temp;    // scaled

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; // Temperature Scaling in [K]
  double V0; // potential scaling in [V]

  int num_basis;

  std::string dof_name;

  bool haveSuffix;

  std::string eqnSetType;

  Teuchos::RCP<const panzer_stk::STK_Interface> dev_mesh;
  Teuchos::RCP<vector<string>> sc_cnts;
  Teuchos::RCP<vector<string>> sc_blks;
  Teuchos::RCP<vector<double>> sc_wf;
  Teuchos::RCP<vector<double>> sc_vapp;
  std::map<string,vector<stk::mesh::Entity>> sch_cnt_nodes;

  Teuchos::RCP<const charon::Names>  m_names;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters(bool haveSuffix) const;

}; // end of class IC_Equilibrium_Potential


}

#endif
