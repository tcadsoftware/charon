

#ifndef CHARON_BAND2BAND_TUNNELING_LOCAL_DECL_HPP
#define CHARON_BAND2BAND_TUNNELING_LOCAL_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::BASIS;
using panzer::Dim;

namespace charon {

//! obtain band2band tunneling rate
template<typename EvalT, typename Traits>
class Band2Band_Tunneling_Local
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Band2Band_Tunneling_Local(
      const Teuchos::ParameterList& p);

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

  // output
  PHX::MDField<ScalarT,Cell,Point> bbt_rate;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_drForce;
  PHX::MDField<const ScalarT,Cell,Point,Dim> hole_drForce;
  PHX::MDField<const ScalarT,Cell,Point,Dim> initial_grad_phi;

  PHX::MDField<const ScalarT, Cell, Point, Dim> curr_dens_e;
  PHX::MDField<const ScalarT, Cell, Point, Dim> curr_dens_h;

  PHX::MDField<const ScalarT, Cell, Point> dens_e;
  PHX::MDField<const ScalarT, Cell, Point> dens_h;
  PHX::MDField<const ScalarT, Cell, Point> intrin_conc;
  PHX::MDField<const ScalarT, Cell, Point> eff_bandgap;
  PHX::MDField<const ScalarT, Cell, Point> net_doping;
  PHX::MDField<const ScalarT, Cell, Point> rel_perm;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double R0;      // recomb./gen. scaling
  double E0;      // electric field scaling
  double C0;      // concentration scaling

  int num_points;
  int num_dims;
  int num_nodes;

  std::size_t basis_index;
  std::string basis_name;

  int int_rule_degree;
  std::size_t int_rule_index;
 
  bool isSGCVFEM;
  bool bAddFactor = false; 

  // Model name
  std::string bbtModel;

  // Kane model parameters
  double A_Kane, B_Kane, gamma_Kane, alpha_Kane, beta_Kane; 

  // Hurkx model parameters
  double A_Hurkx, B_Hurkx, gamma_Hurkx, alpha_Hurkx, beta_Hurkx; 

  // Schenk model parameters
  double A_Schenk, B_Schenk, HW_Schenk;

  // driving force name
  std::string driveForce;

  // Minimum value of the electric field
  double minField;

  double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax; 

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  /**
   * @brief Initialize Band2Band tunneling parameters
   */
  void initBBTParams(const std::string& matName, const Teuchos::ParameterList& bbtParamList);

}; // end of class Band2Band_Tunneling_Local


}

#endif
