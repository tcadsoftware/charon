
#ifndef CHARON_DDLATTICE_CURRENTDENSITY_DECL_HPP
#define CHARON_DDLATTICE_CURRENTDENSITY_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;


namespace charon {

template<typename EvalT, typename Traits>
class DDLattice_CurrentDensity
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLattice_CurrentDensity(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,IP,Dim> current_density;
  PHX::MDField<ScalarT,Cell,IP> ion_curr_dens_x;
  PHX::MDField<ScalarT,Cell,IP> ion_curr_dens_y;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_latt_temp; // gradient of lattice temperature
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_carr_dens; // gradient of carrier density
  PHX::MDField<const ScalarT,Cell,IP,Dim> electric_field; // electric field
  PHX::MDField<const ScalarT,Cell,IP,Dim> velocity;       // ion velocity

  PHX::MDField<const ScalarT,Cell,IP> carr_dens;
  PHX::MDField<const ScalarT,Cell,IP> diff_coeff;
  PHX::MDField<const ScalarT,Cell,IP> mobility;
  PHX::MDField<const ScalarT,Cell,IP> thermodiff_coeff;

  // IP
  std::size_t num_ips;
  std::size_t num_dims;

  // carrier type
  std::string carrType;

  // a sign factor to differ electrons from holes
  double sign;

  bool bTempGrad;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DDLattice_CurrentDensity


}

#endif
