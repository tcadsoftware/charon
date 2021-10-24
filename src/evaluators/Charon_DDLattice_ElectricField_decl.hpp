
#ifndef CHARON_DDLATTICE_ELECTRICFIELD_DECL_HPP
#define CHARON_DDLATTICE_ELECTRICFIELD_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::IP;
using panzer::Dim;
using panzer::Point;


namespace charon {

template<typename EvalT, typename Traits>
class DDLattice_ElectricField
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DDLattice_ElectricField(
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
  PHX::MDField<ScalarT,Cell,IP,Dim> electric_field;

  // input
  PHX::MDField<const ScalarT,Cell,IP,Dim> grad_potential; // gradient of potential \phi

  PHX::MDField<const ScalarT,Cell,Point> bandgap;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double V0; // [V]

  // IP
  std::size_t num_cells;
  std::size_t num_ips;
  std::size_t num_dims;

  // BASIS
  std::string basis_name;
  std::size_t basis_index;
  std::size_t num_basis;

  // carrier type
  std::string carrType;

  // model name
  std::string modelName;

  bool haveBGN;
  double sign;

  // intermediate fields
  Kokkos::DynRankView<ScalarT,PHX::Device> dEg;
  Kokkos::DynRankView<ScalarT,PHX::Device> grad_dEg;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class DDLattice_ElectricField


}

#endif
