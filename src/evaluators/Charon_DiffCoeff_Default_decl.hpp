
#ifndef CHARON_DIFFCOEFF_DEFAULT_DECL_HPP
#define CHARON_DIFFCOEFF_DEFAULT_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"


using panzer::Cell;
using panzer::Point;

namespace charon {

//! obtain default electron and hole diffusion coefficient
template<typename EvalT, typename Traits>
class DiffCoeff_Default
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    DiffCoeff_Default(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> diffcoeff; // scaled

  // input (all fields are scaled, except T0)
  PHX::MDField<const ScalarT,Cell,Point> mobility;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;

  PHX::MDField<const ScalarT,Cell,Point> carr_dens;
  PHX::MDField<const ScalarT,Cell,Point> eff_dos;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0;  // temp. scaling in [K]

  int num_points;
  int num_edges;

  bool isEdgedl;
  bool bUseFD;

  std::string carrType;
  std::string fdFormula;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;
  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > forwardFermiIntegral;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

}; // end of class DiffCoeff_Default


}

#endif
