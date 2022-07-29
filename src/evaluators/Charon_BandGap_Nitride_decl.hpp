#ifndef CHARON_BANDGAP_NITRIDE_DECL_HPP
#define CHARON_BANDGAP_NITRIDE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {
/**
 *  \brief Computes Bandgap for nitrides with respect to mole fraction
 *         R.Passler Phys.Stat.Solidi(b) 216 (1999):975
 */
    
//! calculate temperature-dependent band gap
PHX_EVALUATOR_CLASS(BandGap_Nitride)

private:

  // output
  PHX::MDField<ScalarT,Cell,Point> band_gap; 

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> molefrac;

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double T0; 
  
  int num_points;
  
  // bandgap model parameters

  std::string materialName;
  bool isBinary, isTernary;

    /**
     *  \brief Calculates the bandgap of binary nitrides
     *
     *  \param[in]  arg1  temperature
     * 
     *  \param[in]  arg2  material name
     *
     *  \returns A double of the bandgap energy
     */
  auto binaryBandgap(const ScalarT& lattT, const std::string& materialName) ->ScalarT;
    /**
     *  \brief Calculates the bandgap of ternary nitrides
     *
     *  \param[in]  arg1  temperature
     * 
     *  \param[in]  arg2  material name
     * 
     *  \param[in]  arg3  mole fraction
     *
     *  \returns A double of the bandgap energy
     */
  auto ternaryBandgap(const ScalarT& lattT, const std::string& materialName,
                      const ScalarT& mole_frac) -> ScalarT;
  

private:
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
PHX_EVALUATOR_CLASS_END

}

#endif
