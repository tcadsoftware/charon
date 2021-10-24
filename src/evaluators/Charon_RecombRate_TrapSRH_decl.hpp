
#ifndef CHARON_RECOMBRATE_TRAPSRH_DECL_HPP
#define CHARON_RECOMBRATE_TRAPSRH_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

/**
 * @brief This evaluator uses the SRH recombination expression to model the effects
 * of generic traps. The model neglects the time dynamics of traps. In another word,
 * the dynamics of traps is so fast that they 'instantaneously' reach equilibrium when
 * compared to the carrier dynamics. Note that modeling the effect of traps using
 * the SRH recombination directly implies that the net electron and hole capture rates
 * via the traps are equal. However, the same capture rates do not mean the net charge
 * from the recombination process is zero, so still need to add trap charge to
 * the Poisson equation.
 *
 * The type of traps are characterized by Et (trap energy measured from the conduction
 * band edge), Nt (trap density), S (Huang-Rhys factor), hbar_omega0 (phonon energy), and
 * tau0 (field-independent lifetime) or sigma (field-independent capture cross section).
 * The net electron-hole recombination due to one type of traps is given by
 * R = (n*p-ni^2)/(taup*(n+n1) + taun*(p+p1)), where tau = tau0 / (1 + g(F)), tau0
 * is given by user, or computed from 1 / (sigma*vth*Nt), g(F) is a field-enhancement
 * factor which models the effect of band to trap tunneling. Currently, g(F) can be either
 * set to 0 or computed using the Schenk's model, i.e.,
 * A. Schenk, "A model for the field and temperature dependence of Shockley-Read-Hall
 * lifetimes in silicon," Solid-State Electronics, Vol. 35, No. 11, pp. 1585-1596 (1992).
 *
 * For multiple type of traps, the code sums up the SRH recombination rates. More details
 * on the model are found in the Charon user manual.
 */


template<typename EvalT, typename Traits>
class RecombRate_TrapSRH
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    RecombRate_TrapSRH(
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
  PHX::MDField<ScalarT,Cell,Point> trap_srh_rate;
  PHX::MDField<ScalarT,Cell,Point> trap_srh_charge;
  PHX::MDField<ScalarT,Cell,Point> trap_srh_deriv_e;
  PHX::MDField<ScalarT,Cell,Point> trap_srh_deriv_h;

  // input
  PHX::MDField<const ScalarT,Cell,Point> edensity;
  PHX::MDField<const ScalarT,Cell,Point> hdensity;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point,Dim> elec_field;
  PHX::MDField<const ScalarT,Cell,Point,Dim> hole_field;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double t0; // time scaling, [s]
  double T0; // temperature scaling, [K]
  double E0; // electric field scaling, [V/cm]
  double X0; // coordinate scaling, [cm]
  double C0; // concentration scaling, [cm^(-3)]

  int num_points, num_dims;
  int int_rule_degree;
  int num_nodes;
  std::size_t int_rule_index;
  std::size_t basis_index;
  std::string basis_name; 
  std::string driveForce;     
  bool isSGCVFEM;

  std::vector<double> energyLevel, trapDensity;
  std::vector<double> minXCoord, maxXCoord, minYCoord, maxYCoord, minZCoord, maxZCoord;
  std::vector<double> phononEnergy, hRhysFactor;
  std::vector<double> eLifeTime, hLifeTime, eEffMass, hEffMass;
  std::vector<double> eHJLocation, hHJLocation, eHJBandOffset, hHJBandOffset;
  std::vector<bool> eTimeGiven, hTimeGiven;
  std::vector<std::string> eTunnelModel, hTunnelModel, eTunnelDir, hTunnelDir;
  std::vector<std::string> trapSpcProfile, trapType;

  double hbar, q, m0, pi, kb;  // universal constants
  double fieldSign; 
  
  // The following parameters are used for adaptive trapzoidal integration 
  double field4Adaptive, kbT4Adaptive, bandGap4Adaptive, spcLoc4Adaptive; 
  int itrap4Adaptive; 
  std::string carrType4Adaptive; 
  bool isNewDOSNum = false; 

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  /**
   * @brief Initialize trap parameters
   */
  void initTrapSRHParams(const std::string& matName, const Teuchos::ParameterList& trapSRHParamList);

  /**
   * @brief Evaluate the Schenk field enhancement factor for Tunneling Model = Schenk HighTemp,
   * Schenk LowTemp, Schenk AsymConstFDOS, or Schenk ConstFDOS
   */
  double evalSchenkFieldFactor(const ScalarT& fieldAmp, const ScalarT& kbT, const ScalarT& bandGap,
                               const int& itrap, const std::string& carrType);

  /**
   * @brief Evaluate the Schenk field enhancement factor for Tunneling Model = Schenk HighTemp
   * [Eq.(40)] or Schenk LowTemp [Eq.(49) in Schenk's paper], called by evalSchenkFieldFactor(.)
   */
  double schenkTemperatureApprox(const double& hbarTheta, const double& dkbT, const double& Et,
                                 const int& itrap, const std::string& tunModel);

  /**
   * @brief Evaluate the denominator of the Schenk field enhancement factor, and
   * the denominator remains the same for all the Schenk models that use the energy integration
   */
  double schenkFieldFactorDenominator(const double& dkbT, const double& Et, const double& Eph, const double& z);

  /**
   * @brief Evaluate the numerator of the Schenk field enhancement factor for
   * Tunneling Model = Schenk AsymConstFDOS and Schenk ConstFDOS
   */
  double schenkFieldFactorNumerator(const double& hbarTheta, const double& dkbT, const double& Et,
                                    const double& Eph, const double& z, const std::string& tunModel);

  /**
   * @brief Evaluate the Schenk field enhancement factor for Tunneling Model = Schenk NewDOS
   * using the new DOS model developed by Suzey Gao
   */
  double evalFieldFactorWithNewDOS(const ScalarT& fieldAmp, const ScalarT& kbT, const ScalarT& bandGap,
                                   const int& itrap, const std::string& carrType, const double& spcLoc);

  /**
   * @brief Evaluate the numerator of the Schenk field enhancement factor for Tunneling Model = Schenk NewDOS
   * using the uniform trapezoidal rule
   */
  double fieldFactorWithNewDOSNumerator(const ScalarT& fieldAmp, const ScalarT& kbT, const ScalarT& bandGap,
                                        const int& itrap, const std::string& carrType, const double& spcLoc);

  /**
   * @brief Calculate the tunneling DOS for Tunneling Model = Schenk NewDOS
   */
  double calcTunnelDOSForSchenkNewModel(const double& hbarTheta, const double& bOffset, const double& Et, const double& Eloc,
                                        const double& E, const double& dEx, const double& xloc, const double& alpha, const double& me);

  /**
   * @brief Calculate the tunneling DOS for a step barrier (developed by Suzey Gao) using uniform trapezoidal rule
   */
  double calcDOSForStepBarrier(const double& Enew, const double& dEx, const double& bOffset, const double& xloc, const double& me);

  /**
   * @brief Calculate the tunneling DOS for a linear potential with band offset (developed by Suzey Gao) using uniform trapezoidal rule
   */
  double calcDOSForLinPotWithOffset(const double& hbarTheta, const double& bOffset, const double& Enew,
                                    const double& dEx, const double& xloc, const double& alpha);

  /**
   * @brief Calculate the tunneling DOS for a step barrier using the Gauss-Jacobi quadrature rule
   */
  double calcDOSForStepBarrierGaussQR(const double& Enew, const double& bOffset, const double& xloc, const double& me);

  /**
   * @brief Calculate the tunneling DOS for a linear potential with band offset using the Gauss-Jacobi quadrature rule
   */
  double calcDOSForLinPotWithOffsetGaussQR(const double& hbarTheta, const double& bOffset, const double& Enew,
                                    const double& xloc, const double& alpha);

  /**
   * @brief Use adaptive trapezoidal rule for integration
   */
  double adaptiveIntegrate(const double a, const double b, const double tol);

  /**
   * @brief a subfunction for adaptiveIntegrate(...)
   */
  double adaptlobstp(const double a, const double b, const double fa, const double fb, const double is);

  /**
   * @brief Compute the integrand for adaptive trapezoidal integration
   */
  double fieldFactorIntegrand(const double E);


}; // end of class RecombRate_TrapSRH


}

#endif
