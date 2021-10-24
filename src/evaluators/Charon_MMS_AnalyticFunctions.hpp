
#ifndef CHARON_MMS_ANALYTICFUNCTIONS
#define CHARON_MMS_ANALYTICFUNCTIONS

#include <vector>
#include "Charon_config.hpp"
#include "Charon_Names.hpp"

#include "Charon_Physical_Constants.hpp"
#include "Charon_Scaling_Parameters.hpp"

namespace charon {

/**
 * @brief Class for MMS functions
 *
 * This class contains the analytic functions for MMS problems used in
 * Charon.
 */
class MMS_NLP_GLH_1_AnalyticFunction
{
public:

  double net_doping(double const x)
  {

    // Mathematica function yielding the doping value. The value of the
    // input coordinate needs to be in centimeters
    double value =
      -1.3680544484358358e15/std::exp(11.603171533552503*std::erfc(-5. + 20000.*x)) +
      114213.2903983817*std::exp(11.603171533552503*std::erfc(-5. + 20000.*x)) +
      (8.830264295490818e15 - 3.5321057181963272e19*x)/std::exp(4.e8*(-0.00025 + x)*(-0.00025 + x));

    return value;
  }
};

class MMS_DD_RDH_AnalyticFunction
{
public:

  virtual ~MMS_DD_RDH_AnalyticFunction(){}

  MMS_DD_RDH_AnalyticFunction():
    vl(0.4),
    beta(1e6),
    length(5e-5)
  {
  }

  virtual double edensity(double const& x) const = 0;
  virtual double hdensity(double const& x) const = 0;

  /**
   * @brief The analytic function for electric potential
   *
   * This returns the analytic solution for the electric potential. The
   * value is unscaled. The analytic function is:
   *
   * - @f$\psi(x) = \text{-VL}\ * (\text{atan}(\text{BETA}(x - \text{LENGTH}
   *   * 0.5))\over atan(\text{BETA} * \text{LENGTH} * 0.5))@f$
   *
   */
  double potential(double const& x) const;


  /**
   * @brief The analytic function for the net doping profile
   *
   * This function returns the unscaled value of the net
   * doping profile. The analytic function is:
   *
   *
   * -@f$C(x) = (n(x) - p(x)) + \frac{(2 * \text{BETA}^3 * \epsilon * \text{VL} *
   *      (-\text{LENGTH} * 0.5 + x))}
   *   {(q * (1 + \text{BETA}^2 * (-\text{LENGTH}*0.5 + x)^2)^2 * atan((\text{BETA}*\text{LENGTH})*0.5))}@f$
   *
   *
   */
  std::vector<double> doping(double const& x) const;


protected:
  double dop(double const& x) const;
  double vl;
  double beta;
  double length;
  double lambda;
  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;

};

class MMS_DD_RDH_1_AnalyticFunction:
  public MMS_DD_RDH_AnalyticFunction

{

public:
  MMS_DD_RDH_1_AnalyticFunction(Teuchos::ParameterList const& p):
    nhat(3e3),
    phat(2e3)
  {
    scaleParams = p.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");
    lambda = scaleParams->scale_params.V0;
  }

  /**
   * @brief The analytic function for the electron concentration
   *
   * This function returns the unscaled value of the electron
   * density. The analytic function is:
   *
   * -@f$edensity(x) = \text{NHAT}*e^{\frac{(\psi(x)} {\text{LAMBDA})}}@f$
   *
   */
  double edensity(double const& x) const;

  /**
   * @brief The analytic function for the hole concentration
   *
   * This function returns the unscaled value of the hole
   * density. The analytic function is:
   *
   * -@f$hdensity(x) = \text{PHAT}*e^{-\frac{(\psi(x)} {\text{LAMBDA})}}@f$
   *
   */
  double hdensity(double const& x) const;

private:
  double nhat;
  double phat;
};


class MMS_DD_RDH_2_AnalyticFunction:
  public MMS_DD_RDH_AnalyticFunction

{

public:
  MMS_DD_RDH_2_AnalyticFunction(Teuchos::ParameterList const& p):
    N0(4.6e17)
  {
    scaleParams = p.get< Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");
    lambda = scaleParams->scale_params.V0;
  }

  /**
   * @brief The analytic function for the electron concentration
   *
   * This function returns the unscaled value of the electron
   * density. The analytic function is:
   *
   * -@f$edensity(x) = 6.39145e10*e^(-10.0887*atan(1.0e6*(-0.000025+x)))@f$
   *
   */
  double edensity(double const& x) const;

  /**
   * @brief The analytic function for the hole concentration
   *
   * This function returns the unscaled value of the hole
   * density. The analytic function is:
   *
   * -@f$hdensity(x) = 3.22851e10*e^(-10.0887*atan(1.0e6*(-0.000025+x)))@f$
   *
   */
  double hdensity(double const& x) const;

  /**
   * @brief The analytic function for the donor concentration
   *
   * This function returns the unscaled value of the acceptor and donor
   * concentration. The analytic function is:
   *
   * -@f$Na(x) = \frac{1}{2} * (N + C)@f$
   */
  double Na(double const& x) const;

  /**
   * @brief The analytic function for the acceptor concentration
   *
   * This function returns the unscaled value of the acceptor and donor
   * concentration. The analytic function is:
   *
   * -@f$Nd(x) = \frac{1}{2} * (N - C)@f$
   */
  double Nd(double const& x) const;

private:
  /**
   * @brief The analytic function for the sum of the donor and acceptor
   * concentrations
   *
   * This function returns the unscaled value of the acceptor and donor
   * concentration. The analytic function is:
   *
   * -@f$N(x) = \text{N0}*e^{\frac{-4x} {\text{LENGTH}}}@f$
   */
  double bigN(double const& x) const;
  double N0;
};

}
#endif
