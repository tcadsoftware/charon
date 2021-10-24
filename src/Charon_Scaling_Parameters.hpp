
#ifndef CHARON_SCALING_PARAMETERS_HPP
#define CHARON_SCALING_PARAMETERS_HPP

#include <string>
#include <vector>
#include "Teuchos_ParameterList.hpp"

namespace charon
{

// set up the scaling parameters, instantiate the class in Charon_Main.cpp
// so that the parameters are accessible everywhere in the program

  class Scaling_Parameters
  {
  public:

    // ctor to compute scaling parameters
    Scaling_Parameters (const double& tempScale, const double& concScale, const double& diffScale,
                        const bool& userMeshScaling, const double& meshScaling,
                        const bool& doOutputUnscaling);

    struct ScaleParams
    {
      /* Independent scaling parameters */

      // Coordinate Scaling [cm]
      double X0;

      // Concentration Scaling [cm^-3]
      double C0;

      // Diffusion Coefficient Scaling [cm^2/s]
      double D0;

      // Temperature Scaling [K]
      double T0;

      /* Dependent scaling parameters */

      // Voltage Scaling [V] or Energy Scaling [eV]
      double V0;

      // Mobility Scaling [cm^2/(V.s)]
      double Mu0;

      // Source Scaling [cm^-3.s^-1]
      double R0;

      // Time Scaling [s]
      double t0;

      // Electric Field Scaling [V/cm]
      double E0;

      // Current Density Scaling [A/cm^2]
      double J0;

      // Heat Generation Scaling [W/cm^3]
      double H0;

      // Heat Capacity Scaling [J/(cm^3.K)]
      double cL0;

      // Thermal Conductivity Scaling [W/(K.cm)]
      double kL0;

      // Laplacian Scaling [1]
      double Lambda2;

      // Prefactors [unitless] for the time-derivative and Laplacian terms of the heat equation
      double cLPrefactor;  // for the time-derivative term
      double kLPrefactor;  // for the Laplacian term

    };

    ScaleParams scale_params;

    Teuchos::RCP<std::map<std::string,double> > varScaleFactors;

  private:
    double m_X0;
    double m_C0;
    double m_D0;
    double m_T0;
    double m_V0;
    double m_Mu0;
    double m_R0;
    double m_t0;
    double m_E0;
    double m_J0;
    double m_H0;
    double m_cL0;
    double m_kL0;

  };

} // END "namespace charon"

#endif
