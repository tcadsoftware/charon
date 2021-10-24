
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"


charon::Scaling_Parameters::Scaling_Parameters (const double& tempScale, const double& concScale,
                                                const double& diffScale, const bool& userMeshScaling,
                                                const double& meshScale, const bool& doIOScaling)
{

  if (userMeshScaling)
    m_X0 = meshScale;
  else
    m_X0 = 1.0e-4;  // Default coordinate scaling [cm]

  m_C0 = concScale;
  m_D0 = diffScale;
  m_T0 = tempScale;

  // Independent scaling parameters
  scale_params.X0 = m_X0;
  scale_params.C0 = m_C0;
  scale_params.D0 = m_D0;
  scale_params.T0 = m_T0;

  // Obtain physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;  // [eV/K]
  double eleQ = cpc.q;      // [C]
  double eps0 = cpc.eps0;   // [C/(V.cm)]

  m_V0 = kbBoltz * m_T0 / 1.0;
  m_Mu0 = m_D0 / m_V0;
  m_R0 = m_D0 * m_C0 / pow(m_X0,2.0);
  m_t0 = m_X0 * m_X0 / m_D0;
  m_E0 = m_V0 / m_X0;
  m_J0 = eleQ * m_D0 * m_C0 / m_X0;

  // Suzey: please do not remove the commented Options here !
  // Option 1
//  m_H0 = m_J0 * m_E0;
//  m_cL0 = m_H0 * m_t0 / m_T0;
//  m_kL0 = m_H0 * m_X0 * m_X0 / m_T0;

  // Option 2
//  m_cL0 = 1.0;
//  m_H0 = m_cL0 * m_T0 / m_t0;
//  m_kL0 = m_H0 * m_X0 * m_X0 / m_T0;

  // Option 3
  m_kL0 = 1.0;
  m_H0 = m_kL0 * m_T0 / (m_X0 * m_X0);
  m_cL0 = m_H0 * m_t0 / m_T0;

  // Option 4
//  m_H0 = 3e13;
//  m_cL0 = m_t0 * m_H0 / m_T0;
//  m_kL0 = m_X0 * m_X0 * m_H0 / m_T0;

  // Option 5
//  m_H0 = 3e13;
//  m_cL0 = 2.34;
//  m_kL0 = 1/1.74;

  // Dependent scaling parameters
  scale_params.V0 = m_V0;
  scale_params.Mu0 = m_Mu0;
  scale_params.R0 = m_R0;
  scale_params.t0 = m_t0;
  scale_params.E0 = m_E0;
  scale_params.J0 = m_J0;
  scale_params.H0 = m_H0;
  scale_params.cL0 = m_cL0;
  scale_params.kL0 = m_kL0;

  // Laplacian scaling [unitless]
  scale_params.Lambda2 = m_V0 * eps0 / (eleQ * pow(m_X0,2.0) * m_C0);

  scale_params.cLPrefactor = m_cL0 * m_T0 / (m_t0 * m_H0);
  scale_params.kLPrefactor = m_kL0 * m_T0 / (m_X0 * m_X0 * m_H0);

  // These are used to unscale variables when output to the Exodus
  // file. You manually add names to the map here and they will be
  // unscaled in the Exodus output. If a variable isn't included in this
  // map it will be output scaled to the Exodus file regardless of the
  // input file setting. On the other hand if you add names here then
  // they should be output to the Exodus file unscaled. Note that this
  // also applies to the initial guess if it's read in from an Exodus
  // file. The guess will be scaled right after being read. And scaling
  // the input initial guess automatically turns on unscaling of output
  // variables to the Exodus file.
  varScaleFactors = Teuchos::rcp(new std::map<std::string,double>);
  if (doIOScaling)
  {
    (*varScaleFactors)["ELECTRIC_POTENTIAL"] = m_V0;
    (*varScaleFactors)["ELECTRON_DENSITY"] = m_C0;
    (*varScaleFactors)["HOLE_DENSITY"] = m_C0;

    (*varScaleFactors)["Electron Mobility"] = m_Mu0;
    (*varScaleFactors)["Hole Mobility"] = m_Mu0;

    (*varScaleFactors)["Electron Diffusion Coefficient"] = m_D0;
    (*varScaleFactors)["Hole Diffusion Coefficient"] = m_D0;

    (*varScaleFactors)["GRAD_ELECTRIC_POTENTIALX"] = m_E0;
    (*varScaleFactors)["GRAD_ELECTRIC_POTENTIALY"] = m_E0;
    (*varScaleFactors)["GRAD_ELECTRIC_POTENTIALZ"] = m_E0;
  }

  const int me = Teuchos::GlobalMPISession::getRank();

  if(me==0)
    {
      std::cout << std::endl;
      std::cout << "Coordinate Scaling X0 = " << scale_params.X0 << " [cm]" << std::endl;
      std::cout << "Concentration Scaling C0 = " << scale_params.C0 << " [1/cm^3]" << std::endl;
      std::cout << "Diffusion Coefficient Scaling D0 = " << scale_params.D0 << " [cm^2/s]" << std::endl;
      std::cout << "Temperature Scaling T0 = " << scale_params.T0 << " [K]" << std::endl;
      std::cout << "Voltage Scaling V0 = " << scale_params.V0 << " [V]" << std::endl;
      std::cout << "Mobility Scaling Mu0 = " << scale_params.Mu0 << " [cm^2/(V.s)]" << std::endl;
      std::cout << "Source Scaling R0 = " << scale_params.R0 << " [1/(cm^3.s)]" << std::endl;
      std::cout << "Time Scaling t0 = " << scale_params.t0 << " [s]" << std::endl;
      std::cout << "Electric Field Scaling E0 = " << scale_params.E0 << " [V/cm]" << std::endl;
      std::cout << "Current Density Scaling J0 = " << scale_params.J0 << " [A/cm^2]" << std::endl;
      std::cout << "Heat Generation Scaling H0 = " << scale_params.H0 << " [W/cm^3]" << std::endl;
      std::cout << "Heat Capacity Scaling cL0 = " << scale_params.cL0 << " [J/(cm^3.K)]" << std::endl;
      std::cout << "Thermal Conductivity Scaling kL0 = " << scale_params.kL0 << " [W/(K.cm)]" << std::endl;
      std::cout << "Laplacian Scaling Lambda2 = " << scale_params.Lambda2 << " [unitless]" << std::endl;
      std::cout << std::endl;
    }
}
