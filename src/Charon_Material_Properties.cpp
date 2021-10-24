
#include "Charon_Material_Properties.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"


bool charon::Material_Properties::mFirst = true;


// **********************************************************************
// private constructor to initialize default values
// **********************************************************************
charon::Material_Properties::Material_Properties ()
{
  // All the setup should be done only once.
  if (mFirst)
  {
    // Create the Materials parameter list that holds double values for
    // properties.
    {
      pMaterials.setName("Materials");

      // Set the material-independent properties.
      pMaterials.set<double>("Lattice Temperature", 300.0, "[K]");
      pMaterials.set<double>("Concentration Scaling", 1.25e10, "[cm^-3]");

      // Silicon
      Teuchos::ParameterList& pSilicon    = pMaterials.sublist("Silicon",
        false, "Sublist defining properties for Silicon");
      setSiliconParameters(pSilicon);

      // SiO2
      Teuchos::ParameterList& pSiO2       = pMaterials.sublist("SiO2",
        false, "Sublist defining properties for SiO2.");
      setSiO2Parameters(pSiO2);

      // OxyNitride
      Teuchos::ParameterList& pOxyNitride = pMaterials.sublist("OxyNitride",
        false, "Sublist defining properties for OxyNitride.");
      setOxyNitrideParameters(pOxyNitride);

      // GaAs
      Teuchos::ParameterList& pGaAs       = pMaterials.sublist("GaAs",
        false, "Sublist defining properties for GaAs.");
      setGaAsParameters(pGaAs);

      // InP
      Teuchos::ParameterList& pInP        = pMaterials.sublist("InP",
        false, "Sublist defining properties for InP.");
      setInPParameters(pInP);

      // Al(x)Ga(1-x)As
      Teuchos::ParameterList& pAlGaAs     = pMaterials.sublist("AlGaAs",
        false, "Sublist defining properties for AlGaAs.");
      setAlGaAsParameters(pAlGaAs);

      // In(1-x)Ga(x)As
      Teuchos::ParameterList& pInGaAs     = pMaterials.sublist("InGaAs",
        false, "Sublist defining properties for InGaAs.");
      setInGaAsParameters(pInGaAs);

      // Al(x)In(1-x)As
      Teuchos::ParameterList& pAlInAs     = pMaterials.sublist("AlInAs",
        false, "Sublist defining properties for AlInAs.");
      setAlInAsParameters(pAlInAs);

      // GaAs(x)P(1-x)
      Teuchos::ParameterList& pGaAsP      = pMaterials.sublist("GaAsP",
        false, "Sublist defining properties for GaAsP.");
      setGaAsPParameters(pGaAsP);

      // In(1-x)Ga(x)P
      Teuchos::ParameterList& pInGaP      = pMaterials.sublist("InGaP",
        false, "Sublist defining properties for InGaP.");
      setInGaPParameters(pInGaP);

      // Al(1-x)Ga(x)N
      Teuchos::ParameterList& pAlGaN      = pMaterials.sublist("AlGaN",
        false, "Sublist defining properties for AlGaN.");
      setAlGaNParameters(pAlGaN);

      // AlN
      Teuchos::ParameterList& pAlN        = pMaterials.sublist("AlN",
        false, "Sublist defining properties for AlN.");
      setAlNParameters(pAlN);

      // GaN
      Teuchos::ParameterList& pGaN        = pMaterials.sublist("GaN",
        false, "Sublist defining properties for GaN.");
      setGaNParameters(pGaN);

      // TiO2
      Teuchos::ParameterList& pTiO2       = pMaterials.sublist("TiO2",
        false, "Sublist defining properties for TiO2.");
      setTiO2Parameters(pTiO2);

      // Tantalum
      Teuchos::ParameterList& pTa         = pMaterials.sublist("Tantalum",
        false, "Sublist defining properties for Tantalum.");
      setTantalumParameters(pTa);

      // Platinum
      Teuchos::ParameterList& pPt         = pMaterials.sublist("Platinum",
        false, "Sublist defining properties for Platinum.");
      setPlatinumParameters(pPt);

      // PlatinumSemi
      Teuchos::ParameterList& pPtSC       = pMaterials.sublist("PlatinumSemi",
        false, "Sublist defining properties for PlatinumSemi.");
      setPlatinumSemiParameters(pPtSC);

      // Ta2O5
      Teuchos::ParameterList& pTa2O5      = pMaterials.sublist("Ta2O5",
        false, "Sublist defining properties for Ta2O5.");
      setTa2O5Parameters(pTa2O5);

      // Ta(O)
      Teuchos::ParameterList& pTaO        = pMaterials.sublist("Ta(O)",
        false, "Sublist defining properties for Ta(O).");
      setTaOParameters(pTaO);
    }

    // reset mFirst
    mFirst = false;
  } // end if(mFirst)
}

///////////////////////////////////////////////////////////////////////////////
//
//  setSiliconParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setSiliconParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          11.8,    "[1]"             );  // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.17,    "[eV]"            );  // otherwise noted.
  p.set<double>("Band Gap",                       1.08,    "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.25e10, "[cm^-3]"         );
  p.set<double>("Electron Effective Mass",        0.32,    "[1]"             );
  p.set<double>("Hole Effective Mass",            0.16,    "mlh:[1]"         );
  p.set<double>("Electron Mobility",              1000.,   "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  500.,    "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 25.8,    "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     12.5,    "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-7,    "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-7,    "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,    "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,    "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,    "TPowerLaw:[1]"   );  // From sdevice.
  p.set<double>("Hole SRH TPowerLaw",             -1.5,    "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,    "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,    "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,       "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    0.0,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    2.8e-31,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    9.9e-32,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.17,    "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.08,    "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             4.73e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              636.0,   "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 2.8e19,  "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,     "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1.04e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,     "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  9.e-3, "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17,  "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.5,   "CON_BGN:[1]"   );

  // Slotboom
  p.set<double>("Slotboom BGN V0_BGN",  6.92e-3, "V0_BGN:[V]"    );
  p.set<double>("Slotboom BGN N0_BGN",  1.3e17,  "N0_BGN:[cm^-3]");
  p.set<double>("Slotboom BGN CON_BGN", 0.5,     "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -14.84e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", 0.78e-3,   "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 15.08e-3,  "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 0.74e-3,   "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", 0.0,       "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", -16.27e-3, "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.18e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 18.46e-3,  "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 0.0,       "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", -2.63e-3,  "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 1429.23,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 55.24,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  1.072e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", -2.2,     "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   -3.8,     "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 0.73,     "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 479.37,   "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 49.705,   "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  1.606e17, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", -2.2,     "[1]"         );
  p.set<double>("Analytic Hole xin",   -3.7,     "[1]"         );
  p.set<double>("Analytic Hole alpha", 0.70,     "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    1252.0,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    88.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1.26e17, "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 0.88,    "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     -0.57,   "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     -2.33,   "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     2.4,     "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     -0.146,  "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",    407.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",    54.3,    "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",     2.35e17, "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp", 0.88,    "[1]"         );
  p.set<double>("Arora Hole exp1",     -0.57,   "[1]"         );
  p.set<double>("Arora Hole exp2",     -2.23,   "[1]"         );
  p.set<double>("Arora Hole exp3",     2.4,     "[1]"         );
  p.set<double>("Arora Hole exp4",     -0.146,  "[1]"         );

  //---------------------------------------------------------------------------
  // University of Bologna mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons (As)
  p.set<double>("UniBo Electron(As) mumax",    1441,    "mumax:[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(As) c",        0.07,    "c:[1]"             );
  p.set<double>("UniBo Electron(As) gamma",    -2.45,   "gamma:[1]"         );
  p.set<double>("UniBo Electron(As) mu0d",     55.0,    "mu0d:[cm^2/(V.s)]" );
  p.set<double>("UniBo Electron(As) mu0d_exp", -0.6,    "mu0d_exp:[1]"      );
  p.set<double>("UniBo Electron(As) mu0a",     132.0,   "mu0a:[cm^2/(V.s)]" );
  p.set<double>("UniBo Electron(As) mu0a_exp", -1.3,    "mu0a_exp:[1]"      );
  p.set<double>("UniBo Electron(As) mu1d",     42.4,    "mu1d:[cm^2/(V.s)]" );
  p.set<double>("UniBo Electron(As) mu1d_exp", -0.5,    "mu1d_exp:[1]"      );
  p.set<double>("UniBo Electron(As) mu1a",     73.5,    "mu1a:[cm^2/(V.s)]" );
  p.set<double>("UniBo Electron(As) mu1a_exp", -1.25,   "mu1a_exp:[1]"      );
  p.set<double>("UniBo Electron(As) cr1",      8.9e16,  "cr1:[cm^-3]"       );
  p.set<double>("UniBo Electron(As) cr1_exp",  3.65,    "cr1_exp:[1]"       );
  p.set<double>("UniBo Electron(As) cr2",      1.22e17, "cr2:[cm^-3]"       );
  p.set<double>("UniBo Electron(As) cr2_exp",  2.65,    "cr2_exp:[1]"       );
  p.set<double>("UniBo Electron(As) cs1",      2.9e20,  "cs1:[cm^-3]"       );
  p.set<double>("UniBo Electron(As) cs1_exp",  0.0,     "cs1_exp:[1]"       );
  p.set<double>("UniBo Electron(As) cs2",      7.0e20,  "cs2:[cm^-3]"       );
  p.set<double>("UniBo Electron(As) alpha1",   0.68,    "alpha1:[1]"        );
  p.set<double>("UniBo Electron(As) alpha2",   0.72,    "alpha2:[1]"        );

  // Electrons (P)
  p.set<double>("UniBo Electron(P) mumax",    1441,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(P) c",        0.07,    "[1]"         );
  p.set<double>("UniBo Electron(P) gamma",    -2.45,   "[1]"         );
  p.set<double>("UniBo Electron(P) mu0d",     62.2,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(P) mu0d_exp", -0.7,    "[1]"         );
  p.set<double>("UniBo Electron(P) mu0a",     132.0,   "[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(P) mu0a_exp", -1.3,    "[1]"         );
  p.set<double>("UniBo Electron(P) mu1d",     48.6,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(P) mu1d_exp", -0.7,    "[1]"         );
  p.set<double>("UniBo Electron(P) mu1a",     73.5,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Electron(P) mu1a_exp", -1.25,   "[1]"         );
  p.set<double>("UniBo Electron(P) cr1",      8.5e16,  "[cm^-3]"     );
  p.set<double>("UniBo Electron(P) cr1_exp",  3.65,    "[1]"         );
  p.set<double>("UniBo Electron(P) cr2",      1.22e17, "[cm^-3]"     );
  p.set<double>("UniBo Electron(P) cr2_exp",  2.65,    "[1]"         );
  p.set<double>("UniBo Electron(P) cs1",      4.0e20,  "[cm^-3]"     );
  p.set<double>("UniBo Electron(P) cs1_exp",  0.0,     "[1]"         );
  p.set<double>("UniBo Electron(P) cs2",      7.0e20,  "[cm^-3]"     );
  p.set<double>("UniBo Electron(P) alpha1",   0.68,    "[1]"         );
  p.set<double>("UniBo Electron(P) alpha2",   0.72,    "[1]"         );

  // Holes
  p.set<double>("UniBo Hole mumax",    470.5,   "[cm^2/(V.s)]");
  p.set<double>("UniBo Hole c",        0.0,     "[1]"         );
  p.set<double>("UniBo Hole gamma",    2.16,    "[1]"         );
  p.set<double>("UniBo Hole mu0d",     90.0,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Hole mu0d_exp", -1.3,    "[1]"         );
  p.set<double>("UniBo Hole mu0a",     44.0,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Hole mu0a_exp", -0.7,    "[1]"         );
  p.set<double>("UniBo Hole mu1d",     28.2,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Hole mu1d_exp", -2.0,    "[1]"         );
  p.set<double>("UniBo Hole mu1a",     28.2,    "[cm^2/(V.s)]");
  p.set<double>("UniBo Hole mu1a_exp", -0.8,    "[1]"         );
  p.set<double>("UniBo Hole cr1",      1.3e18,  "[cm^-3]"     );
  p.set<double>("UniBo Hole cr1_exp",  2.2,     "[1]"         );
  p.set<double>("UniBo Hole cr2",      2.45e17, "[cm^-3]"     );
  p.set<double>("UniBo Hole cr2_exp",  3.1,     "[1]"         );
  p.set<double>("UniBo Hole cs1",      1.1e18,  "[cm^-3]"     );
  p.set<double>("UniBo Hole cs1_exp",  6.2,     "[1]"         );
  p.set<double>("UniBo Hole cs2",      6.1e20,  "[cm^-3]"     );
  p.set<double>("UniBo Hole alpha1",   0.77,    "[1]"         );
  p.set<double>("UniBo Hole alpha2",   0.719,   "[1]"         );

  //---------------------------------------------------------------------------
  // Masetti mobility model parameters (low-field).
  //--------------------------------------------------------------------------

  // Electrons (As)
  p.set<double>("Masetti Electron(As) mumax",  1417,    "mumax:[cm^2/(V.s)]" );
  p.set<double>("Masetti Electron(As) mumin1", 52.2,    "mumin1:[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(As) mumin2", 52.2,    "mumin2:[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(As) mu1",    43.4,    "mu1:[cm^2/(V.s)]"   );
  p.set<double>("Masetti Electron(As) gamma",  -2.45,   "gamma:[1]"          );
  p.set<double>("Masetti Electron(As) pc",     0,       "pc:[cm^-3]"         );
  p.set<double>("Masetti Electron(As) cr",     9.68e16, "cr:[cm^-3]"         );
  p.set<double>("Masetti Electron(As) cs",     3.43e20, "cs:[cm^-3]"         );
  p.set<double>("Masetti Electron(As) alpha",  0.680,   "alpha:[1]"          );
  p.set<double>("Masetti Electron(As) beta",   2.0,     "beta:[1]"           );

  // Electrons (P)
  p.set<double>("Masetti Electron(P) mumax",  1414,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(P) mumin1", 68.5,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(P) mumin2", 68.5,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(P) mu1",    56.1,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Electron(P) gamma",  -2.45,   "[1]"         );
  p.set<double>("Masetti Electron(P) pc",     0,       "[cm^-3]"     );
  p.set<double>("Masetti Electron(P) cr",     9.20e16, "[cm^-3]"     );
  p.set<double>("Masetti Electron(P) cs",     3.41e20, "[cm^-3]"     );
  p.set<double>("Masetti Electron(P) alpha",  0.711,   "[1]"         );
  p.set<double>("Masetti Electron(P) beta",   1.98,    "[1]"         );

  // Holes
  p.set<double>("Masetti Hole mumax",  470.5,   "[cm^2/(V.s)]");
  p.set<double>("Masetti Hole mumin1", 44.9,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Hole mumin2", 0,       "[cm^2/(V.s)]");
  p.set<double>("Masetti Hole mu1",    29.0,    "[cm^2/(V.s)]");
  p.set<double>("Masetti Hole gamma",  2.2,     "[1]"         );
  p.set<double>("Masetti Hole pc",     9.23e16, "[cm^-3]"     );
  p.set<double>("Masetti Hole cr",     2.23e17, "[cm^-3]"     );
  p.set<double>("Masetti Hole cs",     6.10e20, "[cm^-3]"     );
  p.set<double>("Masetti Hole alpha",  0.719,   "[1]"         );
  p.set<double>("Masetti Hole beta",   2.0,     "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons (As)
  p.set<double>("Philips Electron(As) mumax", 1417,    "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin", 52.2,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma", -2.285,  "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",  9.68e16, "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha", 0.68,    "alpha:[1]"         );

  // Electrons (P)
  p.set<double>("Philips Electron(P) mumax", 1414,   "[cm^2/(V.s)]");
  p.set<double>("Philips Electron(P) mumin", 68.5,   "[cm^2/(V.s)]");
  p.set<double>("Philips Electron(P) gamma", -2.285, "[1]"         );
  p.set<double>("Philips Electron(P) nref",  9.2e16, "[cm^-3]"     );
  p.set<double>("Philips Electron(P) alpha", 0.711,  "[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 470.5,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 44.9,    "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", -2.247,  "[1]"         );
  p.set<double>("Philips Hole nref",  2.23e17, "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 0.719,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    4e20,   "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 7.2e20, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    0.21,   "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 0.5,    "cref_a:[1]"    );

  //---------------------------------------------------------------------------
  // Klaassen mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Klaassen Electron mumax", 1417,    "mumax:[cm^2/(V.s)]");
  p.set<double>("Klaassen Electron mumin", 52.2,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Klaassen Electron gamma", -2.285,  "gamma:[1]"         );
  p.set<double>("Klaassen Electron nref",  9.68e16, "nref:[cm^-3]"      );
  p.set<double>("Klaassen Electron theta",2.285);
  p.set<double>("Klaassen Electron alpha",0.68);
  p.set<double>("Klaassen Electron ZNrefD",4.0e20);
  p.set<double>("Klaassen Electron ZcdD",0.21);
  p.set<double>("Klaassen Electron Vsat Alpha",2.4e7);
  p.set<double>("Klaassen Electron Vsat Beta",2.0);
  p.set<double>("Klaassen Electron Vsat Theta",0.8);
  p.set<double>("Klaassen Electron Vsat TNominal",600.0);

  // Holes
  p.set<double>("Klaassen Hole mumax", 470.5,   "[cm^2/(V.s)]");
  p.set<double>("Klaassen Hole mumin", 44.9,    "[cm^2/(V.s)]");
  p.set<double>("Klaassen Hole gamma", -2.247,  "[1]"         );
  p.set<double>("Klaassen Hole nref",  2.23e17, "[cm^-3]"     );
  p.set<double>("Klaassen Hole theta",2.247);
  p.set<double>("Klaassen Hole alpha",0.719);
  p.set<double>("Klaassen Hole ZNrefA",7.2e20);
  p.set<double>("Klaassen Hole ZcdA",0.5);
  p.set<double>("Klaassen Hole Vsat Alpha",2.4e7);
  p.set<double>("Klaassen Hole Vsat Beta",1.0);
  p.set<double>("Klaassen Hole Vsat Theta",0.8);
  p.set<double>("Klaassen Hole Vsat TNominal",600.0);

  //Generic
  p.set<double>("Klaassen s1", 0.89233);
  p.set<double>("Klaassen s2", 0.41372);
  p.set<double>("Klaassen s3", 0.19778);
  p.set<double>("Klaassen s4", 0.28227);
  p.set<double>("Klaassen s5", 0.005978);
  p.set<double>("Klaassen s6", 1.80618);
  p.set<double>("Klaassen s7", 0.72169);

  p.set<double>("Klaassen r1", 0.7643);
  p.set<double>("Klaassen r2", 2.2999);
  p.set<double>("Klaassen r3", 6.5502);
  p.set<double>("Klaassen r4", 2.3670);
  p.set<double>("Klaassen r5", -0.8552);
  p.set<double>("Klaassen r6", 0.6478);

  p.set<double>("Klaassen Pfc", 2.459);
  p.set<double>("Klaassen Pfb", 3.828);

  p.set<double>("Klaassen FME", 1.0);
  p.set<double>("Klaassen FMH", 1.258);


  //---------------------------------------------------------------------------
  // Klaassen mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Saturation Velocity", 7.7e6,    "[cm/s]");
  p.set<double>("Hole Saturation Velocity", 7.7e6,    "[cm/s]");

  //---------------------------------------------------------------------------
  // Shirahata mobility model parameters, carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Shirahata Electron muo", 1430,    "muo:[cm^2/(V.s)]");
  p.set<double>("Shirahata Electron theta",2.285,    "theta:[1]");
  p.set<double>("Shirahata Electron E1",6.3e3,    "E1:[V/cm]");
  p.set<double>("Shirahata Electron E2",7.7e5,    "E2:[V/cm]");
  p.set<double>("Shirahata Electron P1",0.28,    "P1:[1]");
  p.set<double>("Shirahata Electron P2",2.9,    "P2:[2]");

  // Holes
  p.set<double>("Shirahata Hole muo", 500,    "muo:[cm^2/(V.s)]");
  p.set<double>("Shirahata Hole theta",2.247,    "theta:[1]");
  p.set<double>("Shirahata Hole E1",8.0e3,    "E1:[V/cm]" );
  p.set<double>("Shirahata Hole E2",3.9e5,    "E2:[V/cm]");
  p.set<double>("Shirahata Hole P1",0.3,    "P1:[1]" );
  p.set<double>("Shirahata Hole P2",1.0,    "P2:[1]");

  //---------------------------------------------------------------------------
  // Selberherr avalanche model, which is the only one implemented in Charon 1.0.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Selberherr Electron a0",        7.03e5,     "a0_e:[1/cm] = N_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Electron a1",        0.0,        "a1_e:[1/(cm.K)] = N_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Electron a2",        0.0,        "a2_e:[1/(cm.K^2)] = N_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Electron delta",     1.0,        "delta_e:[1] = EXN_II in Charon 1.0"               );
  p.set<double>("Selberherr Electron E0",        1.231e6,    "E0_e:[V/cm] = ECN.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Electron lambda300", 10.4542e-7, "lambda300_e:[cm] = LAN300 in Charon 1.0"          );
  p.set<double>("Selberherr Electron hbarOmega", 0.063,      "hbarOmega_e:[eV] = OP_PH_EN in Charon 1.0"        );

  // Holes
  p.set<double>("Selberherr Hole a0",        1.528e6,    "a0_h:[1/cm] = P_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Hole a1",        0.0,        "a1_h:[1/(cm.K)] = P_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Hole a2",        0.0,        "a2_h:[1/(cm.K^2)] = P_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Hole delta",     1.0,        "delta_h:[1] = EXP_II in Charon 1.0"               );
  p.set<double>("Selberherr Hole E0",        2.036e6,    "E0_h:[V/cm] = ECP.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Hole lambda300", 6.32079e-7, "lambda300_h:[cm] = LAP300 in Charon 1.0"          );
  p.set<double>("Selberherr Hole hbarOmega", 0.063,      "hbarOmega_h:[eV] = OP_PH_EN in Charon 1.0"        );

  //---------------------------------------------------------------------------
  // Crowell-Sze avalanche model.
  //---------------------------------------------------------------------------

  p.set<double>("Crowell-Sze E_opt_ph",            0.063,   "E_opt_ph:[eV]");
  // Electrons
  p.set<double>("Crowell-Sze Electron lambda300",  6.2e-7,  "lambda300_e:[cm]");
  p.set<double>("Crowell-Sze Electron Ei",         1.1,     "Ei_e:[eV]");
  
  // Hole
  p.set<double>("Crowell-Sze Hole lambda300",      3.8e-7,  "lambda300_h:[cm]");
  p.set<double>("Crowell-Sze Hole Ei",             1.8,     "Ei_h:[eV]");

  //---------------------------------------------------------------------------
  // vanOverstraeten avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("vanOverstraeten Electron alow",  7.03e5,  "al_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron blow",  1.231e6, "bl_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron ahigh", 7.03e5,  "ah_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron bhigh", 1.231e6, "bh_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron E0",    4e5,     "E0_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron hbarOmega", 0.063, "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("vanOverstraeten Hole alow",      1.582e6, "al_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole blow",      2.036e6, "bl_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole ahigh",     6.71e5,  "ah_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole bhigh",     1.693e6, "bh_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole E0",        4e5,     "E0_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // Okuto-Crowell impact ionization model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Okuto Electron a",     0.426,   "a_e:[1/V]"  );
  p.set<double>("Okuto Electron b",     4.81e5,  "b_e:[V/cm]" );
  p.set<double>("Okuto Electorn c",     3.05e-4, "c_e:[1/K]"  );
  p.set<double>("Okuto Electorn d",     6.86e-4, "d_e:[1/K]"  );
  p.set<double>("Okuto Electorn gamma", 1,       "gamma_e:[1]");
  p.set<double>("Okuto Electorn delta", 2,       "delta_e:[1]");

  // Holes
  p.set<double>("Okuto Hole a",     0.243,   "a_h:[1/V]"  );
  p.set<double>("Okuto Hole b",     6.53e5,  "b_h:[V/cm]" );
  p.set<double>("Okuto Hole c",     5.35e-4, "c_h:[1/K]"  );
  p.set<double>("Okuto Hole d",     5.67e-4, "d_h:[1/K]"  );
  p.set<double>("Okuto Hole gamma", 1,       "gamma_h:[1]");
  p.set<double>("Okuto Hole delta", 2,       "delta_h:[1]");

  //---------------------------------------------------------------------------
  // Lackner avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Lackner Electron a",         1.316e6, "a_e:[1/cm]"      );
  p.set<double>("Lackner Electron b",         1.474e6, "b_e:[V/cm]"      );
  p.set<double>("Lackner Electron hbarOmega", 0.063,   "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("Lackner Hole a",         1.818e6, "a_h:[1/cm]"      );
  p.set<double>("Lackner Hole b",         2.036e6, "b_h:[V/cm]"      );
  p.set<double>("Lackner Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBo Electron a0", 4.3383,    "a0_e:[V]"   );
  p.set<double>("UniBo Electron a1", -2.42e-12, "a1_e:[V]"   );
  p.set<double>("UniBo Electron a2", 4.1233,    "a2_e:[1]"   );
  p.set<double>("UniBo Electron b0", 0.235,     "b0_e:[V]"   );
  p.set<double>("UniBo Electron b1", 0.0,       "b1_e:[1]"   );
  p.set<double>("UniBo Electron c0", 1.6831e4,  "c0_e:[V/cm]");
  p.set<double>("UniBo Electron c1", 4.3796,    "c1_e:[V/cm]");
  p.set<double>("UniBo Electron c2", 0.13005,   "c2_e:[V/cm]");
  p.set<double>("UniBo Electron d0", 1.2337e6,  "d0_e:[V/cm]");
  p.set<double>("UniBo Electron d1", 1.2039e3,  "d1_e:[V/cm]");
  p.set<double>("UniBo Electron d2", 0.56703,   "d2_e:[V/cm]");

  // Holes
  p.set<double>("UniBo Hole a0", 2.376,     "a0_h:[V]"   );
  p.set<double>("UniBo Hole a1", 1.033e-2,  "a1_h:[V]"   );
  p.set<double>("UniBo Hole a2", 0.0,       "a2_h:[1]"   );
  p.set<double>("UniBo Hole b0", 0.17714,   "b0_h:[V]"   );
  p.set<double>("UniBo Hole b1", -2.178e-3, "b1_h:[1]"   );
  p.set<double>("UniBo Hole c0", 9.47e-3,   "c0_h:[V/cm]");
  p.set<double>("UniBo Hole c1", 2.4924,    "c1_h:[1]"   );
  p.set<double>("UniBo Hole c2", 0.0,       "c2_h:[1]"   );
  p.set<double>("UniBo Hole d0", 1.4043e6,  "d0_h:[V/cm]");
  p.set<double>("UniBo Hole d1", 2.9744e3,  "d1_h:[V/cm]");
  p.set<double>("UniBo Hole d2", 1.4829,    "d2_h:[V/cm]");

  //---------------------------------------------------------------------------
  // New University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBoNew Electron a0",  4.65403,      "a0_e:[V]"   );
  p.set<double>("UniBoNew Electron a1",  -8.76031e-3,  "a1_e:[V]"   );
  p.set<double>("UniBoNew Electron a2",  1.34037e-5,   "a2_e:[V]"   );
  p.set<double>("UniBoNew Electron a3",  -2.75108e-9,  "a3_e:[V]"   );
  p.set<double>("UniBoNew Electron b0",  -0.128302,    "b0_e:[V]"   );
  p.set<double>("UniBoNew Electron b1",  4.45552e-3,   "b1_e:[V]"   );
  p.set<double>("UniBoNew Electron b2",  -1.0866e-5,   "b2_e:[V]"   );
  p.set<double>("UniBoNew Electron b3",  9.23119e-9,   "b3_e:[V]"   );
  p.set<double>("UniBoNew Electron b4",  -1.82482e-12, "b4_e:[V]"   );
  p.set<double>("UniBoNew Electron b5",  -4.82689e-15, "b5_e:[V]"   );
  p.set<double>("UniBoNew Electron b6",  1.09402e-17,  "b6_e:[V]"   );
  p.set<double>("UniBoNew Electron b7",  -1.24961e-20, "b7_e:[V]"   );
  p.set<double>("UniBoNew Electron b8",  7.55584e-24,  "b8_e:[V]"   );
  p.set<double>("UniBoNew Electron b9",  -2.28615e-27, "b9_e:[V]"   );
  p.set<double>("UniBoNew Electron b10", 2.73344e-31,  "b10_e:[V]"  );
  p.set<double>("UniBoNew Electron c0",  7.76221e3,    "c0_e:[V/cm]");
  p.set<double>("UniBoNew Electron c1",  25.18888,     "c1_e:[V/cm]");
  p.set<double>("UniBoNew Electron c2",  -1.37417e-3,  "c2_e:[V/cm]");
  p.set<double>("UniBoNew Electron c3",  1.59525e-4,   "c3_e:[V/cm]");
  p.set<double>("UniBoNew Electron d0",  7.10481e5,    "d0_e:[V/cm]");
  p.set<double>("UniBoNew Electron d1",  3.98594e3,    "d1_e:[V/cm]");
  p.set<double>("UniBoNew Electron d2",  -7.19956,     "d2_e:[V/cm]");
  p.set<double>("UniBoNew Electron d3",  6.96431e-3,   "d3_e:[V/cm]");

  // Holes
  p.set<double>("UniBoNew Hole a0",  2.26018,      "a0_h:[V]"   );
  p.set<double>("UniBoNew Hole a1",  0.0134001,    "a1_h:[V]"   );
  p.set<double>("UniBoNew Hole a2",  -5.87724e-6,  "a2_h:[V]"   );
  p.set<double>("UniBoNew Hole a3",  -1.14021e-9,  "a3_h:[V]"   );
  p.set<double>("UniBoNew Hole b0",  0.058547,     "b0_h:[V]"   );
  p.set<double>("UniBoNew Hole b1",  -1.95755e-4,  "b1_h:[V]"   );
  p.set<double>("UniBoNew Hole b2",  2.44357e-7,   "b2_h:[V]"   );
  p.set<double>("UniBoNew Hole b3",  -1.33202e-10, "b3_h:[V]"   );
  p.set<double>("UniBoNew Hole b4",  2.68082e-14,  "b4_h:[V]"   );
  p.set<double>("UniBoNew Hole b5",  0.0,          "b5_h:[V]"   );
  p.set<double>("UniBoNew Hole b6",  0.0,          "b6_h:[V]"   );
  p.set<double>("UniBoNew Hole b7",  0.0,          "b7_h:[V]"   );
  p.set<double>("UniBoNew Hole b8",  0.0,          "b8_h:[V]"   );
  p.set<double>("UniBoNew Hole b9",  0.0,          "b9_h:[V]"   );
  p.set<double>("UniBoNew Hole b10", 0.0,          "b10_h:[V]"  );
  p.set<double>("UniBoNew Hole c0",  1.95399e4,    "c0_h:[V/cm]");
  p.set<double>("UniBoNew Hole c1",  -104.441,     "c1_h:[V/cm]");
  p.set<double>("UniBoNew Hole c2",  0.498768,     "c2_h:[V/cm]");
  p.set<double>("UniBoNew Hole c3",  0.0,          "c3_h:[V/cm]");
  p.set<double>("UniBoNew Hole d0",  2.07712e6,    "d0_h:[V/cm]");
  p.set<double>("UniBoNew Hole d1",  993.153,      "d1_h:[V/cm]");
  p.set<double>("UniBoNew Hole d2",  7.77769,      "d2_h:[V/cm]");
  p.set<double>("UniBoNew Hole d3",  0.0,          "d3_h:[V/cm]");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  // Parameters for the "TempDep" heat capacity model
  p.set<double>("Heat Capacity a",        1.63,    "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,       "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,       "c:[J/(K^3.cm^3)]");

  // Parameters for the "PowerLawTempDep" heat capacity model
  p.set<double>("Mass Density",        2.33,  "rho:[g/cm^3]"   );
  p.set<double>("Heat Capacity c300",  0.711, "c300:[J/(K.g)]" );
  p.set<double>("Heat Capacity c1",    0.255, "c1:[J/(K.g)]"   );
  p.set<double>("Heat Capacity beta",  1.85,  "beta: [1]"      );

  // Parameters for the "TempDep" thermal conductivity model -------------------
  p.set<double>("Thermal Conductivity a", 0.03,    "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 1.56e-3, "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 1.65e-6, "c:[cm/(W.K)]"    );

  // Parameters for the "PowerLawTempDep" thermal conductivity model -----------
  p.set<double>("Thermal Conductivity kappa300", 1.48, "kappa300:[W/(K.cm)]");
  p.set<double>("Thermal Conductivity alpha",    -1.3, "alpha:[1]"          );
} // end of setSiliconParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setSiO2Parameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setSiO2Parameters(Teuchos::ParameterList& p)
{
  p.set<std::string>("Material Type", "Insulator", "");
  p.set<double>("Relative Permittivity", 3.9, "[1]" );
  p.set<double>("Electron Affinity",     1.0, "[eV]");
  p.set<double>("Band Gap",              9.0, "[eV]");

  p.set<double>("Mass Density",          2.65,"[g/cm^3]"   );
} // end of setSiO2Parameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setOxyNitrideParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setOxyNitrideParameters(
  Teuchos::ParameterList& p)
{
  p.set<std::string>("Material Type", "Insulator", "");
  p.set<double>("Relative Permittivity", 3.9, "[1]" );
  p.set<double>("Electron Affinity",     1.0, "[eV]");
  p.set<double>("Band Gap",              9.0, "[eV]");
} // end of setOxyNitrideParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setGaAsParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setGaAsParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          13.1,  "[1]"             );    // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.07,  "[eV]"            );    // otherwise noted.
  p.set<double>("Band Gap",                       1.424, "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.8e6, "[cm^-3]"         );
  p.set<double>("Electron Effective Mass",        0.067, "[1]"             );
  p.set<double>("Hole Effective Mass",            0.082, "mlh:[1]"         );
  p.set<double>("Electron Mobility",              8500., "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  400.,  "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 220.,  "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     10.,   "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-9,  "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-9,  "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Hole SRH TPowerLaw",             -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,  "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,  "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,     "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.0e-10,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    1.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    1.0e-30,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.07,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.424,    "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 4.7e17, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     7.0e18, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,    "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 8500.0,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 0.0,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  1.69e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", -1.0,    "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,     "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 0.436,   "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 400.0,   "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,     "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  2.75e17, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", -2.1,    "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,     "[1]"         );
  p.set<double>("Analytic Hole alpha", 0.395,   "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    8500.0,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    0.0,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1.26e17, "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 1.0,     "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     -0.57,   "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     0.0,     "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     0.0,     "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     0.0,     "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",    400.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",    0.0,     "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",     2.35e17, "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp", 1.0,     "[1]"         );
  p.set<double>("Arora Hole exp1",     0.0,     "[1]"         );
  p.set<double>("Arora Hole exp2",     0.0,     "[1]"         );
  p.set<double>("Arora Hole exp3",     0.0,     "[1]"         );
  p.set<double>("Arora Hole exp4",     0.0,     "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax", 8500.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin", 0.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma", 0.0,    "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",  1e30,   "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha", 1.0,    "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 400.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );

  //---------------------------------------------------------------------------
  // Field-dependent mobility parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Saturation Velocity", 7.7e6, "[cm/s]");
  p.set<double>("Electron Saturation Field",    4e3,   "[V/cm]");
  p.set<double>("Hole Saturation Velocity",     7.7e6, "[cm/s]");
  p.set<double>("Hole Saturation Field",        4e3,   "[V/cm]");

  //---------------------------------------------------------------------------
  // Selberherr avalanche model
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Selberherr Electron a0",        7.03e5,     "a0_e:[1/cm] = N_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Electron a1",        0.0,        "a1_e:[1/(cm.K)] = N_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Electron a2",        0.0,        "a2_e:[1/(cm.K^2)] = N_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Electron delta",     1.0,        "delta_e:[1] = EXN_II in Charon 1.0"               );
  p.set<double>("Selberherr Electron E0",        1.231e6,    "E0_e:[V/cm] = ECN.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Electron lambda300", 10.4542e-7, "lambda300_e:[cm] = LAN300 in Charon 1.0"          );
  p.set<double>("Selberherr Electron hbarOmega", 0.036,      "hbarOmega_e:[eV] = OP_PH_EN in Charon 1.0"        );

  // Holes
  p.set<double>("Selberherr Hole a0",        1.582e6,    "a0_h:[1/cm] = P_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Hole a1",        0.0,        "a1_h:[1/(cm.K)] = P_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Hole a2",        0.0,        "a2_h:[1/(cm.K^2)] = P_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Hole delta",     1.0,        "delta_h:[1] = EXP_II in Charon 1.0"               );
  p.set<double>("Selberherr Hole E0",        2.036e6,    "E0_h:[V/cm] = ECP.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Hole lambda300", 6.32079e-7, "lambda300_h:[cm] = LAP300 in Charon 1.0"          );
  p.set<double>("Selberherr Hole hbarOmega", 0.036,      "hbarOmega_h:[eV] = OP_PH_EN in Charon 1.0"        );

  //---------------------------------------------------------------------------
  // vanOverstraeten avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("vanOverstraeten Electron alow",  7.03e5,  "al_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron blow",  1.231e6, "bl_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron ahigh", 7.03e5,  "ah_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron bhigh", 1.231e6, "bh_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron E0",    4e5,     "E0_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron hbarOmega", 0.036, "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("vanOverstraeten Hole alow",      1.582e6, "al_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole blow",      2.036e6, "bl_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole ahigh",     6.71e5,  "ah_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole bhigh",     1.693e6, "bh_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole E0",        4e5,     "E0_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole hbarOmega", 0.036,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // Okuto-Crowell impact ionization model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Okuto Electron a",     0.426,   "a_e:[1/V]"  );
  p.set<double>("Okuto Electron b",     4.81e5,  "b_e:[V/cm]" );
  p.set<double>("Okuto Electorn c",     3.05e-4, "c_e:[1/K]"  );
  p.set<double>("Okuto Electorn d",     6.86e-4, "d_e:[1/K]"  );
  p.set<double>("Okuto Electorn gamma", 1,       "gamma_e:[1]");
  p.set<double>("Okuto Electorn delta", 2,       "delta_e:[1]");

  // Holes
  p.set<double>("Okuto Hole a",     0.243,   "a_h:[1/V]"  );
  p.set<double>("Okuto Hole b",     6.53e5,  "b_h:[V/cm]" );
  p.set<double>("Okuto Hole c",     5.35e-4, "c_h:[1/K]"  );
  p.set<double>("Okuto Hole d",     5.67e-4, "d_h:[1/K]"  );
  p.set<double>("Okuto Hole gamma", 1,       "gamma_h:[1]");
  p.set<double>("Okuto Hole delta", 2,       "delta_h:[1]");

  //---------------------------------------------------------------------------
  // Lackner avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Lackner Electron a",         1.316e6, "a_e:[1/cm]"      );
  p.set<double>("Lackner Electron b",         1.474e6, "b_e:[V/cm]"      );
  p.set<double>("Lackner Electron hbarOmega", 0.063,   "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("Lackner Hole a",         1.818e6, "a_h:[1/cm]"      );
  p.set<double>("Lackner Hole b",         2.036e6, "b_h:[V/cm]"      );
  p.set<double>("Lackner Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBo Electron a0", 4.3383,    "a0_e:[V]"   );
  p.set<double>("UniBo Electron a1", -2.42e-12, "a1_e:[V]"   );
  p.set<double>("UniBo Electron a2", 4.1233,    "a2_e:[1]"   );
  p.set<double>("UniBo Electron b0", 0.235,     "b0_e:[V]"   );
  p.set<double>("UniBo Electron b1", 0.0,       "b1_e:[1]"   );
  p.set<double>("UniBo Electron c0", 1.6831e4,  "c0_e:[V/cm]");
  p.set<double>("UniBo Electron c1", 4.3796,    "c1_e:[V/cm]");
  p.set<double>("UniBo Electron c2", 0.13005,   "c2_e:[V/cm]");
  p.set<double>("UniBo Electron d0", 1.2337e6,  "d0_e:[V/cm]");
  p.set<double>("UniBo Electron d1", 1.2039e3,  "d1_e:[V/cm]");
  p.set<double>("UniBo Electron d2", 0.56703,   "d2_e:[V/cm]");

  // Holes
  p.set<double>("UniBo Hole a0", 2.376,     "a0_h:[V]"   );
  p.set<double>("UniBo Hole a1", 1.033e-2,  "a1_h:[V]"   );
  p.set<double>("UniBo Hole a2", 0.0,       "a2_h:[1]"   );
  p.set<double>("UniBo Hole b0", 0.17714,   "b0_h:[V]"   );
  p.set<double>("UniBo Hole b1", -2.178e-3, "b1_h:[1]"   );
  p.set<double>("UniBo Hole c0", 9.47e-3,   "c0_h:[V/cm]");
  p.set<double>("UniBo Hole c1", 2.4924,    "c1_h:[1]"   );
  p.set<double>("UniBo Hole c2", 0.0,       "c2_h:[1]"   );
  p.set<double>("UniBo Hole d0", 1.4043e6,  "d0_h:[V/cm]");
  p.set<double>("UniBo Hole d1", 2.9744e3,  "d1_h:[V/cm]");
  p.set<double>("UniBo Hole d2", 1.4829,    "d2_h:[V/cm]");

  //---------------------------------------------------------------------------
  // New University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBoNew Electron a0",  4.65403,      "a0_e:[V]"   );
  p.set<double>("UniBoNew Electron a1",  -8.76031e-3,  "a1_e:[V]"   );
  p.set<double>("UniBoNew Electron a2",  1.34037e-5,   "a2_e:[V]"   );
  p.set<double>("UniBoNew Electron a3",  -2.75108e-9,  "a3_e:[V]"   );
  p.set<double>("UniBoNew Electron b0",  -0.128302,    "b0_e:[V]"   );
  p.set<double>("UniBoNew Electron b1",  4.45552e-3,   "b1_e:[V]"   );
  p.set<double>("UniBoNew Electron b2",  -1.0866e-5,   "b2_e:[V]"   );
  p.set<double>("UniBoNew Electron b3",  9.23119e-9,   "b3_e:[V]"   );
  p.set<double>("UniBoNew Electron b4",  -1.82482e-12, "b4_e:[V]"   );
  p.set<double>("UniBoNew Electron b5",  -4.82689e-15, "b5_e:[V]"   );
  p.set<double>("UniBoNew Electron b6",  1.09402e-17,  "b6_e:[V]"   );
  p.set<double>("UniBoNew Electron b7",  -1.24961e-20, "b7_e:[V]"   );
  p.set<double>("UniBoNew Electron b8",  7.55584e-24,  "b8_e:[V]"   );
  p.set<double>("UniBoNew Electron b9",  -2.28615e-27, "b9_e:[V]"   );
  p.set<double>("UniBoNew Electron b10", 2.73344e-31,  "b10_e:[V]"  );
  p.set<double>("UniBoNew Electron c0",  7.76221e3,    "c0_e:[V/cm]");
  p.set<double>("UniBoNew Electron c1",  25.18888,     "c1_e:[V/cm]");
  p.set<double>("UniBoNew Electron c2",  -1.37417e-3,  "c2_e:[V/cm]");
  p.set<double>("UniBoNew Electron c3",  1.59525e-4,   "c3_e:[V/cm]");
  p.set<double>("UniBoNew Electron d0",  7.10481e5,    "d0_e:[V/cm]");
  p.set<double>("UniBoNew Electron d1",  3.98594e3,    "d1_e:[V/cm]");
  p.set<double>("UniBoNew Electron d2",  -7.19956,     "d2_e:[V/cm]");
  p.set<double>("UniBoNew Electron d3",  6.96431e-3,   "d3_e:[V/cm]");

  // Holes
  p.set<double>("UniBoNew Hole a0",  2.26018,      "a0_h:[V]"   );
  p.set<double>("UniBoNew Hole a1",  0.0134001,    "a1_h:[V]"   );
  p.set<double>("UniBoNew Hole a2",  -5.87724e-6,  "a2_h:[V]"   );
  p.set<double>("UniBoNew Hole a3",  -1.14021e-9,  "a3_h:[V]"   );
  p.set<double>("UniBoNew Hole b0",  0.058547,     "b0_h:[V]"   );
  p.set<double>("UniBoNew Hole b1",  -1.95755e-4,  "b1_h:[V]"   );
  p.set<double>("UniBoNew Hole b2",  2.44357e-7,   "b2_h:[V]"   );
  p.set<double>("UniBoNew Hole b3",  -1.33202e-10, "b3_h:[V]"   );
  p.set<double>("UniBoNew Hole b4",  2.68082e-14,  "b4_h:[V]"   );
  p.set<double>("UniBoNew Hole b5",  0.0,          "b5_h:[V]"   );
  p.set<double>("UniBoNew Hole b6",  0.0,          "b6_h:[V]"   );
  p.set<double>("UniBoNew Hole b7",  0.0,          "b7_h:[V]"   );
  p.set<double>("UniBoNew Hole b8",  0.0,          "b8_h:[V]"   );
  p.set<double>("UniBoNew Hole b9",  0.0,          "b9_h:[V]"   );
  p.set<double>("UniBoNew Hole b10", 0.0,          "b10_h:[V]"  );
  p.set<double>("UniBoNew Hole c0",  1.95399e4,    "c0_h:[V/cm]");
  p.set<double>("UniBoNew Hole c1",  -104.441,     "c1_h:[V/cm]");
  p.set<double>("UniBoNew Hole c2",  0.498768,     "c2_h:[V/cm]");
  p.set<double>("UniBoNew Hole c3",  0.0,          "c3_h:[V/cm]");
  p.set<double>("UniBoNew Hole d0",  2.07712e6,    "d0_h:[V/cm]");
  p.set<double>("UniBoNew Hole d1",  993.153,      "d1_h:[V/cm]");
  p.set<double>("UniBoNew Hole d2",  7.77769,      "d2_h:[V/cm]");
  p.set<double>("UniBoNew Hole d3",  0.0,          "d3_h:[V/cm]");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  // Parameters for the "PowerLawTempDep" heat capacity model
  p.set<double>("Mass Density",        5.32,  "rho:[g/cm^3]"   );
  p.set<double>("Heat Capacity c300",  0.322, "c300:[J/(K.g)]" );
  p.set<double>("Heat Capacity c1",    0.05,  "c1:[J/(K.g)]"   );
  p.set<double>("Heat Capacity beta",  1.6,   "beta: [1]"      );

  // Parameters for the "PowerLawTempDep" thermal conductivity model -----------
  p.set<double>("Thermal Conductivity kappa300", 0.46,  "kappa300:[W/(K.cm)]");
  p.set<double>("Thermal Conductivity alpha",    -1.25, "alpha:[1]"          );

} // end of setGaAsParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setInPParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setInPParameters(Teuchos::ParameterList& p)
{
  p.set<std::string>("Material Type", "Semiconductor", "");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  // Parameters for the "PowerLawTempDep" heat capacity model
  p.set<double>("Mass Density",        4.81,  "rho:[g/cm^3]"   );
  p.set<double>("Heat Capacity c300",  0.41,  "c300:[J/(K.g)]" );
  p.set<double>("Heat Capacity c1",    0.05,  "c1:[J/(K.g)]"   );
  p.set<double>("Heat Capacity beta",  2.05,  "beta: [1]"      );

  // Parameters for the "PowerLawTempDep" thermal conductivity model -----------
  p.set<double>("Thermal Conductivity kappa300", 0.68, "kappa300:[W/(K.cm)]");
  p.set<double>("Thermal Conductivity alpha",    -1.4, "alpha:[1]"          );

} // end of setInPParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setAlGaAsParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setAlGaAsParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          13.1,  "[1]"             );    // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.07,  "[eV]"            );    // otherwise noted.
  p.set<double>("Band Gap",                       1.424, "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.8e6, "[cm^-3]"         );
  p.set<double>("Electron Mobility",              9890., "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  400.,  "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 220.,  "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     10.,   "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-9,  "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-9,  "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Hole SRH TPowerLaw",             -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,  "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,  "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,     "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.0e-10,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    1.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    1.0e-30,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.07,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.424,    "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 4.7e17, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     7.0e18, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,    "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 9890.0,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 2.37e3,  "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  3.63e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", 0.0,     "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,     "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 400.0,   "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,     "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  2.75e17, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", 0.0,     "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,     "[1]"         );
  p.set<double>("Analytic Hole alpha", 0.395,   "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    9890.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    0.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1e20,   "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 1.0,    "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     0.0,    "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     0.0,    "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     0.0,    "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     0.0,    "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",     400.0, "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",     0.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",      1e20,  "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp",  1.0,   "[1]"         );
  p.set<double>("Arora Hole exp1",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp2",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp3",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp4",      0.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax", 9.892e3, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin", 0.0,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma", 0.0,     "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",  1e30,    "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha", 1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 400.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );
} // end of setAlGaAsParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setInGaAsParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setInGaAsParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          14.0,     "[1]"          );    // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.6,      "[eV]"         );    // otherwise noted.
  p.set<double>("Band Gap",                       0.75,     "[eV]"         );
  p.set<double>("Intrinsic Concentration",        7.418e11, "[cm^-3]"      );
  p.set<double>("Electron Mobility",              3000.,    "[cm^2/(V.s)]" );
  p.set<double>("Hole Mobility",                  76.,      "[cm^2/(V.s)]" );
  p.set<double>("Electron Diffusion Coefficient", 220.,     "[cm^2/s]"     );
  p.set<double>("Hole Diffusion Coefficient",     10.,      "[cm^2/s]"     );
  p.set<double>("Electron Lifetime",              1e-9,     "Tau0:[s]"     );
  p.set<double>("Hole Lifetime",                  1e-9,     "Tau0:[s]"     );
  p.set<double>("Electron SRH Conc",              5e16,     "Nsrh:[cm^-3]" );
  p.set<double>("Hole SRH Conc",                  5e16,     "Nsrh:[cm^-3]" );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,     "TPowerLaw:[1]");
  p.set<double>("Hole SRH TPowerLaw",             -1.5,     "TPowerLaw:[1]");
  p.set<double>("SRH Etrap",                      0,        "[eV]"         );
  p.set<double>(
    "Electron SRH TExponential",
    2.55,
    "TExponential:[1]");
  p.set<double>(
    "Hole SRH TExponential",
    2.55,
    "TExponential:[1]");
  p.set<double>(
    "Radiative Recombination Coefficient",
    5.0e-11,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    5.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    2.0e-29,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.6,      "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          0.75,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 2.0833e17, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,       "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS Exponent",     1.5,       "Nv_F:[1]"     );
  p.set<double>(
    "Hole Effective DOS at 300 K",
    1.02408e19,
    "Nv300:[cm^-3]");

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 2.73e4,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 4e3,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  3.63e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", 0.0,     "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,     "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 480.0,  "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,    "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  1.0e30, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", 0.0,    "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,    "[1]"         );
  p.set<double>("Analytic Hole alpha", 1.0,    "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    2.73e4, "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    0.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1e20,   "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 1.0,    "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     0.0,    "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     0.0,    "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     0.0,    "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     0.0,    "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",    480.0, "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",    0.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",     1e20,  "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp", 1.0,   "[1]"         );
  p.set<double>("Arora Hole exp1",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp2",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp3",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp4",     0.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax",  2.725e4, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin",  0.0,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma",  0.0,     "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",   1e30,    "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha",  1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 400.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );
} // end of setInGaAsParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setAlInAsParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setAlInAsParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          12.5,   "[1]"             );   // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.1,    "[eV]"            );   // otherwise noted.
  p.set<double>("Band Gap",                       1.46,   "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.39e6, "[cm^-3]"         );
  p.set<double>("Electron Mobility",              1700.,  "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  27.,    "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 220.,   "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     10.,    "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-9,   "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-9,   "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,   "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,   "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,   "TPowerLaw:[1]"   );
  p.set<double>("Hole SRH TPowerLaw",             -1.5,   "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,   "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,   "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,      "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.0e-10,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    1.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    1.0e-30,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.1,      "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.46,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 5.02e17,   "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,       "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1.2249e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,       "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 2.41e4, "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 497.0,  "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  1.0e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", 0.0,    "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,    "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 1.0,    "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 480.0,  "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,    "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  1.0e30, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", 0.0,    "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,    "[1]"         );
  p.set<double>("Analytic Hole alpha", 1.0,    "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    2.41e4, "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    0.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1e20,   "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 1.0,    "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     0.0,    "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     0.0,    "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     0.0,    "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     0.0,    "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",    480.0, "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",    0.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",     1e20,  "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp", 1.0,   "[1]"         );
  p.set<double>("Arora Hole exp1",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp2",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp3",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp4",     0.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax",  2.414e4, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin",  0.0,     "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma",  0.0,     "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",   1e30,    "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha",  1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 400.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );
} // end of setAlInAsParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setGaAsPParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setGaAsPParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          13.1,  "[1]"             );    // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.07,  "[eV]"            );    // otherwise noted.
  p.set<double>("Band Gap",                       1.424, "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.8e6, "[cm^-3]"         );
  p.set<double>("Electron Mobility",              200.0, "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  150.0, "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 220.,  "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     10.,   "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-9,  "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-9,  "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Hole SRH TPowerLaw",             -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,  "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,  "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,     "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.0e-10,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    1.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    1.0e-30,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.07,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.424,    "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 4.7e17, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     7.0e18, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,    "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 200.0,   "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 95.0,    "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  3.63e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", 0.0,     "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,     "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 1.0,     "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 150.0,  "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,    "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  1.0e30, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", 0.0,    "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,    "[1]         ");
  p.set<double>("Analytic Hole alpha", 1.0,    "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",    200.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",    0.0,   "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",     1e20,  "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp", 1.0,   "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",     0.0,   "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",     0.0,   "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",     0.0,   "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",     0.0,   "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",    150.0, "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",    0.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",     1e20,  "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp", 1.0,   "[1]"         );
  p.set<double>("Arora Hole exp1",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp2",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp3",     0.0,   "[1]"         );
  p.set<double>("Arora Hole exp4",     0.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax",  200.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin",  0.0,   "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma",  0.0,   "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",   1e30,  "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha",  1.0,   "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 150.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );
} // end of setGaAsPParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setInGaPParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setInGaPParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          11.8,  "[1]"             );    // The following are from Charon 1.0 unless
  p.set<double>("Electron Affinity",              4.04,  "[eV]"            );    // otherwise noted.
  p.set<double>("Band Gap",                       1.86,  "[eV]"            );
  p.set<double>("Intrinsic Concentration",        1.2e3, "[cm^-3]"         );
  p.set<double>("Electron Effective Mass",        0.067, "[1]"             );
  p.set<double>("Hole Effective Mass",            0.082, "mlh:[1]"         );
  p.set<double>("Electron Mobility",              5400., "[cm^2/(V.s)]"    );
  p.set<double>("Hole Mobility",                  200.0, "[cm^2/(V.s)]"    );
  p.set<double>("Electron Diffusion Coefficient", 130.0, "[cm^2/s]"        );
  p.set<double>("Hole Diffusion Coefficient",     5.0,   "[cm^2/s]"        );
  p.set<double>("Electron Lifetime",              1e-9,  "Tau0:[s]"        );
  p.set<double>("Hole Lifetime",                  1e-9,  "Tau0:[s]"        );
  p.set<double>("Electron SRH Conc",              5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Hole SRH Conc",                  5e16,  "Nsrh:[cm^-3]"    );
  p.set<double>("Electron SRH TPowerLaw",         -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Hole SRH TPowerLaw",             -1.5,  "TPowerLaw:[1]"   );
  p.set<double>("Electron SRH TExponential",      2.55,  "TExponential:[1]");
  p.set<double>("Hole SRH TExponential",          2.55,  "TExponential:[1]");
  p.set<double>("SRH Etrap",                      0,     "[eV]"            );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.0e-10,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    1.0e-30,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    1.0e-30,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.04,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          1.86,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             5.405e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              204.0,    "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 4.7e17, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     7.0e18, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,    "Nv_F:[1]"     );

  // Old Slotboom
  p.set<double>("Old Slotboom BGN V0_BGN",  0.0,  "V0_BGN:[V]"    );
  p.set<double>("Old Slotboom BGN N0_BGN",  1e17, "N0_BGN:[cm^-3]");
  p.set<double>("Old Slotboom BGN CON_BGN", 0.0,  "CON_BGN:[1]"   );

  // Persson
  p.set<double>("Persson ANC_BGN", -16.30e-3, "ANC_BGN:[eV]");
  p.set<double>("Persson BNC_BGN", 0.0,       "BNC_BGN:[eV]");
  p.set<double>("Persson CNC_BGN", -18.13e-3, "CNC_BGN:[eV]");
  p.set<double>("Persson ANV_BGN", 0.0,       "ANV_BGN:[eV]");
  p.set<double>("Persson BNV_BGN", 7.47e-3,   "BNV_BGN:[eV]");
  p.set<double>("Persson CNV_BGN", 72.52e-3,  "CNV_BGN:[eV]");
  p.set<double>("Persson APC_BGN", -9.71e-3,  "APC_BGN:[eV]");
  p.set<double>("Persson BPC_BGN", 0.0,       "BPC_BGN:[eV]");
  p.set<double>("Persson CPC_BGN", -0.47e-3,  "CPC_BGN:[eV]");
  p.set<double>("Persson APV_BGN", 0.0,       "APV_BGN:[eV]");
  p.set<double>("Persson BPV_BGN", 12.19e-3,  "BPV_BGN:[eV]");
  p.set<double>("Persson CPV_BGN", 3.41e-3,   "CPV_BGN:[eV]");

  //---------------------------------------------------------------------------
  // Analytic mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Analytic Electron mumax", 200.0,  "mumax:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron mumin", 0.95,   "mumin:[cm^2/(V.s)]");
  p.set<double>("Analytic Electron nref",  1.0e17, "nref:[cm^-3]"      );
  p.set<double>("Analytic Electron gamma", 0.0,    "gamma:[1]"         );
  p.set<double>("Analytic Electron xin",   0.0,    "xin:[1]"           );
  p.set<double>("Analytic Electron alpha", 1.0,    "alpha:[1]"         );

  // Holes
  p.set<double>("Analytic Hole mumax", 150.0,  "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole mumin", 0.0,    "[cm^2/(V.s)]");
  p.set<double>("Analytic Hole nref",  1.0e30, "[cm^-3]"     );
  p.set<double>("Analytic Hole gamma", 0.0,    "[1]"         );
  p.set<double>("Analytic Hole xin",   0.0,    "[1]"         );
  p.set<double>("Analytic Hole alpha", 1.0,    "[1]"         );

  //---------------------------------------------------------------------------
  // Arora mobility model parameters (low-field).
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Arora Electron mumax",     200.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Arora Electron mumin",     0.0,   "mumin:[cm^2/(V.s)]");
  p.set<double>("Arora Electron nref",      1e20,  "nref:[cm^-3]"      );
  p.set<double>("Arora Electron nref_exp",  1.0,   "nref_exp:[1]"      );
  p.set<double>("Arora Electron exp1",      0.0,   "exp1:[1]"          );
  p.set<double>("Arora Electron exp2",      0.0,   "exp2:[1]"          );
  p.set<double>("Arora Electron exp3",      0.0,   "exp3:[1]"          );
  p.set<double>("Arora Electron exp4",      0.0,   "exp4:[1]"          );

  // Holes
  p.set<double>("Arora Hole mumax",     150.0, "[cm^2/(V.s)]");
  p.set<double>("Arora Hole mumin",     0.0,   "[cm^2/(V.s)]");
  p.set<double>("Arora Hole nref",      1e20,  "[cm^-3]"     );
  p.set<double>("Arora Hole nref_exp",  1.0,   "[1]"         );
  p.set<double>("Arora Hole exp1",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp2",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp3",      0.0,   "[1]"         );
  p.set<double>("Arora Hole exp4",      0.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-dependent.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Philips Electron(As) mumax",  200.0, "mumax:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) mumin",  0.0,   "mumin:[cm^2/(V.s)]");
  p.set<double>("Philips Electron(As) gamma",  0.0,   "gamma:[1]"         );
  p.set<double>("Philips Electron(As) nref",   1e30,  "nref:[cm^-3]"      );
  p.set<double>("Philips Electron(As) alpha",  1.0,   "alpha:[1]"         );

  // Holes
  p.set<double>("Philips Hole mumax", 150.0, "[cm^2/(V.s)]");
  p.set<double>("Philips Hole mumin", 0.0,   "[cm^2/(V.s)]");
  p.set<double>("Philips Hole gamma", 0.0,   "[1]"         );
  p.set<double>("Philips Hole nref",  1e30,  "[cm^-3]"     );
  p.set<double>("Philips Hole alpha", 1.0,   "[1]"         );

  //---------------------------------------------------------------------------
  // Philips mobility model parameters (low-field), carrier-independent.
  //---------------------------------------------------------------------------
  p.set<double>("Philips Donor nref_d",    1e30, "nref_d:[cm^-3]");
  p.set<double>("Philips Acceptor nref_a", 1e30, "nref_a:[cm^-3]");
  p.set<double>("Philips Donor cref_d",    1e30, "cref_d:[1]"    );
  p.set<double>("Philips Acceptor cref_a", 1e30, "cref_a:[1]"    );

  //---------------------------------------------------------------------------
  // Field-dependent mobility parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Saturation Velocity", 7.7e6, "[cm/s]");
  p.set<double>("Electron Saturation Field",    4e3,   "[V/cm]");
  p.set<double>("Hole Saturation Velocity",     7.7e6, "[cm/s]");
  p.set<double>("Hole Saturation Field",        4e3,   "[V/cm]");

  //---------------------------------------------------------------------------
  // Selberherr avalanche model
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Selberherr Electron a0",        7.03e5,     "a0_e:[1/cm] = N_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Electron a1",        0.0,        "a1_e:[1/(cm.K)] = N_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Electron a2",        0.0,        "a2_e:[1/(cm.K^2)] = N_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Electron delta",     1.0,        "delta_e:[1] = EXN_II in Charon 1.0"               );
  p.set<double>("Selberherr Electron E0",        1.231e6,    "E0_e:[V/cm] = ECN.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Electron lambda300", 10.4542e-7, "lambda300_e:[cm] = LAN300 in Charon 1.0"          );
  p.set<double>("Selberherr Electron hbarOmega", 0.063,      "hbarOmega_e:[eV] = OP_PH_EN in Charon 1.0"        );

  // Holes
  p.set<double>("Selberherr Hole a0",        1.528e6,    "a0_h:[1/cm] = P_IONIZATION in Charon 1.0"         );
  p.set<double>("Selberherr Hole a1",        0.0,        "a1_h:[1/(cm.K)] = P_ION_1 in Charon 1.0"          );
  p.set<double>("Selberherr Hole a2",        0.0,        "a2_h:[1/(cm.K^2)] = P_ION_2 in Charon 1.0"        );
  p.set<double>("Selberherr Hole delta",     1.0,        "delta_h:[1] = EXP_II in Charon 1.0"               );
  p.set<double>("Selberherr Hole E0",        2.036e6,    "E0_h:[V/cm] = ECP.II in Selberherr, not in Charon 1.0");
  p.set<double>("Selberherr Hole lambda300", 6.32079e-7, "lambda300_h:[cm] = LAP300 in Charon 1.0"          );
  p.set<double>("Selberherr Hole hbarOmega", 0.063,      "hbarOmega_h:[eV] = OP_PH_EN in Charon 1.0"        );

  //---------------------------------------------------------------------------
  // vanOverstraeten avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("vanOverstraeten Electron alow",  7.03e5,  "al_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron blow",  1.231e6, "bl_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron ahigh", 7.03e5,  "ah_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron bhigh", 1.231e6, "bh_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron E0",    4e5,     "E0_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron hbarOmega", 0.036, "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("vanOverstraeten Hole alow",      1.582e6, "al_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole blow",      2.036e6, "bl_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole ahigh",     6.71e5,  "ah_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole bhigh",     1.693e6, "bh_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole E0",        4e5,     "E0_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole hbarOmega", 0.036,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // Okuto-Crowell avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Okuto Electron a",     0.426,   "a_e:[1/V]"  );
  p.set<double>("Okuto Electron b",     4.81e5,  "b_e:[V/cm]" );
  p.set<double>("Okuto Electorn c",     3.05e-4, "c_e:[1/K]"  );
  p.set<double>("Okuto Electorn d",     6.86e-4, "d_e:[1/K]"  );
  p.set<double>("Okuto Electorn gamma", 1,       "gamma_e:[1]");
  p.set<double>("Okuto Electorn delta", 2,       "delta_e:[1]");

  // Holes
  p.set<double>("Okuto Hole a",     0.243,   "a_h:[1/V]"  );
  p.set<double>("Okuto Hole b",     6.53e5,  "b_h:[V/cm]" );
  p.set<double>("Okuto Hole c",     5.35e-4, "c_h:[1/K]"  );
  p.set<double>("Okuto Hole d",     5.67e-4, "d_h:[1/K]"  );
  p.set<double>("Okuto Hole gamma", 1,       "gamma_h:[1]");
  p.set<double>("Okuto Hole delta", 2,       "delta_h:[1]");

  //---------------------------------------------------------------------------
  // Lackner avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Lackner Electron a",         1.316e6, "a_e:[1/cm]"      );
  p.set<double>("Lackner Electron b",         1.474e6, "b_e:[V/cm]"      );
  p.set<double>("Lackner Electron hbarOmega", 0.063,   "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("Lackner Hole a",         1.818e6, "a_h:[1/cm]"      );
  p.set<double>("Lackner Hole b",         2.036e6, "b_h:[V/cm]"      );
  p.set<double>("Lackner Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBo Electron a0", 4.3383,    "a0_e:[V]"   );
  p.set<double>("UniBo Electron a1", -2.42e-12, "a1_e:[V]"   );
  p.set<double>("UniBo Electron a2", 4.1233,    "a2_e:[1]"   );
  p.set<double>("UniBo Electron b0", 0.235,     "b0_e:[V]"   );
  p.set<double>("UniBo Electron b1", 0.0,       "b1_e:[1]"   );
  p.set<double>("UniBo Electron c0", 1.6831e4,  "c0_e:[V/cm]");
  p.set<double>("UniBo Electron c1", 4.3796,    "c1_e:[V/cm]");
  p.set<double>("UniBo Electron c2", 0.13005,   "c2_e:[V/cm]");
  p.set<double>("UniBo Electron d0", 1.2337e6,  "d0_e:[V/cm]");
  p.set<double>("UniBo Electron d1", 1.2039e3,  "d1_e:[V/cm]");
  p.set<double>("UniBo Electron d2", 0.56703,   "d2_e:[V/cm]");

  // Holes
  p.set<double>("UniBo Hole a0", 2.376,     "a0_h:[V]"   );
  p.set<double>("UniBo Hole a1", 1.033e-2,  "a1_h:[V]"   );
  p.set<double>("UniBo Hole a2", 0.0,       "a2_h:[1]"   );
  p.set<double>("UniBo Hole b0", 0.17714,   "b0_h:[V]"   );
  p.set<double>("UniBo Hole b1", -2.178e-3, "b1_h:[1]"   );
  p.set<double>("UniBo Hole c0", 9.47e-3,   "c0_h:[V/cm]");
  p.set<double>("UniBo Hole c1", 2.4924,    "c1_h:[1]"   );
  p.set<double>("UniBo Hole c2", 0.0,       "c2_h:[1]"   );
  p.set<double>("UniBo Hole d0", 1.4043e6,  "d0_h:[V/cm]");
  p.set<double>("UniBo Hole d1", 2.9744e3,  "d1_h:[V/cm]");
  p.set<double>("UniBo Hole d2", 1.4829,    "d2_h:[V/cm]");

  //---------------------------------------------------------------------------
  // New University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBoNew Electron a0",  4.65403,      "a0_e:[V]"   );
  p.set<double>("UniBoNew Electron a1",  -8.76031e-3,  "a1_e:[V]"   );
  p.set<double>("UniBoNew Electron a2",  1.34037e-5,   "a2_e:[V]"   );
  p.set<double>("UniBoNew Electron a3",  -2.75108e-9,  "a3_e:[V]"   );
  p.set<double>("UniBoNew Electron b0",  -0.128302,    "b0_e:[V]"   );
  p.set<double>("UniBoNew Electron b1",  4.45552e-3,   "b1_e:[V]"   );
  p.set<double>("UniBoNew Electron b2",  -1.0866e-5,   "b2_e:[V]"   );
  p.set<double>("UniBoNew Electron b3",  9.23119e-9,   "b3_e:[V]"   );
  p.set<double>("UniBoNew Electron b4",  -1.82482e-12, "b4_e:[V]"   );
  p.set<double>("UniBoNew Electron b5",  -4.82689e-15, "b5_e:[V]"   );
  p.set<double>("UniBoNew Electron b6",  1.09402e-17,  "b6_e:[V]"   );
  p.set<double>("UniBoNew Electron b7",  -1.24961e-20, "b7_e:[V]"   );
  p.set<double>("UniBoNew Electron b8",  7.55584e-24,  "b8_e:[V]"   );
  p.set<double>("UniBoNew Electron b9",  -2.28615e-27, "b9_e:[V]"   );
  p.set<double>("UniBoNew Electron b10", 2.73344e-31,  "b10_e:[V]"  );
  p.set<double>("UniBoNew Electron c0",  7.76221e3,    "c0_e:[V/cm]");
  p.set<double>("UniBoNew Electron c1",  25.18888,     "c1_e:[V/cm]");
  p.set<double>("UniBoNew Electron c2",  -1.37417e-3,  "c2_e:[V/cm]");
  p.set<double>("UniBoNew Electron c3",  1.59525e-4,   "c3_e:[V/cm]");
  p.set<double>("UniBoNew Electron d0",  7.10481e5,    "d0_e:[V/cm]");
  p.set<double>("UniBoNew Electron d1",  3.98594e3,    "d1_e:[V/cm]");
  p.set<double>("UniBoNew Electron d2",  -7.19956,     "d2_e:[V/cm]");
  p.set<double>("UniBoNew Electron d3",  6.96431e-3,   "d3_e:[V/cm]");

  // Holes
  p.set<double>("UniBoNew Hole a0",  2.26018,      "a0_h:[V]"   );
  p.set<double>("UniBoNew Hole a1",  0.0134001,    "a1_h:[V]"   );
  p.set<double>("UniBoNew Hole a2",  -5.87724e-6,  "a2_h:[V]"   );
  p.set<double>("UniBoNew Hole a3",  -1.14021e-9,  "a3_h:[V]"   );
  p.set<double>("UniBoNew Hole b0",  0.058547,     "b0_h:[V]"   );
  p.set<double>("UniBoNew Hole b1",  -1.95755e-4,  "b1_h:[V]"   );
  p.set<double>("UniBoNew Hole b2",  2.44357e-7,   "b2_h:[V]"   );
  p.set<double>("UniBoNew Hole b3",  -1.33202e-10, "b3_h:[V]"   );
  p.set<double>("UniBoNew Hole b4",  2.68082e-14,  "b4_h:[V]"   );
  p.set<double>("UniBoNew Hole b5",  0.0,          "b5_h:[V]"   );
  p.set<double>("UniBoNew Hole b6",  0.0,          "b6_h:[V]"   );
  p.set<double>("UniBoNew Hole b7",  0.0,          "b7_h:[V]"   );
  p.set<double>("UniBoNew Hole b8",  0.0,          "b8_h:[V]"   );
  p.set<double>("UniBoNew Hole b9",  0.0,          "b9_h:[V]"   );
  p.set<double>("UniBoNew Hole b10", 0.0,          "b10_h:[V]"  );
  p.set<double>("UniBoNew Hole c0",  1.95399e4,    "c0_h:[V/cm]");
  p.set<double>("UniBoNew Hole c1",  -104.441,     "c1_h:[V/cm]");
  p.set<double>("UniBoNew Hole c2",  0.498768,     "c2_h:[V/cm]");
  p.set<double>("UniBoNew Hole c3",  0.0,          "c3_h:[V/cm]");
  p.set<double>("UniBoNew Hole d0",  2.07712e6,    "d0_h:[V/cm]");
  p.set<double>("UniBoNew Hole d1",  993.153,      "d1_h:[V/cm]");
  p.set<double>("UniBoNew Hole d2",  7.77769,      "d2_h:[V/cm]");
  p.set<double>("UniBoNew Hole d3",  0.0,          "d3_h:[V/cm]");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  // Parameters for the "PowerLawTempDep" heat capacity model
  // Linear interpolation b.t.w. InP and GaP for x = 0.49
  p.set<double>("Mass Density",        4.4807,  "rho:[g/cm^3]"   );
  p.set<double>("Heat Capacity c300",  0.4634,  "c300:[J/(K.g)]" );
  p.set<double>("Heat Capacity c1",    0.0500,  "c1:[J/(K.g)]"   );
  p.set<double>("Heat Capacity beta",  2.3195,  "beta: [1]"      );

  // Parameters for the "PowerLawTempDep" thermal conductivity model -----------
  // Linear interpolation b.t.w. InP and GaP for x = 0.49
  p.set<double>("Thermal Conductivity kappa300", 0.72,  "kappa300:[W/(K.cm)]" );
  p.set<double>("Thermal Conductivity alpha",    -1.4,  "alpha:[1]"           );

} // end of setInGaPParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setAlNParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setAlNParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          8.5,        "[1]"         );
  p.set<double>("Electron Affinity",              0.6,        "[eV]"        ); 
  p.set<double>("Band Gap",                       6.28,       "[eV]"        );
  p.set<double>("Intrinsic Concentration",        9.4e-34,    "[cm^-3]"     );
  p.set<double>("Electron Mobility",              300.0,      "[cm^2/(V.s)]");
  p.set<double>("Hole Mobility",                  14.0,       "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 7.0,        "[cm^2/s]"    );
  p.set<double>("Hole Diffusion Coefficient",     0.3,        "[cm^2/s]"    );
  p.set<double>(
    "Radiative Recombination Coefficient",
    0.4e-10,
    "Coefficient:[cm^3.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 0.6,        "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          6.28,       "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             1.799e-3,  "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              1462,       "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K",  6.3e18, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent",  1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",      4.8e20, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",      1.5,    "Nv_F:[1]"     );
  //---------------------------------------------------------------------------
  // Polarization Parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Lattice Constant",          3.112,  "A0:[Å]"     );
  p.set<double>("Elastic Constant 33",       373.0,  "C33:[GPa]"  );
  p.set<double>("Elastic Constant 13",       108.0,  "C13:[GPa]"  );
  p.set<double>("Spontaneous Polarization",  -0.081, "Psp:[C/m^2]");
  p.set<double>("Piezoelectric Constant 33", 1.46,   "E33:[C/m^2]");
  p.set<double>("Piezoelectric Constant 31", -0.60,  "E31:[C/m^2]");

} // end of setAlNParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setAlGaNParameters()
//  Al(x)Ga(1-x)N, where x=0.3
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setAlGaNParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          8.78,        "[1]"         );
  p.set<double>("Electron Affinity",              2.3,         "[eV]"        ); 
  p.set<double>("Band Gap",                       3.984,       "[eV]"        );
  p.set<double>("Intrinsic Concentration",        9.4e-34,     "[cm^-3]"     );
  p.set<double>("Electron Mobility",              790,         "[cm^2/(V.s)]");
  p.set<double>("Hole Mobility",                  144.2,       "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 19.6,        "[cm^2/s]"    );
  p.set<double>("Hole Diffusion Coefficient",     3.59,        "[cm^2/s]"    );

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 2.3,        "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          6.28,       "Eg300:[eV]"  );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K",  4.02e18,"Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent",  1.5,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",      2.66e19,"Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",      1.5,    "Nv_F:[1]"     );

} // end of setAlGaNParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setGaNParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setGaNParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          8.9,       "[1]"         );    // The following are from the NSM archive for
  p.set<double>("Electron Affinity",              4.1,       "[eV]"        );    // Wurzite.
  p.set<double>("Band Gap",                       3.39,      "[eV]"        );
  p.set<double>("Intrinsic Concentration",        3.327e-10, "[cm^-3]"     );
  p.set<double>("Electron Mobility",              1000.0,    "[cm^2/(V.s)]");
  p.set<double>("Hole Mobility",                  200.0,     "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 25.0,      "[cm^2/s]"    );
  p.set<double>("Hole Diffusion Coefficient",     5.0,       "[cm^2/s]"    );
  /*
  p.set<double>("Relative Permittivity",          9.7,       "[1]"         );    // The entries that are commented out are from
  p.set<double>("Electron Affinity",              4.1,       "[eV]"        );    // the NSM archive for Zinc Blende.
  p.set<double>("Band Gap",                       3.2,       "[eV]"        );
  p.set<double>("Intrinsic Concentration",        9.1e-9,    "[cm^-3]"     );
  p.set<double>("Electron Mobility",              1000.0,    "[cm^2/(V.s)]");
  p.set<double>("Hole Mobility",                  350.0,     "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 25.0,      "[cm^2/s]"    );
  p.set<double>("Hole Diffusion Coefficient",     9.0,       "[cm^2/s]"    );
  */
  p.set<double>("Electron Lifetime",              2e-9,      "Tau0:[s]"    );
  p.set<double>("Hole Lifetime",                  6.5e-9,    "Tau0:[s]"    );
  p.set<double>(
    "Radiative Recombination Coefficient",
    1.1e-8,
    "Coefficient:[cm^3.s^-1]");
  p.set<double>(
    "Electron Auger Coefficient",
    2.0e-31,
    "Electron Auger Coefficient:[cm^6.s^-1]");
  p.set<double>(
    "Hole Auger Coefficient",
    2.0e-31,
    "Hole Auger Coefficient:[cm^6.s^-1]");

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.1,    "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          3.39,   "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             7.7e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              600,    "beta:[K]"    );
  /*
  p.set<double>("Electron Affinity at 300 K", 4.1,    "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          3.2,    "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             7.7e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              600,    "beta:[K]"    );
  */

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 2.234e18, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,      "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     4.625e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,      "Nv_F:[1]"     );
  /*
  p.set<double>("Electron Effective DOS at 300 K", 1.2e18,   "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1.5,      "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     4.1e19,   "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1.5,      "Nv_F:[1]"     );
  */

  //---------------------------------------------------------------------------
  // Polarization Parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Lattice Constant",          3.189,  "A0:[Å]"     );
  p.set<double>("Elastic Constant 33",       405.0,  "C33:[GPa]"  );
  p.set<double>("Elastic Constant 13",       103.0,  "C13:[GPa]"  );
  p.set<double>("Spontaneous Polarization",  -0.029, "Psp:[C/m^2]");
  p.set<double>("Piezoelectric Constant 33", 0.73,   "E33:[C/m^2]");
  p.set<double>("Piezoelectric Constant 31", -0.49,  "E31:[C/m^2]");

  //---------------------------------------------------------------------------
  // Albrecht mobility model parameters (low-field).
  //---------------------------------------------------------------------------
  p.set<double>("Albrecht AN",  2.61e-4,  "AN:[(V.s)/cm^2]");
  p.set<double>("Albrecht BN",  2.9e-4,   "BN:[(V.s)/cm^2]");
  p.set<double>("Albrecht CN",  170.0e-4, "CN:[(V.s)/cm^2]");
  p.set<double>("Albrecht N0N", 1.0e17,   "N0N:[cm^-3]"    );
  p.set<double>("Albrecht T0N", 300.0,    "T0N:[K]"        );
  p.set<double>("Albrecht T1N", 1065.0,   "T1N:[K]"        );
  p.set<double>("Albrecht AP",  2.61e-4,  "AP:[(V.s)/cm^2]");
  p.set<double>("Albrecht BP",  2.9e-4,   "BP:[(V.s)/cm^2]");
  p.set<double>("Albrecht CP",  170.0e-4, "CP:[(V.s)/cm^2]");
  p.set<double>("Albrecht N0P", 1.0e17,   "N0P:[cm^-3]"    );
  p.set<double>("Albrecht T0P", 300.0,    "T0P:[K]"        );
  p.set<double>("Albrecht T1P", 1065.0,   "T1P:[K]"        );

  //---------------------------------------------------------------------------
  // Farahmand Modified Caughey Thomas (low-field).
  //---------------------------------------------------------------------------
  p.set<double>("Farahmand MU1N",   295.0,  "MU1N:[(cm^2/V.s)]");
  p.set<double>("Farahmand MU1P",   295.0,  "MU1P:[(cm^2/V.s)]");
  p.set<double>("Farahmand MU2N",   1460.7, "MU2N:[(cm^2/V.s)]");
  p.set<double>("Farahmand MU2P",   1460.7, "MU2P:[(cm^2/V.s)]");
  p.set<double>("Farahmand ALPHAN", 0.66,   "ALPHAN:[1]"       );
  p.set<double>("Farahmand ALPHAP", 0.66,   "ALPHAP:[1]"       );
  p.set<double>("Farahmand BETAN",  -1.02,  "BETAN:[1]"        );
  p.set<double>("Farahmand BETAP",  -1.02,  "BETAP:[1]"        );
  p.set<double>("Farahmand DELTAN", -3.84,  "DELTAN:[1]"       );
  p.set<double>("Farahmand DELTAP", -3.84,  "DELTAP:[1]"       );
  p.set<double>("Farahmand GAMMAN", 3.02,   "GAMMAN:[1]"       );
  p.set<double>("Farahmand GAMMAP", 3.02,   "GAMMAP:[1]"       );
  p.set<double>("Farahmand EPSN",   0.81,   "EPSN:[1]"         );
  p.set<double>("Farahmand EPSP",   0.81,   "EPSP:[1]"         );
  p.set<double>("Farahmand NCRITN", 1.0e17, "NCRITN:[cm^-3]"   );
  p.set<double>("Farahmand NCRITP", 1.0e17, "NCRITP:[cm^-3]"   );

  //---------------------------------------------------------------------------
  // High field mobility parameters.
  //---------------------------------------------------------------------------
  p.set<double>("N1N GANSAT", 7.2044,     "[1]"         );
  p.set<double>("N1P GANSAT", 1.0,        "[1]"         );
  p.set<double>("N2N GANSAT", 0.7857,     "[1]"         );
  p.set<double>("N2P GANSAT", 1.0,        "[1]"         );
  p.set<double>("ANN GANSAT", 6.1973,     "[1]"         );
  p.set<double>("ANP GANSAT", 1.0,        "[1]"         );
  p.set<double>("ECN GANSAT", 220.8936e3, "[V/cm]"      );
  p.set<double>("ECP GANSAT", 1.0,        "[V/cm]"      );
  p.set<double>("VSATN",      1.9064e7,   "VSATN:[cm/s]");
  p.set<double>("VSATP",      1.9064e7,   "VSATN:[cm/s]");

  //---------------------------------------------------------------------------
  // Selberherr avalanche model
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>(
    "Selberherr Electron a0",
    3.1e7,
    "a0_e:[1/cm] = N_IONIZATION in Charon 1.0");
  p.set<double>(
    "Selberherr Electron a1",
    0.0,
    "a1_e:[1/(cm.K)] = N_ION_1 in Charon 1.0");
  p.set<double>(
    "Selberherr Electron a2",
    0.0,
    "a2_e:[1/(cm.K^2)] = N_ION_2 in Charon 1.0");
  p.set<double>(
    "Selberherr Electron delta",
    1.0,
    "delta_e:[1] = EXN_II in Charon 1.0");
  p.set<double>(
    "Selberherr Electron E0",
    3.5e7,
    "E0_e:[V/cm] = ECN.II in Selberherr, not in Charon 1.0");
  p.set<double>(
    "Selberherr Electron lambda300",
    10.4542e-7,
    "lambda300_e:[cm] = LAN300 in Charon 1.0");
  p.set<double>(
    "Selberherr Electron hbarOmega",
    0.063,
    "hbarOmega_e:[eV] = OP_PH_EN in Charon 1.0");

  // Holes
  p.set<double>(
    "Selberherr Hole a0",
    3.1e7,
    "a0_h:[1/cm] = P_IONIZATION in Charon 1.0");
  p.set<double>(
    "Selberherr Hole a1",
    0.0,
    "a1_h:[1/(cm.K)] = P_ION_1 in Charon 1.0");
  p.set<double>(
    "Selberherr Hole a2",
    0.0,
    "a2_h:[1/(cm.K^2)] = P_ION_2 in Charon 1.0");
  p.set<double>(
    "Selberherr Hole delta",
    1.0,
    "delta_h:[1] = EXP_II in Charon 1.0");
  p.set<double>(
    "Selberherr Hole E0",
    3.5e7,
    "E0_h:[V/cm] = ECP.II in Selberherr, not in Charon 1.0");
  p.set<double>(
    "Selberherr Hole lambda300",
    6.32079e-7,
    "lambda300_h:[cm] = LAP300 in Charon 1.0");
  p.set<double>(
    "Selberherr Hole hbarOmega",
    0.063,
    "hbarOmega_h:[eV] = OP_PH_EN in Charon 1.0");

  //---------------------------------------------------------------------------
  // vanOverstraeten avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("vanOverstraeten Electron alow",  7.03e5,  "al_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron blow",  1.231e6, "bl_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron ahigh", 7.03e5,  "ah_e:[1/cm]");
  p.set<double>("vanOverstraeten Electron bhigh", 1.231e6, "bh_e:[V/cm]");
  p.set<double>("vanOverstraeten Electron E0",    4e5,     "E0_e:[V/cm]");
  p.set<double>(
    "vanOverstraeten Electron hbarOmega",
    0.063,
    "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("vanOverstraeten Hole alow",      1.582e6, "al_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole blow",      2.036e6, "bl_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole ahigh",     6.71e5,  "ah_h:[1/cm]"     );
  p.set<double>("vanOverstraeten Hole bhigh",     1.693e6, "bh_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole E0",        4e5,     "E0_h:[V/cm]"     );
  p.set<double>("vanOverstraeten Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // Okuto-Crowell impact ionization model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Okuto Electron a",     0.426,   "a_e:[1/V]"  );
  p.set<double>("Okuto Electron b",     4.81e5,  "b_e:[V/cm]" );
  p.set<double>("Okuto Electorn c",     3.05e-4, "c_e:[1/K]"  );
  p.set<double>("Okuto Electorn d",     6.86e-4, "d_e:[1/K]"  );
  p.set<double>("Okuto Electorn gamma", 1,       "gamma_e:[1]");
  p.set<double>("Okuto Electorn delta", 2,       "delta_e:[1]");

  // Holes
  p.set<double>("Okuto Hole a",     0.243,   "a_h:[1/V]"  );
  p.set<double>("Okuto Hole b",     6.53e5,  "b_h:[V/cm]" );
  p.set<double>("Okuto Hole c",     5.35e-4, "c_h:[1/K]"  );
  p.set<double>("Okuto Hole d",     5.67e-4, "d_h:[1/K]"  );
  p.set<double>("Okuto Hole gamma", 1,       "gamma_h:[1]");
  p.set<double>("Okuto Hole delta", 2,       "delta_h:[1]");

  //---------------------------------------------------------------------------
  // Lackner avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("Lackner Electron a",         1.316e6, "a_e:[1/cm]"      );
  p.set<double>("Lackner Electron b",         1.474e6, "b_e:[V/cm]"      );
  p.set<double>("Lackner Electron hbarOmega", 0.063,   "hbarOmega_e:[eV]");

  // Holes
  p.set<double>("Lackner Hole a",         1.818e6, "a_h:[1/cm]"      );
  p.set<double>("Lackner Hole b",         2.036e6, "b_h:[V/cm]"      );
  p.set<double>("Lackner Hole hbarOmega", 0.063,   "hbarOmega_h:[eV]");

  //---------------------------------------------------------------------------
  // University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBo Electron a0", 4.3383,    "a0_e:[V]"   );
  p.set<double>("UniBo Electron a1", -2.42e-12, "a1_e:[V]"   );
  p.set<double>("UniBo Electron a2", 4.1233,    "a2_e:[1]"   );
  p.set<double>("UniBo Electron b0", 0.235,     "b0_e:[V]"   );
  p.set<double>("UniBo Electron b1", 0.0,       "b1_e:[1]"   );
  p.set<double>("UniBo Electron c0", 1.6831e4,  "c0_e:[V/cm]");
  p.set<double>("UniBo Electron c1", 4.3796,    "c1_e:[V/cm]");
  p.set<double>("UniBo Electron c2", 0.13005,   "c2_e:[V/cm]");
  p.set<double>("UniBo Electron d0", 1.2337e6,  "d0_e:[V/cm]");
  p.set<double>("UniBo Electron d1", 1.2039e3,  "d1_e:[V/cm]");
  p.set<double>("UniBo Electron d2", 0.56703,   "d2_e:[V/cm]");

  // Holes
  p.set<double>("UniBo Hole a0", 2.376,     "a0_h:[V]"   );
  p.set<double>("UniBo Hole a1", 1.033e-2,  "a1_h:[V]"   );
  p.set<double>("UniBo Hole a2", 0.0,       "a2_h:[1]"   );
  p.set<double>("UniBo Hole b0", 0.17714,   "b0_h:[V]"   );
  p.set<double>("UniBo Hole b1", -2.178e-3, "b1_h:[1]"   );
  p.set<double>("UniBo Hole c0", 9.47e-3,   "c0_h:[V/cm]");
  p.set<double>("UniBo Hole c1", 2.4924,    "c1_h:[1]"   );
  p.set<double>("UniBo Hole c2", 0.0,       "c2_h:[1]"   );
  p.set<double>("UniBo Hole d0", 1.4043e6,  "d0_h:[V/cm]");
  p.set<double>("UniBo Hole d1", 2.9744e3,  "d1_h:[V/cm]");
  p.set<double>("UniBo Hole d2", 1.4829,    "d2_h:[V/cm]");

  //---------------------------------------------------------------------------
  // New University of Bologna avalanche model.
  //---------------------------------------------------------------------------

  // Electrons
  p.set<double>("UniBoNew Electron a0",  4.65403,      "a0_e:[V]"   );
  p.set<double>("UniBoNew Electron a1",  -8.76031e-3,  "a1_e:[V]"   );
  p.set<double>("UniBoNew Electron a2",  1.34037e-5,   "a2_e:[V]"   );
  p.set<double>("UniBoNew Electron a3",  -2.75108e-9,  "a3_e:[V]"   );
  p.set<double>("UniBoNew Electron b0",  -0.128302,    "b0_e:[V]"   );
  p.set<double>("UniBoNew Electron b1",  4.45552e-3,   "b1_e:[V]"   );
  p.set<double>("UniBoNew Electron b2",  -1.0866e-5,   "b2_e:[V]"   );
  p.set<double>("UniBoNew Electron b3",  9.23119e-9,   "b3_e:[V]"   );
  p.set<double>("UniBoNew Electron b4",  -1.82482e-12, "b4_e:[V]"   );
  p.set<double>("UniBoNew Electron b5",  -4.82689e-15, "b5_e:[V]"   );
  p.set<double>("UniBoNew Electron b6",  1.09402e-17,  "b6_e:[V]"   );
  p.set<double>("UniBoNew Electron b7",  -1.24961e-20, "b7_e:[V]"   );
  p.set<double>("UniBoNew Electron b8",  7.55584e-24,  "b8_e:[V]"   );
  p.set<double>("UniBoNew Electron b9",  -2.28615e-27, "b9_e:[V]"   );
  p.set<double>("UniBoNew Electron b10", 2.73344e-31,  "b10_e:[V]"  );
  p.set<double>("UniBoNew Electron c0",  7.76221e3,    "c0_e:[V/cm]");
  p.set<double>("UniBoNew Electron c1",  25.18888,     "c1_e:[V/cm]");
  p.set<double>("UniBoNew Electron c2",  -1.37417e-3,  "c2_e:[V/cm]");
  p.set<double>("UniBoNew Electron c3",  1.59525e-4,   "c3_e:[V/cm]");
  p.set<double>("UniBoNew Electron d0",  7.10481e5,    "d0_e:[V/cm]");
  p.set<double>("UniBoNew Electron d1",  3.98594e3,    "d1_e:[V/cm]");
  p.set<double>("UniBoNew Electron d2",  -7.19956,     "d2_e:[V/cm]");
  p.set<double>("UniBoNew Electron d3",  6.96431e-3,   "d3_e:[V/cm]");

  // Holes
  p.set<double>("UniBoNew Hole a0",  2.26018,      "a0_h:[V]"   );
  p.set<double>("UniBoNew Hole a1",  0.0134001,    "a1_h:[V]"   );
  p.set<double>("UniBoNew Hole a2",  -5.87724e-6,  "a2_h:[V]"   );
  p.set<double>("UniBoNew Hole a3",  -1.14021e-9,  "a3_h:[V]"   );
  p.set<double>("UniBoNew Hole b0",  0.058547,     "b0_h:[V]"   );
  p.set<double>("UniBoNew Hole b1",  -1.95755e-4,  "b1_h:[V]"   );
  p.set<double>("UniBoNew Hole b2",  2.44357e-7,   "b2_h:[V]"   );
  p.set<double>("UniBoNew Hole b3",  -1.33202e-10, "b3_h:[V]"   );
  p.set<double>("UniBoNew Hole b4",  2.68082e-14,  "b4_h:[V]"   );
  p.set<double>("UniBoNew Hole b5",  0.0,          "b5_h:[V]"   );
  p.set<double>("UniBoNew Hole b6",  0.0,          "b6_h:[V]"   );
  p.set<double>("UniBoNew Hole b7",  0.0,          "b7_h:[V]"   );
  p.set<double>("UniBoNew Hole b8",  0.0,          "b8_h:[V]"   );
  p.set<double>("UniBoNew Hole b9",  0.0,          "b9_h:[V]"   );
  p.set<double>("UniBoNew Hole b10", 0.0,          "b10_h:[V]"  );
  p.set<double>("UniBoNew Hole c0",  1.95399e4,    "c0_h:[V/cm]");
  p.set<double>("UniBoNew Hole c1",  -104.441,     "c1_h:[V/cm]");
  p.set<double>("UniBoNew Hole c2",  0.498768,     "c2_h:[V/cm]");
  p.set<double>("UniBoNew Hole c3",  0.0,          "c3_h:[V/cm]");
  p.set<double>("UniBoNew Hole d0",  2.07712e6,    "d0_h:[V/cm]");
  p.set<double>("UniBoNew Hole d1",  993.153,      "d1_h:[V/cm]");
  p.set<double>("UniBoNew Hole d2",  7.77769,      "d2_h:[V/cm]");
  p.set<double>("UniBoNew Hole d3",  0.0,          "d3_h:[V/cm]");
} // end of setGaNParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setTiO2Parameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setTiO2Parameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          10,       "[1]"         );
  p.set<double>("Electron Affinity",              4.0,      "[eV]"        );
  p.set<double>("Band Gap",                       3.0,      "[eV]"        );
  p.set<double>("Electron Mobility",              1,        "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 0.0258,   "[cm^2/s]"    );
  p.set<double>("Hole Mobility",                  0.5,      "[cm^2/(V.s)]");
  p.set<double>("Hole Diffusion Coefficient",     0.0129,   "[cm^2/s]"    );
  p.set<double>("Ion Mobility",                   1e-10,    "[cm^2/(V.s)]");
  p.set<double>("Ion Diffusion Coefficient",      2.58e-12, "[cm^2/s]"    );

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.0,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          3.0,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             4.73e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              636.0,   "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 5e20, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1,    "Nv_F:[1]"     );

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        1.63,    "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,       "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,       "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 0.03,    "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 1.56e-3, "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 1.65e-6, "c:[cm/(W.K)]"    );

  //---------------------------------------------------------------------------
  // Energy barrier for the soret coefficient.
  //---------------------------------------------------------------------------
  p.set<double>("Soret Energy Barrier", 1.2, "[eV]");

  //---------------------------------------------------------------------------
  // RigidPointIon mobility model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Ion Escape Frequency",  1e13,    "[Hz]");
  p.set<double>("Ion Hopping Distance",  0.15e-7, "[cm]");
  p.set<double>("Ion Activation Energy", 1.2,     "[eV]");                       // Updated November, 2014.
} // end of setTiO2Parameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setTantalumParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setTantalumParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          30,       "[1]"         );
  p.set<double>("Electron Affinity",              4.0,      "[eV]"        );
  p.set<double>("Band Gap",                       4.0,      "[eV]"        );
  p.set<double>("Electron Mobility",              1,        "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 0.0258,   "[cm^2/s]"    );
  p.set<double>("Hole Mobility",                  0.5,      "[cm^2/(V.s)]");
  p.set<double>("Hole Diffusion Coefficient",     0.0129,   "[cm^2/s]"    );
  p.set<double>("Ion Mobility",                   1e-10,    "[cm^2/(V.s)]");
  p.set<double>("Ion Diffusion Coefficient",      2.58e-12, "[cm^2/s]"    );
  p.set<double>("Intrinsic Concentration",        2.15e-14, "[cm^2/s]"    );
 

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.0,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          4.0,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             4.73e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              636.0,   "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 1e21, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1,    "Nv_F:[1]"     );

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        2.34, "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,    "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,    "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 1.74, "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 0,    "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 0,    "c:[cm/(W.K)]"    );

  //---------------------------------------------------------------------------
  // Energy barrier for the soret coefficient.
  //---------------------------------------------------------------------------
  p.set<double>("Soret Energy Barrier", 1.2, "[eV]");

  //---------------------------------------------------------------------------
  // RigidPointIon mobility model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Ion Escape Frequency",  1e13,    "[Hz]");
  p.set<double>("Ion Hopping Distance",  0.15e-7, "[cm]");
  p.set<double>("Ion Activation Energy", 1.2,     "[eV]");
} // end of setTantalumParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setPlatinumParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setPlatinumParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Metal", "");
  p.set<double>("Work Function", 4.0, "[eV]");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        2.84, "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,    "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,    "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 1.4,  "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 0,    "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 0,    "c:[cm/(W.K)]"    );
} // end of setPlatinumParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setPlatinumSemiParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setPlatinumSemiParameters(
  Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          30,       "[1]"         );
  p.set<double>("Electron Affinity",              4.0,      "[eV]"        );
  p.set<double>("Band Gap",                       4.0,      "[eV]"        );
  p.set<double>("Electron Mobility",              1,        "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 0.0258,   "[cm^2/s]"    );
  p.set<double>("Hole Mobility",                  0.5,      "[cm^2/(V.s)]");
  p.set<double>("Hole Diffusion Coefficient",     0.0129,   "[cm^2/s]"    );
  p.set<double>("Ion Mobility",                   1e-10,    "[cm^2/(V.s)]");
  p.set<double>("Ion Diffusion Coefficient",      2.58e-12, "[cm^2/s]"    );

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.0,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          4.0,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             4.73e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              636.0,   "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 1e21, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1,    "Nv_F:[1]"     );

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        2.84, "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,    "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,    "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 1.4,  "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 0,    "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 0,    "c:[cm/(W.K)]"    );
} // end of setPlatinumSemiParameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setTa2O5Parameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setTa2O5Parameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Insulator", "");
  p.set<double>("Relative Permittivity", 25,  "[1]" );
  p.set<double>("Electron Affinity",     4.0, "[eV]");
  p.set<double>("Band Gap",              4.0, "[eV]");

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        1.63, "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,    "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,    "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 300., "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 0,    "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 0,    "c:[cm/(W.K)]"    );
} // end of setTa2O5Parameters()

///////////////////////////////////////////////////////////////////////////////
//
//  setTaOParameters()
//
///////////////////////////////////////////////////////////////////////////////
void charon::Material_Properties::setTaOParameters(Teuchos::ParameterList& p)
{
  //---------------------------------------------------------------------------
  // General default parameters.
  //---------------------------------------------------------------------------
  p.set<std::string>("Material Type", "Semiconductor", "");
  p.set<double>("Relative Permittivity",          30,       "[1]"         );
  p.set<double>("Electron Affinity",              4.0,      "[eV]"        );
  p.set<double>("Band Gap",                       4.0,      "[eV]"        );
  p.set<double>("Electron Mobility",              1,        "[cm^2/(V.s)]");
  p.set<double>("Electron Diffusion Coefficient", 0.0258,   "[cm^2/s]"    );
  p.set<double>("Hole Mobility",                  0.5,      "[cm^2/(V.s)]");
  p.set<double>("Hole Diffusion Coefficient",     0.0129,   "[cm^2/s]"    );
  p.set<double>("Ion Mobility",                   1e-10,    "[cm^2/(V.s)]");
  p.set<double>("Ion Diffusion Coefficient",      2.58e-12, "[cm^2/s]"    );
  p.set<double>("Intrinsic Concentration",        1.52e-14, "[cm^2/s]"    );

  //---------------------------------------------------------------------------
  // Bandgap and affinity temperature-dependent model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Affinity at 300 K", 4.0,     "Chi300:[eV]" );
  p.set<double>("Band Gap at 300 K",          4.0,     "Eg300:[eV]"  );
  p.set<double>("Band Gap alpha",             4.73e-4, "alpha:[eV/K]");
  p.set<double>("Band Gap beta",              636.0,   "beta:[K]"    );

  //---------------------------------------------------------------------------
  // Intrinsic concentration models parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Electron Effective DOS at 300 K", 5e20, "Nc300:[cm^-3]");
  p.set<double>("Electron Effective DOS Exponent", 1,    "Nc_F:[1]"     );
  p.set<double>("Hole Effective DOS at 300 K",     1e19, "Nv300:[cm^-3]");
  p.set<double>("Hole Effective DOS Exponent",     1,    "Nv_F:[1]"     );

  //---------------------------------------------------------------------------
  // Heat capacity and thermal conductivity parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Heat Capacity a",        2.34, "a:[J/(K.cm^3)]"  );
  p.set<double>("Heat Capacity b",        0,    "b:[J/(K^2.cm^3)]");
  p.set<double>("Heat Capacity c",        0,    "c:[J/(K^3.cm^3)]");
  p.set<double>("Thermal Conductivity a", 1.74, "a:[cm.K/W]"      );
  p.set<double>("Thermal Conductivity b", 0,    "b:[cm/W]"        );
  p.set<double>("Thermal Conductivity c", 0,    "c:[cm/(W.K)]"    );

  //---------------------------------------------------------------------------
  // Energy barrier for the soret coefficient.
  //---------------------------------------------------------------------------
  p.set<double>("Soret Energy Barrier", 1.2, "[eV]");

  //---------------------------------------------------------------------------
  // RigidPointIon mobility model parameters.
  //---------------------------------------------------------------------------
  p.set<double>("Ion Escape Frequency",  1e13,    "[Hz]");
  p.set<double>("Ion Hopping Distance",  0.15e-7, "[cm]");
  p.set<double>("Ion Activation Energy", 1.2,     "[eV]");
} // end of setTaOParameters()

// *****************************************************************************
// method to obtain singleton instance of Material_Properties
// *****************************************************************************
charon::Material_Properties& charon::Material_Properties::getInstance ()
{
   static Material_Properties instance;
   return instance;
}


// *****************************************************************************
// validate the material name, if not found, throw Teuchos::Exceptions
// *****************************************************************************
void charon::Material_Properties::validateMaterialName(const std::string& materialName)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isSublist(materialName),
       Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid material name "
                     << materialName <<  std::endl);
}


// *****************************************************************************
// set the property value for a given material and property name
// *****************************************************************************
void charon::Material_Properties::setPropertyValue
(const std::string& materialName, const std::string& propertyName, double propertyValue)
{
  // check if the material name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isSublist(materialName),
         Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid material name "
                     << materialName <<  std::endl);

  Teuchos::ParameterList& matList = pMaterials.sublist(materialName);

  // check if the property name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!matList.isParameter(propertyName),
                     Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid property name "
                     << propertyName <<  std::endl);

  // if the material and property names both exist, set the value
  matList.set<double>(propertyName, propertyValue);
}


// *****************************************************************************
// set the property value for a material-independent property name
// *****************************************************************************
void charon::Material_Properties::setPropertyValue
(const std::string& propertyName, double propertyValue)
{
  // check if the property name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isParameter(propertyName),
                     Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid property name "
                     << propertyName <<  std::endl);

  // if the property name exists, set the value
  pMaterials.set<double>(propertyName, propertyValue);
}


// *****************************************************************************
// get the property value for a given material and property name
// *****************************************************************************
double charon::Material_Properties::getPropertyValue
(const std::string& materialName, const std::string& propertyName)
{
  // check if the material name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isSublist(materialName),
                     Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid material name "
                     << materialName <<  std::endl);

  Teuchos::ParameterList& matList = pMaterials.sublist(materialName);

  // check if the property name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!matList.isParameter(propertyName),
                     Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! " << materialName << " does not have the property name of "
                     << propertyName <<  std::endl);

  // if the material and property names both exist, return the value
  return matList.get<double>(propertyName);
}


// *****************************************************************************
// get the type for a given material name
// *****************************************************************************
std::string charon::Material_Properties::getMaterialType(const std::string& materialName)
{
  // check if the material name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isSublist(materialName),
               Teuchos::Exceptions::InvalidParameter, std::endl
               << "Material_Properties Error! Invalid material name "
               << materialName <<  std::endl);

  Teuchos::ParameterList& matList = pMaterials.sublist(materialName);

  return matList.get<std::string>("Material Type");
}


// *****************************************************************************
// get the property value for a material-independent property name
// *****************************************************************************
double charon::Material_Properties::getPropertyValue(const std::string& propertyName)
{
  // check if the property name exists
  TEUCHOS_TEST_FOR_EXCEPTION(!pMaterials.isParameter(propertyName),
                     Teuchos::Exceptions::InvalidParameter, std::endl
                     << "Material_Properties Error! Invalid property name "
                     << propertyName <<  std::endl);

  // if the material-independent property exists, return the value
  return pMaterials.get<double>(propertyName);
}
