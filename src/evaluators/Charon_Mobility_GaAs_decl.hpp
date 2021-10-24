
#ifndef CHARON_MOBILITY_GAAS_DECL_HPP
#define CHARON_MOBILITY_GAAS_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

/**
 * @brief The low-field mobility for electron and hole can be a constant default value from
 * Material_Properties, or a user-defined constant value, or reads from a file.
 *
 * When reading from a file, the code assumes the file is a two-column space-separated
 * text file with the first column being the doping in [1/cm^3] and the second column
 * being the mobility in [cm^2/(V.s)]. The code first sorts the table according to
 * doping values in an ascending order and then linearly interpolates an interval
 * in the table for the doping at a given spatial location to obtain the mobility.
 *
 * The temperature dependence is included in the low field mobility in addition to the doping
 * dependence using lfMob(TL) = lfMob(doping)*(300/TL)^tempExp.
 *
 * tempExp is given by specifying "Temperature Exponent" in the input, and it is set to 0 when
 * not specified by the user.
 *
 * The high-field dependence of electron mobility for GaAs comes from J. J. Barnes,
 * R. J. Lomax, and G. I. Haddad, IEEE Transactions on Electron Devices ED-23, 1042 (1976).
 * It takes the following form:
 *   ScalarT fieldRatio = std::pow(hiField, 3.) / std::pow(Fsat, 4.);
 *   hfMob = (lfMob + vsat * fieldRatio) / (1. + hiField * fieldRatio);
 *
 * The temperature dependence in "vsat" is modeled as
 * vsat(TL) = vsat300 / (1 - vsatTempCoeff + vsatTempCoeff * TL/300.0),
 * according to http://www.iue.tuwien.ac.at/phd/quay/node39.html.
 * vsat300 and vsatTempCoeff are specified respectively by "Saturation Velocity"
 * and "Vsat Temperature Coefficient".
 *
 * The high-field dependence of hole mobility for GaAs comes from
 * D. M. Caughey and R. E. Thomas, Proceedings of the IEEE 55, 2192 (1967).
 * It takes the following form:
 *  hfMob = lfMob / (1. + lfMob * hiField / vsat);
 *
 * hiField is computed by one of the two methods:
 *
 * (1) when Driving Force = ElectricField (default), hiField = abs(effective electric field) that
 * includes BGN contribution, i.e., Fn/p,eff = grad(Ei/V0 -/+ 0.5*dEg/V0) (scaled).
 *
 * (2) when Driving Force = GradQuasiFermi, hiField = abs(gradient of quasi-fermi
 * potential) that includes both drift and diffusion contributions, i.e.,
 * grad_qfp = -grad(n)/n-Fn,eff for n and = grad(p)/p-Fp,eff for p (scaled), valid
 * for isothermal simulation.
 *
 * For the SUPG-FEM formulation, hiField lives at IPs; while for CVFEM-SG and EFFPG-FEM,
 * hiField lives at centers of primary edges.
 *
 * An example of using the model is given below:
 *
 * <ParameterList name="Electron Mobility">
 *    <Parameter name="Value" type="string" value="GaAs" />
 *    <!-- <Parameter name="Low Field Mobility Value" type="double" value="5000.0" />  !-->
 *    <Parameter name="Low Field Mobility File" type="string" value="GaAs_LowField_EMob.txt" />
 *    <Parameter name="High Field" type="string" value="On" />
 *    <Parameter name="Driving Force" type="string" value="ElectricField" />
 *    <Parameter name="Saturation Velocity" type="double" value="0.85e7" />
 *    <Parameter name="Saturation Field" type="double" value="4e3" />
 *    <Parameter name="Temperature Exponent" type="double" value="2.0" />
 *    <Parameter name="Vsat Temperature Coefficient" type="double" value="0.56" />
 * </ParameterList>
 *
 * <ParameterList name="Hole Mobility">
 *    <Parameter name="Value" type="string" value="GaAs" />
 *    <!-- <Parameter name="Low Field Mobility Value" type="double" value="500.0" />  !-->
 *    <Parameter name="Low Field Mobility File" type="string" value="GaAs_LowField_HMob.txt" />
 *    <Parameter name="High Field" type="string" value="On" />
 *    <Parameter name="Driving Force" type="string" value="ElectricField" />
 *    <Parameter name="Saturation Velocity" type="double" value="0.85e7" />
 *    <Parameter name="Saturation Field" type="double" value="4e3" />
 * </ParameterList>
 */

template<typename EvalT, typename Traits>
class Mobility_GaAs
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_GaAs(
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

  // initialize mobility parameters
  void initMobilityParams(const std::string& matName, const Teuchos::ParameterList& mobParamList);

  // evaluate the low field mobility
  ScalarT evaluateLowFieldMobility(const ScalarT& Na, const ScalarT& Nd);

  // compute GaAs mobility for CVFEM-SG and EFFPG-FEM at primary edge center
  ScalarT evaluateMobilityForEdgedl(const std::size_t& cell, const int& edge, const ScalarT& edgelfMob,
    const Kokkos::DynRankView<double,PHX::Device>& edgePoints, const ScalarT& edgeLatt);

  // compute GaAs mobility for SUPG-FEM at IP
  ScalarT evaluateMobilityForIPdl(const std::size_t& cell, const int& point,
    const ScalarT& lfMob, const ScalarT& latt);

  // read the low field mobility from a file
  void readLowFieldMobilityFile(const Teuchos::ParameterList& mobParamList);

  // output
  PHX::MDField<ScalarT,Cell,Point> mobility; // @ IPs or Edge

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp; // scaled

  PHX::MDField<const ScalarT,Cell,Point> acceptor;  // scaled
  PHX::MDField<const ScalarT,Cell,Point> donor;     // scaled
  PHX::MDField<const ScalarT,Cell,Point> edensity;  // scaled
  PHX::MDField<const ScalarT,Cell,Point> hdensity;  // scaled

  PHX::MDField<const ScalarT,Cell,Point> intrin_fermi; // intrinsic fermi energy in [eV]
  PHX::MDField<const ScalarT,Cell,Point> bandgap;      // band gap w/o BGN in [eV]
  PHX::MDField<const ScalarT,Cell,Point> eff_bandgap;  // effective band gap w/ BGN in [eV]

  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_qfp;  // gradient of quasi-fermi potential at IP, scaled
  PHX::MDField<const ScalarT,Cell,Point,Dim> eff_field; // effective electric field at IP, scaled

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double Mu0; // mobility scaling, [cm^3/(V.s)]
  double C0;  // conc. scaling, [cm^-3]
  double X0;  // length scaling, [cm]
  double E0;  // electric field scaling, [V/cm]
  double T0;  // temperature scaling, [K]

  int num_ips;
  int num_nodes;
  int num_points;
  int num_dims;
  int num_edges;

  std::string basis_name;
  std::size_t basis_index;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  // carrier type
  std::string carrType;

  // driving force name for high field
  std::string driveForce;

  double sign;   // + for n and - for p

  // saturation velocity in [cm/s]
  double vsat300;

  // saturation field in [V/cm]
  double Fsat;

  // low field mobility in [cm^2/(V.s)]
  double lowMob;

  // turn on/off high field dependence
  bool hiFieldOn;

  // want mobility at the edge data layout ?
  bool isEdgedl;

  // read low field mobility from a file
  bool isMobFromFile;

  // temperature exponent for the low field mobility
  double tempExponent;    // unitless

  // temperature coefficient for the saturation velocity
  double vsatTempCoeff;   // unitless

  // use a struct to sort the doping-mobility pair according to doping values in ascending order
  struct dopMobStruct
  {
    double dop, mob;

    inline bool operator < (const dopMobStruct &dms) const
    { return ( dop < dms.dop); }

    inline bool operator == (const dopMobStruct &dms) const
    { return ( dop == dms.dop); }
  };

  std::map<ScalarT,ScalarT> dopMobMap;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_GaAs


}

#endif
