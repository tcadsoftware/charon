
#ifndef CHARON_MOBILITY_MOSFET_DECL_HPP
#define CHARON_MOBILITY_MOSFET_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;
using panzer::Edge;

namespace charon {
 

/**
 * @brief 
 *
 * An example of using the model is given below:
 *
 * <ParameterList name="Electron Mobility">
 *    <Parameter name="Value" type="string" value="MOSFET" />
 *    <!-- <Parameter name="Low Field Mobility Value" type="double" value="5000.0" />  !-->
 *    <Parameter name="Low Field Mobility File" type="string" value="MOSFET_LowField_EMob.txt" />
 *    <Parameter name="High Field" type="string" value="On" />
 *    <Parameter name="Driving Force" type="string" value="ElectricField" />
 *    <Parameter name="Saturation Velocity" type="double" value="0.85e7" />
 *    <Parameter name="Saturation Field" type="double" value="4e3" />
 *    <Parameter name="Temperature Exponent" type="double" value="2.0" />
 *    <Parameter name="Vsat Temperature Coefficient" type="double" value="0.56" />
 * </ParameterList>
 *
 * <ParameterList name="Hole Mobility">
 *    <Parameter name="Value" type="string" value="MOSFET" />
 *    <!-- <Parameter name="Low Field Mobility Value" type="double" value="500.0" />  !-->
 *    <Parameter name="Low Field Mobility File" type="string" value="MOSFET_LowField_HMob.txt" />
 *    <Parameter name="High Field" type="string" value="On" />
 *    <Parameter name="Driving Force" type="string" value="ElectricField" />
 *    <Parameter name="Saturation Velocity" type="double" value="0.85e7" />
 *    <Parameter name="Saturation Field" type="double" value="4e3" />
 * </ParameterList>
 */

template<typename EvalT, typename Traits>
class Mobility_MOSFET
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Mobility_MOSFET(
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

  PHX::MDField<const ScalarT,Cell,Point> bulk_mobility; // bulk low field mobility
  PHX::MDField<const ScalarT,Cell,Point> perp_mobility; // bulk low field mobility

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


  //The bulk mobility model
  std::string bulkMobilityModel;

  //The perpendicular Mobility Model
  std::string perpMobilityModel;

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

  //Homotopy/scaling stuff
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > mobilityScaling;


  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class Mobility_MOSFET


}

#endif
