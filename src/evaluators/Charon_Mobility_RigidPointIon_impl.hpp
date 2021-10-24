
#ifndef CHARON_MOBILITY_RIGIDPOINTION_IMPL_HPP
#define CHARON_MOBILITY_RIGIDPOINTION_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
The RigidPointIon model comes from the paper by D. B. Strukov and R. S. Williams,
Appl. Phys. A 94, 515-519 (2009) and references in the paper. According to the
original paper by N. Cabrera and N. F. Mott, Rep. Prog. Phys. 12, 163 (1948),
eq.(2) in Strukov's paper is missing a factor of 2. Hence,

v = 2*f*a*exp(-Ua/kbT)*sinh(qEa/(2kbT)) = mu*Ec*sinh(E/Ec), where
Ec = 2*kbT/(q*a), and mu = q*f*a*a/kbT * exp(-Ua/kbT).

Then the diffusion coefficient D is given by D = kbT/q*mu = f*a*a*exp(-Ua/kbT).
The D expression is consistent with D. B. Strukov, F. Alibart, and R. S. Williams,
Appl. Phys. A. 107, 509-518 (2012).

f = Ion Escape Frequency in [Hz], a = Ion Hopping Distance in [cm], and
Ua = Ion Activation Energy in [eV].

Note that the paper by Sungho Kim, ShinHyun Choi, and Wei Lu, ACSNano. 8, 2369-2376
(2014) is missing a factor of 2 for v, and has an extra factor of 1/2 for D.

<ParameterList name="Ion Mobility">
     <Parameter name="Value" type="string" value="RigidPointIon"/>
     <Parameter name="Escape Frequency" type="double" value="1e12" />
     <Parameter name="Hopping Distance" type="double" value="0.1e-7" />
     <Parameter name="Activation Energy" type="double" value="0.8" />
</ParameterList>

When Escape Frequency, Hopping Distance, or Activation Energy is not specified,
the code will look for their values from Material_Properties.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_RigidPointIon<EvalT, Traits>::
Mobility_RigidPointIon(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  // Obtain the sign of the ion charge
  int ion_charge = p.get<int>("Ion Charge");
  if (ion_charge > 0)
    sign = 1.0;   // positive charge
  else if (ion_charge < 0)
    sign = -1.0;  // negative charge
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
      << "Error: Ion Charge cannot be zero !");

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize the mobility parameters
  initMobilityParams(matName, mobParamList);

  // Set up data layout
  RCP<DataLayout> input_scalar = scalar;
  RCP<DataLayout> input_vector = vector;
  RCP<DataLayout> output_scalar = scalar;
  RCP<DataLayout> output_vector = vector;

  // retrieve edge data layout if isEdgedl = true
  isEdgedl = p.get<bool>("Is Edge Data Layout");
  if (isEdgedl)
  {
    RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
    basis_name = basis->name();
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    input_scalar = basis->functional;
    output_scalar = cellTopoInfo->edge_scalar;
    num_edges = output_scalar->dimension(1);
    cellType = cellTopoInfo->getCellTopology();
  }

  // Evaluated fields
  mobility = MDField<ScalarT,Cell,Point>(n.field.ion_mobility, output_scalar);
  this->addEvaluatedField(mobility);

  if (isEdgedl)  // edge velocity (scalar)
  {
    edge_velocity = MDField<ScalarT,Cell,Point>(n.field.ion_velocity, output_scalar);
    this->addEvaluatedField(edge_velocity);
  }
  else  // velocity at IPs (vector)
  {
    ip_velocity = MDField<ScalarT,Cell,Point,Dim>(n.field.ion_velocity, output_vector);
    this->addEvaluatedField(ip_velocity);
  }

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  E0 = scaleParams->scale_params.E0;
  T0 = scaleParams->scale_params.T0;
  C0 = scaleParams->scale_params.C0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,input_scalar);
  carr_dens = MDField<const ScalarT,Cell,Point>(n.dof.iondensity,input_scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(carr_dens);

  if (isEdgedl)
  {
    potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,input_scalar);
    this->addDependentField(potential);
  }
  else
  {
    ip_ionfield = MDField<const ScalarT,Cell,Point,Dim>(n.field.ion_efield,input_vector);
    this->addDependentField(ip_ionfield);
  }

  std::string name = "RigidPointIon_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_RigidPointIon<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if (isEdgedl)
    basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_RigidPointIon<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;       // Boltzmann constant in [eV/K]

 //----------------------------------------------------------------------------
 // compute the ion mobility and velocity at the center of a primary edge
 //----------------------------------------------------------------------------
 if (isEdgedl)
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int edge = 0; edge < num_edges; ++edge)
    {
      // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
      int node0 = cellType->getNodeMap(1,edge,0);
      int node1 = cellType->getNodeMap(1,edge,1);

      // obtain temperature in [K] at the center of a primary edge
      ScalarT lattT = (latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0 * T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      ScalarT kbT = lattT*kb;  // [eV]

      // get ion density in [cm^-3] at the center of a primary edge
      ScalarT ionDens = (carr_dens(cell,node0) + carr_dens(cell,node1)) / 2.0 * C0;

      // compute ion mobility
      ScalarT ionMob = computeIonMobility(kbT, ionDens);

      // assign the value to edge ion mobility
      mobility(cell,edge) = ionMob / Mu0;  // scaled

      // get local node coordinates
      double x0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,0);
      double x1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,0);
      double y0 = 0.0, y1 = 0.0;
      double z0 = 0.0, z1 = 0.0;
      if (num_dims > 1)  // 2D or 3D
      {
        y0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,1);
        y1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,1);
      }
      if (num_dims > 2)  // 3D
      {
        z0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,2);
        z1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,2);
      }

      // compute the primary cell edge length [scaled]
      double edgeLen = std::sqrt(pow((x0-x1),2.0) + pow((y0-y1),2.0) + pow((z0-z1),2.0) );

      // compute the primary cell edge field in [V/cm]
      ScalarT ionF = (potential(cell,node0) - potential(cell,node1)) / edgeLen * E0;

      // compute the ion velocity at a primary cell edge in [cm/s]
      ScalarT ionVel = computeIonVelocity(kbT, ionMob, ionF);

      // assign the value to edge ion velocity
      edge_velocity(cell,edge) = ionVel / (Mu0 * E0);  // scaled

    }  // end of loop over edges
  }  // end of loop over cells

 }  // end of if (isEdgedl) block


 //----------------------------------------------------------------------------
 // Compute the ion mobility and velocity at IP
 //----------------------------------------------------------------------------
 else
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature in [K]
      ScalarT lattT = latt_temp(cell,point) * T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;
      ScalarT kbT = lattT*kb;  // [eV]

      // get ion density in [cm^-3]
      ScalarT ionDens = carr_dens(cell,point) * C0;

      // compute ion mobility
      ScalarT ionMob = computeIonMobility(kbT, ionDens);  // [cm^2/(V.s)]

      // assign the value to ion mobility
      mobility(cell,point) = ionMob / Mu0;  // scaled

      // compute ion velocity
      for (int dim = 0; dim < num_dims; ++dim)
      {
        ScalarT ionF = ip_ionfield(cell,point,dim) * E0;  // [V/cm]
        ScalarT ionVel = computeIonVelocity(kbT,ionMob,ionF);  // [cm/s]
        ip_velocity(cell,point,dim) = ionVel / (Mu0 * E0);  // scaled
      }

    }  // end of loop over ip points
  }  // end of loop over cells

 }  // end of else (isEdgedl) block

}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_RigidPointIon<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for the mobility model
  if (mobParamList.isParameter("Escape Frequency"))
    escFreq = mobParamList.get<double>("Escape Frequency");
  else
    escFreq = matProperty.getPropertyValue(matName, "Ion Escape Frequency");

  if (mobParamList.isParameter("Hopping Distance"))
    hopDist = mobParamList.get<double>("Hopping Distance");
  else
    hopDist = matProperty.getPropertyValue(matName, "Ion Hopping Distance");

  if (mobParamList.isParameter("Activation Energy"))
    actE = mobParamList.get<double>("Activation Energy");
  else
    actE = matProperty.getPropertyValue(matName, "Ion Activation Energy");

  maxIonDens = mobParamList.get<double>("Maximum Ion Density");

  bSetMaxDens = mobParamList.get<bool>("Enforce Maximum Ion Density");

  velMultiplier = 1.0;
  if (mobParamList.isParameter("Velocity Multiplier"))
    velMultiplier = mobParamList.get<double>("Velocity Multiplier");
}


///////////////////////////////////////////////////////////////////////////////
//
//  computeIonMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits> typename Mobility_RigidPointIon<EvalT,Traits>::ScalarT
Mobility_RigidPointIon<EvalT, Traits>::computeIonMobility(const ScalarT& kbT, const ScalarT& ionDens)
{
  // compute temperature-dependent ion mobility in [cm^2/(V.s)]
  ScalarT tempIonMob = (1.0*escFreq*hopDist*hopDist/kbT)*std::exp(-actE/kbT);  // where 1.0 represents 1[e]

  // compute ion density ratio
  ScalarT iratio = ionDens / maxIonDens;

  ScalarT ionMob = 0.0;

  if (bSetMaxDens)   // apply the (1-iratio) factor
  {
    // note that ratio could be negative since iondens can be negative during solving
    if (Sacado::ScalarValue<ScalarT>::eval(iratio) <= 0.0)
      ionMob = tempIonMob;

    else if (Sacado::ScalarValue<ScalarT>::eval(iratio) >= 1.0)
      ionMob = 0.0;

    else  // 0. < ratio < 1.
      ionMob = tempIonMob * (1. - iratio);
      // ionMob = tempIonMob;
  }
  else  // do not apply the (1-iratio) factor
    ionMob = tempIonMob;

  return ionMob;   // [cm^/(V.s)]
}


///////////////////////////////////////////////////////////////////////////////
//
//  computeIonVelocity()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits> typename Mobility_RigidPointIon<EvalT,Traits>::ScalarT
Mobility_RigidPointIon<EvalT, Traits>::computeIonVelocity(const ScalarT& kbT, const ScalarT& ionMob, const ScalarT& ionF)
{
  // characteristic field in [V/cm]
  ScalarT Fcrit = 2.0 * kbT / (1.0 * hopDist);

  ScalarT ratio = std::fabs(ionF) / Fcrit;
  ScalarT ionVel = 0.0;

  // when ratio is too large, sinh(ratio) goes to infinity, leading to velocity = nan,
  // hence limit ratio to <= 200.0 for the sinh(.) calculation, sinh(200.0) = 3.6e86.
  // Two extremes for large ratio: (1) Fcrit is too small due to unphysically small positive T and
  // ionF has a reasonable value. In this case, ionMob = 0, so velocity should be 0.
  // (2) Fcrit has a reasonable value while |ionF| >> Fcrit, i.e., in the extremely strong field regime.
  // In this case, the RigidPointIon model loses its validity.

  if (Sacado::ScalarValue<ScalarT>::eval(ratio) >= 200.0)
    ionVel = ionMob * ionF * sign;  // [cm/s]
  else
    ionVel = ionMob * Fcrit * std::sinh(ionF/Fcrit) * sign;  // [cm/s]

  ionVel = ionVel * velMultiplier;

  return ionVel;  // [cm/s]
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_RigidPointIon<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<std::string>("Material Name", "?");
  p->set<int>("Ion Charge", 1);
  p->set<bool>("Is Edge Data Layout", false);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "RigidPointIon", "RigidPointIon mobility model");
  p->sublist("Mobility ParameterList").set<double>("Escape Frequency", 0., "[1/s]");
  p->sublist("Mobility ParameterList").set<double>("Hopping Distance", 0., "[cm]");
  p->sublist("Mobility ParameterList").set<double>("Activation Energy", 0., "[eV]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Ion Density", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<bool>("Enforce Maximum Ion Density", false, "Enforce Maximum Ion Density");

  p->sublist("Mobility ParameterList").set<double>("Velocity Multiplier", 0., "[1]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
