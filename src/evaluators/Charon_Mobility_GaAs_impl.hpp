
#ifndef CHARON_MOBILITY_GAAS_IMPL_HPP
#define CHARON_MOBILITY_GAAS_IMPL_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_GaAs<EvalT, Traits>::
Mobility_GaAs(
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

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // BASIS
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis_scalar->dimension(1);
  basis_name = basis->name();

  // Check if edge data layout is specified; if yes, isEdgedl = true
  isEdgedl = p.get<bool>("Is Edge Data Layout");

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize the GaAs mobility parameters, must be called here since
  // it uses the above variables in the code
  initMobilityParams(matName, mobParamList);

  // Set data layouts for the fields
  RCP<DataLayout> output_scalar = ip_scalar;  // default for SUPG-FEM
  RCP<DataLayout> input_scalar = ip_scalar;
  RCP<DataLayout> input_vector = ip_vector;
  num_points = num_ips;

  // Retrieve edge data layout and cellType if isEdgedl = true (for CVFEM-SG and EFFPG-FEM)
  if (isEdgedl)
  {
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
    RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
    num_edges = edge_scalar->dimension(1);
    cellType = cellTopoInfo->getCellTopology();

    // Overwrite the data layouts
    output_scalar = edge_scalar;
    input_scalar = basis_scalar;
    input_vector = edge_vector;
    num_points = num_nodes;
  }

  // Evaluated field
  if (carrType == "Electron")
    mobility = MDField<ScalarT,Cell,Point>(n.field.elec_mobility, output_scalar);
  else if (carrType == "Hole")
    mobility = MDField<ScalarT,Cell,Point>(n.field.hole_mobility, output_scalar);
  this->addEvaluatedField(mobility);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  C0 = scaleParams->scale_params.C0;
  X0 = scaleParams->scale_params.X0;
  E0 = scaleParams->scale_params.E0;
  T0 = scaleParams->scale_params.T0;

  // Dependent field
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,input_scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,input_scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);

  if (isEdgedl)  // for CVFEM-SG and EFFPG-FEM
  {
    // Always need these fields when isEdgedl = true
    intrin_fermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi, input_scalar);
    bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap, input_scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, input_scalar);

    this->addDependentField(intrin_fermi);
    this->addDependentField(bandgap);
    this->addDependentField(eff_bandgap);

    if ((hiFieldOn) && (driveForce == "GradQuasiFermi"))
    {
      edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
      hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);
      this->addDependentField(edensity);
      this->addDependentField(hdensity);
    }
  }  // end of if (isEdgedl)

  else  // for SUPG-FEM
  {
   if (hiFieldOn)
   {
    // driveForce == GradQuasiFermi
    if (driveForce == "GradQuasiFermi")
    {
      if (carrType == "Electron")
        grad_qfp = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_grad_qfp, input_vector);
      else if (carrType == "Hole")
        grad_qfp = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_grad_qfp, input_vector);
      this->addDependentField(grad_qfp);
    }

    // driveForce == ElectricField
    else if (driveForce == "ElectricField")
    {
      if (carrType == "Electron")
        eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield, input_vector);
      else if (carrType == "Hole")
        eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield, input_vector);
      this->addDependentField(eff_field);
    }
   }
  }  // end of outer else block

  std::string name = "GaAs_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_GaAs<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_GaAs<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold low-field mobility
  Kokkos::DynRankView<ScalarT,PHX::Device> lowFieldMob = Kokkos::createDynRankView(mobility.get_static_view(),"lowFieldMob",workset.num_cells, num_points);

  // Compute the low-field mobility at IP or BASIS depending on isEdgedl
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Obtain temperature [K]
      ScalarT latt = latt_temp(cell,point)*T0;

      // obtain doping and carrier concentrations
      ScalarT Na = acceptor(cell,point) * C0;
      ScalarT Nd = donor(cell,point) * C0;

      // set numerical limits
      if (Na <= 1.0) Na = 1.0;
      if (Nd <= 1.0) Nd = 1.0;

      // add the possible doping dependence
      ScalarT lfMob = evaluateLowFieldMobility(Na, Nd);

      // add the temperature dependence
      lowFieldMob(cell,point) = lfMob * std::pow(300.0/latt, tempExponent);  // in [cm2/V.s]
    }
  }

 // ----------------------------------------------------------------------------
 // Want mobility at the primary edge center for CVFEM-SG and EFFPG-FEM
 // ----------------------------------------------------------------------------

 if (isEdgedl)
 {
   // compute high-field mobility at the center of primary edge
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
   {
     for (int edge = 0; edge < num_edges; ++edge)
     {
       // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
       int node0 = cellType->getNodeMap(1,edge,0);
       int node1 = cellType->getNodeMap(1,edge,1);

       // primary edge lattice temperature in [K]
       ScalarT edgeLatt = (latt_temp(cell,node0) + latt_temp(cell,node1))*T0 / 2.0;

       // primary edge low-field mobility in [cm2/V.s]
       ScalarT edgelfMob = (lowFieldMob(cell,node0) + lowFieldMob(cell,node1)) / 2.0;
       ScalarT edgehfMob = edgelfMob;

       // enable the high field dependence when hiFieldOn = true
       if (hiFieldOn)
       {
         // get local nodes coordinates of a primary edge (each edge has 2 points)
         Kokkos::DynRankView<double,PHX::Device> edgePoints("edgePoints",2,num_dims);
         for (int dim = 0; dim < num_dims; ++dim)
         {
           edgePoints(0,dim) = (workset.bases[basis_index])->basis_coordinates(cell,node0,dim);
           edgePoints(1,dim) = (workset.bases[basis_index])->basis_coordinates(cell,node1,dim);
         }

         // compute the Philips-Thomas mobility for a primary edge
         edgehfMob = evaluateMobilityForEdgedl(cell, edge, edgelfMob, edgePoints, edgeLatt);

       }

       // primary edge mobility (scaled)
       mobility(cell,edge) = edgehfMob / Mu0;
     }
   }

 }  // end of if (isEdgedl)


 // ----------------------------------------------------------------------------
 // ----- Want mobility available at IP for SUPG-FEM ---------------------------
 // ----------------------------------------------------------------------------

 else
 {
   // loop over cells
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
   {
     for (int point = 0; point < num_points; ++point)
     {
       // Obtain temperature [K]
       ScalarT latt = latt_temp(cell,point)*T0;
       ScalarT lfMob = lowFieldMob(cell,point);
       ScalarT hfMob = lfMob;

       // enable the high field dependence when hiFieldOn = true
       if (hiFieldOn)  hfMob = evaluateMobilityForIPdl(cell, point, lfMob, latt);

       // mobility at IP (scaled)
       mobility(cell,point) = hfMob / Mu0;
     }
   }

 }  // end of else block

}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_GaAs<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up default parameters for the GaAs mobility model
  if (carrType == "Electron")
  {
    sign = 1.0;   // + for n and - for p
    lowMob = matProperty.getPropertyValue(matName, "Electron Mobility");  // default
    vsat300 = matProperty.getPropertyValue(matName, "Electron Saturation Velocity");
    Fsat = matProperty.getPropertyValue(matName, "Electron Saturation Field");
  }
  else if (carrType == "Hole")
  {
    sign = -1.0;   // + for n and - for p
    lowMob = matProperty.getPropertyValue(matName, "Hole Mobility");  // default
    vsat300 = matProperty.getPropertyValue(matName, "Hole Saturation Velocity");
    Fsat = matProperty.getPropertyValue(matName, "Hole Saturation Field");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite the default parameters when given in the input xml
  if (mobParamList.isParameter("Saturation Velocity"))
    vsat300 = mobParamList.get<double>("Saturation Velocity");

  if (mobParamList.isParameter("Saturation Field"))
    Fsat = mobParamList.get<double>("Saturation Field");

  isMobFromFile = false;

  if (mobParamList.isParameter("Low Field Mobility Value") && mobParamList.isParameter("Low Field Mobility File"))
  {
    std::stringstream msg;
    msg << "Error: Both Low Field Mobility Value and Low Field Mobility File are specified !"
        << "Please remove one of them !" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  else if (mobParamList.isParameter("Low Field Mobility Value"))
    lowMob = mobParamList.get<double>("Low Field Mobility Value");

  else if (mobParamList.isParameter("Low Field Mobility File"))
  {
    readLowFieldMobilityFile(mobParamList);
    isMobFromFile = true;
  }

  // read in the temperature exponent
  if (mobParamList.isParameter("Temperature Exponent"))
    tempExponent = mobParamList.get<double>("Temperature Exponent");
  else
    tempExponent = 0.0;

  // read in the temperature coefficient for the saturation velocity
  if (mobParamList.isParameter("Vsat Temperature Coefficient"))
    vsatTempCoeff = mobParamList.get<double>("Vsat Temperature Coefficient");
  else
    vsatTempCoeff = 0.0;

   // std::cout << "carrType=" << carrType <<", lowMob=" << lowMob << ", vsat=" << vsat << ", Fsat=" << Fsat << std::endl;

  // Is high field dependence on ?
  hiFieldOn = false;   // high field off by default
  if (mobParamList.isParameter("High Field"))
  {
    std::string highField = mobParamList.get<std::string>("High Field");
    if (highField == "Off")
      hiFieldOn = false;
    else if (highField == "On")
      hiFieldOn = true;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid High Field value ! Must be either Off or On !");
  }

  // Driving Force = ElectricField or GradQuasiFermi. It lives on IP for SUPG-FEM,
  // but lives on primary edge center for CVFEM-SG and EFFPG-FEM.
  if (mobParamList.isParameter("Driving Force"))
    driveForce = mobParamList.get<std::string>("Driving Force");
  else
    driveForce = "ElectricField";  // includes BGN contribution

  if ((driveForce != "ElectricField") && (driveForce != "GradQuasiFermi"))
     TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid Driving Force ! Must be either ElectricField or GradQuasiFermi !");
}


///////////////////////////////////////////////////////////////////////////////
//
//  readLowFieldMobilityFile()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_GaAs<EvalT, Traits>::readLowFieldMobilityFile
(const Teuchos::ParameterList& mobParamList)
{
  using std::string;
  using std::ifstream;

  string fileName = mobParamList.get<string>("Low Field Mobility File");
  ifstream inf(fileName.c_str());

  if (!inf)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error ! Cannot read the mobility file '" << fileName << "'" << std::endl);

  dopMobStruct dms;
  std::vector<dopMobStruct> dopMobVec;

  while (inf)
  {
    inf >> dms.dop >> dms.mob;
    dopMobVec.push_back(dms);
  }

  // sort the doping-mobility pairs according to doping values in ascending order
  std::sort(dopMobVec.begin(), dopMobVec.end());

  // eliminate duplicate (dop,mob) entries
  dopMobVec.resize( unique (dopMobVec.begin(), dopMobVec.end()) - dopMobVec.begin() );

  // assign the doping-mobility pairs to a map for easy search
  for (std::size_t i = 0; i< dopMobVec.size(); i++)
  {
    double dop = dopMobVec[i].dop;
    dopMobMap[dop] = dopMobVec[i].mob;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateLowFieldMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_GaAs<EvalT,Traits>::ScalarT
Mobility_GaAs<EvalT, Traits>::evaluateLowFieldMobility(const ScalarT& Na, const ScalarT& Nd)
{
  using std::map;

  ScalarT lfMob = 0.0;

  if (!isMobFromFile)  // use a constant low field mobility
    lfMob = lowMob;

  else  // use a doping-dependent mobility table read from a file
  {
    ScalarT Ntot = Na + Nd;

    typename map<ScalarT,ScalarT>::const_iterator d_beg = dopMobMap.begin();
    typename map<ScalarT,ScalarT>::const_iterator d_end = dopMobMap.end();

    // lower_bound() actually returns the point in the data that is
    // equal to or greater than the input Ntot. In that sense
    // it's really the "upper bound" as far as the search is concerned.
    typename map<ScalarT,ScalarT>::const_iterator d_ub = dopMobMap.lower_bound(Ntot);

    if (d_ub == d_beg)   // Before or at the first point
      lfMob = d_beg->second;

    else if (d_ub == d_end)  // After or at the last point
    {
      --d_ub;
      lfMob = d_ub->second;
    }
    else  // Somewhere between two points, do a linear interpolation
    {
      typename map<ScalarT,ScalarT>::const_iterator d_lb = d_ub;
      --d_lb;
      ScalarT slope = (d_ub->second - d_lb->second) / (d_ub->first - d_lb->first);
      lfMob = slope * (Ntot - d_lb->first) + d_lb->second;
    }
  }

  return lfMob;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateMobilityForEdgedl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_GaAs<EvalT,Traits>::ScalarT
Mobility_GaAs<EvalT, Traits>::evaluateMobilityForEdgedl (const std::size_t& cell, const int& edge,
  const ScalarT& edgelfMob, const Kokkos::DynRankView<double,PHX::Device>& edgePoints, const ScalarT& edgeLatt)
{
  // compute primary cell edge length and tangent
  Kokkos::DynRankView<double,PHX::Device> distance("distance",num_dims);
  double edgeLen = 0.0;

  for (int dim = 0; dim < num_dims; ++dim)
  {
    // get distance vector between the two local nodes of a primary edge
    distance(dim) = edgePoints(1,dim) - edgePoints(0,dim);

    // primary cell edge length
    edgeLen += distance(dim) * distance(dim);
  }
  edgeLen = std::sqrt(edgeLen);

  // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
  int node0 = cellType->getNodeMap(1,edge,0);
  int node1 = cellType->getNodeMap(1,edge,1);

  // compute the effective potential at local nodes
  ScalarT dEg0 = bandgap(cell,node0) - eff_bandgap(cell,node0) ;  // [eV]
  ScalarT dEg1 = bandgap(cell,node1) - eff_bandgap(cell,node1) ;
  ScalarT iEf0 = intrin_fermi(cell,node0);  // [eV]
  ScalarT iEf1 = intrin_fermi(cell,node1);

  ScalarT effPot0 = (sign*0.5*dEg0 - iEf0); // 1.0 converts [eV] to [V]
  ScalarT effPot1 = (sign*0.5*dEg1 - iEf1);

  // primary edge effective electric field in [V/cm], Eij = -(Potj-Poti)/hij
  ScalarT edgeField = -(effPot1-effPot0) / (edgeLen*X0);

  // default high field value
  ScalarT hiField = 1e-20;  // 0.0 initialization could cause nan in the Jacobian matrix;

  // driveForce == ElectricField, equivalent to the SGFVM high field calc. in charon1,
  // because primary edge current and edge field are parallel to each other
  if (driveForce == "ElectricField")
    hiField = std::abs(edgeField);

  // driveForce == GradQuasiFermi
  else if (driveForce == "GradQuasiFermi")
  {
    // get carrier density at local nodes
    ScalarT dens0 = 0.0;
    ScalarT dens1 = 0.0;
    if (carrType == "Electron")
    {
      dens0 = edensity(cell, node0);  // scaled
      dens1 = edensity(cell, node1);
    }
    else if (carrType == "Hole")
    {
      dens0 = hdensity(cell, node0);
      dens1 = hdensity(cell, node1);
    }

    // compute grad(dens)/dens in [1/cm] at the primary edge center, density scaling cancels out
    ScalarT gdens_over_dens = (dens1-dens0)/(edgeLen*X0) * 2./(dens1+dens0);

    // lattice temperature in [K]
    ScalarT lattT = edgeLatt;

    // obtain kb*T/q in [V]
    charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
    double kbBoltz = cpc.kb;          // Boltzmann constant in [eV/K]
    ScalarT kbT = kbBoltz*lattT/1.0;  // [V], 1.0 converts [eV] to [V]

    // gradient of quasi fermi potential at primary edge center in [V/cm]
    ScalarT gqfp = -sign*kbT*gdens_over_dens - edgeField;
    hiField = std::abs(gqfp);
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
         << "Invalid Driving Force ! Must be either ElectricField or GradQuasiFermi !");

  // initialize edgehfMob
  ScalarT edgehfMob = edgelfMob;

  // compute the temperature-dependent saturation velocity in [cm/s]
  ScalarT vsat = vsat300 / (1. - vsatTempCoeff + vsatTempCoeff * (edgeLatt/300.));

  if (carrType == "Electron")
  {
    ScalarT fieldRatio = std::pow(hiField, 3.) / std::pow(Fsat, 4.);
    edgehfMob = (edgelfMob + vsat * fieldRatio) / (1. + hiField * fieldRatio);
  }
  else if (carrType == "Hole")
    edgehfMob = edgelfMob / (1. + edgelfMob * hiField / vsat);

  return edgehfMob;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateMobilityForIPdl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_GaAs<EvalT,Traits>::ScalarT
Mobility_GaAs<EvalT, Traits>::evaluateMobilityForIPdl
(const std::size_t& cell, const int& point, const ScalarT& lfMob, const ScalarT& latt)
{
  // default high field value
  ScalarT hiField = 1e-20;  // 0.0 initialization could cause nan in the Jacobian matrix;

  // driveForce == "ElectricField", include BGN contribution (see charon::FEM_ElectricField)
  if (driveForce == "ElectricField")
  {
    for (int dim = 0; dim < num_dims; dim++)
    {
      const ScalarT& efield = eff_field(cell, point, dim);
      hiField += efield * efield;
    }
    hiField = std::sqrt(hiField) * E0;  // in [V/cm]
  }

  // driveForce == "GradQuasiFermi", include drift and diffusion contributions (see charon::FEM_ElectricField)
  else if (driveForce == "GradQuasiFermi")
  {
    for (int dim = 0; dim < num_dims; dim++)
    {
      const ScalarT& gqfp = grad_qfp(cell, point, dim);
      hiField += gqfp * gqfp;
    }
    hiField = std::sqrt(hiField) * E0;  // in [V/cm]
  }

  // initialize hfMob
  ScalarT hfMob = lfMob;

  // compute the temperature-dependent saturation velocity in [cm/s]
  ScalarT vsat = vsat300 / (1. - vsatTempCoeff + vsatTempCoeff * (latt/300.));

  if (carrType == "Electron")
  {
    ScalarT fieldRatio = std::pow(hiField, 3.) / std::pow(Fsat, 4.);
    hfMob = (lfMob + vsat * fieldRatio) / (1. + hiField * fieldRatio);
  }
  else if (carrType == "Hole")
    hfMob = lfMob / (1. + lfMob * hiField / vsat);

  return hfMob;

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_GaAs<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "GaAs", "GaAs mobility model");

  // for low-field part
  p->sublist("Mobility ParameterList").set<double>("Low Field Mobility Value", 0., "Specify a constant low field mobility in [cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<std::string>("Low Field Mobility File", "?", "Read the doping-dependent low field mobility from a file");
  p->sublist("Mobility ParameterList").set<double>("Temperature Exponent", 0., "Specify the temperature exponent for the low field mobility");

  // for high-field part
  p->sublist("Mobility ParameterList").set<std::string>("Driving Force", "ElectricField", "Different high field calculation methods");
  p->sublist("Mobility ParameterList").set<std::string>("High Field", "On", "Turn on/off high field dependence");
  p->sublist("Mobility ParameterList").set<double>("Saturation Velocity", 0., "[cm/s]");
  p->sublist("Mobility ParameterList").set<double>("Saturation Field", 0., "[V/cm]");
  p->sublist("Mobility ParameterList").set<double>("Vsat Temperature Coefficient", 0., "[1]");

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
