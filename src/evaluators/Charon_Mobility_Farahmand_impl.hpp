
#ifndef CHARON_MOBILITY_FARAHMAND_IMPL_HPP
#define CHARON_MOBILITY_FARAHMAND_IMPL_HPP

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


/*
The low-field electron mobility for GaN comes from
Farahmand, M. IEEE Trans. Elec. Devices vol48, No3. (Mar 2001) 535-542

For the SUPG-FEM formulation, hiField lives at IPs; while for CVFEM-SG and EFFPG-FEM,
hiField lives at centers of primary edges.

An example of using the model is given below:

<ParameterList name="Electron Mobility">
    <Parameter name="Value" type="string" value="Farahmand" />
</ParameterList>

<ParameterList name="Hole Mobility">
    <Parameter name="Value" type="string" value="Farahmand" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Farahmand<EvalT, Traits>::
Mobility_Farahmand(
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

  // Initialize the Farahmand mobility parameters, must be called here since
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
  //acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor,input_scalar);
  //donor = MDField<const ScalarT,Cell,Point>(n.field.donor,input_scalar);
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

  std::string name = "Farahmand_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Farahmand<EvalT, Traits>::
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
Mobility_Farahmand<EvalT, Traits>::
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
      ScalarT Na = acceptor(cell,point);
      ScalarT Nd = donor(cell,point);
      //ScalarT Ntot = (Na + Nd) * C0; // unscaled conc.


      ScalarT lfMob = evaluateLowFieldMobility(Na, Nd, latt);
      lowFieldMob(cell,point) = lfMob ;  // in [cm2/V.s]
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
       ScalarT lfMob = lowFieldMob(cell,point);
       ScalarT hfMob = lfMob;

       // enable the high field dependence when hiFieldOn = true
       if (hiFieldOn)  hfMob = evaluateMobilityForIPdl(cell, point, lfMob);

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
void Mobility_Farahmand<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up default parameters for the Farahmand mobility model
  if (carrType == "Electron")
  {
    sign = 1.0;   // + for n and - for p
    mu1 = matProperty.getPropertyValue(matName, "Farahmand MU1N");
    mu2 = matProperty.getPropertyValue(matName, "Farahmand MU2N");
    alpha = matProperty.getPropertyValue(matName, "Farahmand ALPHAN");
    beta = matProperty.getPropertyValue(matName, "Farahmand BETAN");
    delta = matProperty.getPropertyValue(matName, "Farahmand DELTAN");
    gamma = matProperty.getPropertyValue(matName, "Farahmand GAMMAN");
    eps = matProperty.getPropertyValue(matName, "Farahmand EPSN");
    ncrit = matProperty.getPropertyValue(matName, "Farahmand NCRITN");
  }
  else if (carrType == "Hole")
  {
    sign = -1.0;   // + for n and - for p
    mu1 = matProperty.getPropertyValue(matName, "Farahmand MU1P");
    mu2 = matProperty.getPropertyValue(matName, "Farahmand MU2P");
    alpha = matProperty.getPropertyValue(matName, "Farahmand ALPHAP");
    beta = matProperty.getPropertyValue(matName, "Farahmand BETAP");
    delta = matProperty.getPropertyValue(matName, "Farahmand DELTAP");
    gamma = matProperty.getPropertyValue(matName, "Farahmand GAMMAP");
    eps = matProperty.getPropertyValue(matName, "Farahmand EPSP");
    ncrit = matProperty.getPropertyValue(matName, "Farahmand NCRITP");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Overwrite the default parameters when given in the input xml
  if (mobParamList.isParameter("mu1"))
    mu1 = mobParamList.get<double>("mu1");
  if (mobParamList.isParameter("mu2"))
    mu2 = mobParamList.get<double>("mu2");
  if (mobParamList.isParameter("alpha"))
    alpha = mobParamList.get<double>("alpha");
  if (mobParamList.isParameter("beta"))
    beta = mobParamList.get<double>("beta");
  if (mobParamList.isParameter("delta"))
    delta = mobParamList.get<double>("delta");
  if (mobParamList.isParameter("gamma"))
    gamma = mobParamList.get<double>("gamma");
  if (mobParamList.isParameter("eps"))
    eps = mobParamList.get<double>("eps");
  if (mobParamList.isParameter("ncrit"))
    ncrit = mobParamList.get<double>("ncrit");

  // Is high field dependence on ?
  hiFieldOn = false;   // high field off by default
  if (mobParamList.isParameter("High Field"))
  {
    std::string highField = mobParamList.get<std::string>("High Field");
    if (highField == "Off")
      hiFieldOn = false;
    else if (highField == "On")
    {
      hiFieldOn = true;
      if (carrType == "Electron")
      {
        n1 = matProperty.getPropertyValue(matName, "N1N GANSAT");
        n2 = matProperty.getPropertyValue(matName, "N2N GANSAT");
        an = matProperty.getPropertyValue(matName, "ANN GANSAT");
        ec = matProperty.getPropertyValue(matName, "ECN GANSAT");
        vsat = matProperty.getPropertyValue(matName, "VSATN");
      }
      else if (carrType == "Hole")
      {
        n1 = matProperty.getPropertyValue(matName, "N1P GANSAT");
        n2 = matProperty.getPropertyValue(matName, "N2P GANSAT");
        an = matProperty.getPropertyValue(matName, "ANP GANSAT");
        ec = matProperty.getPropertyValue(matName, "ECP GANSAT");
        vsat = matProperty.getPropertyValue(matName, "VSATP");
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
          << "Invalid Carrier Type ! Must be either Electron or Hole !");

      if (mobParamList.isParameter("n1"))
        n1 = mobParamList.get<double>("n1");
      if (mobParamList.isParameter("n2"))
        n2 = mobParamList.get<double>("n2");
      if (mobParamList.isParameter("an"))
        an = mobParamList.get<double>("an");
      if (mobParamList.isParameter("ec"))
        ec = mobParamList.get<double>("ec");
      if (mobParamList.isParameter("vsat"))
        vsat = mobParamList.get<double>("vsat");

    }
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
//  evaluateLowFieldMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Farahmand<EvalT,Traits>::ScalarT
Mobility_Farahmand<EvalT, Traits>::evaluateLowFieldMobility(const ScalarT& Na, const ScalarT& Nd, const ScalarT& latt)
{

  ScalarT lfMob = 1e-20;
  ScalarT Ntot = (Na + Nd) * C0;
  ScalarT latratio = latt / 300.0;
  ScalarT nratio = Ntot / (ncrit * std::pow(latratio,gamma));
  ScalarT powbit = alpha * std::pow(latratio, eps);

  lfMob = mu1 * std::pow(latratio,beta) +
    ((mu2 - mu1) * std::pow(latratio,delta)) / (1 + std::pow(nratio,powbit));

  return lfMob;
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateMobilityForEdgedl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Farahmand<EvalT,Traits>::ScalarT
Mobility_Farahmand<EvalT, Traits>::evaluateMobilityForEdgedl (const std::size_t& cell, const int& edge,
  const ScalarT& edgelfMob, const Kokkos::DynRankView<double,PHX::Device>& edgePoints, const ScalarT& edgeLatt)
{
  ScalarT edgehfMob = 1e-20;
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
    ScalarT kbT = kbBoltz*lattT;

    // gradient of quasi fermi potential at primary edge center in [V/cm]
    ScalarT gqfp = -sign*kbT*gdens_over_dens - edgeField;
    hiField = std::abs(gqfp);
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
         << "Invalid Driving Force ! Must be either ElectricField or GradQuasiFermi !");

  // initialize edgehfMob

  edgehfMob = evaluateHighFieldMobility(edgelfMob, hiField);

  return edgehfMob;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateMobilityForIPdl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Farahmand<EvalT,Traits>::ScalarT
Mobility_Farahmand<EvalT, Traits>::evaluateMobilityForIPdl
(const std::size_t& cell, const int& point, const ScalarT& lfMob)
{
  // default high field value
  ScalarT hiField = 1e-20;  // 0.0 initialization could cause nan in the Jacobian matrix;
  ScalarT hfMob = 1e-20;

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

  hfMob = evaluateHighFieldMobility(lfMob, hiField);

  return hfMob;

}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateHighFieldMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Farahmand<EvalT,Traits>::ScalarT
Mobility_Farahmand<EvalT, Traits>::evaluateHighFieldMobility
(const ScalarT& lfMob, const ScalarT& hiField)
{
  ScalarT hfMob = lfMob;
  ScalarT eratio = 0.0;
  eratio = hiField / ec;
  hfMob = (lfMob + vsat * std::pow(hiField,(n1 - 1.0)) / std::pow(ec,n1)) /
    ((1.0 + an * std::pow(eratio,n2)) + (std::pow(eratio,n1)));

  return hfMob;
}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Farahmand<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Farahmand", "Farahmand mobility model");

  // for low-field part
  //p->sublist("Mobility ParameterList").set<std::string>("Low Field Mobility Model", "?", "Farahmand");
  p->sublist("Mobility ParameterList").set<double>("mu1", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mu2", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("alpha", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("beta", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("delta", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("gamma", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("eps", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("ncrit", 0., "[cm^-3]");
  // for high-field part
  p->sublist("Mobility ParameterList").set<std::string>("Driving Force", "ElectricField", "Different high field calculation methods");
  p->sublist("Mobility ParameterList").set<std::string>("High Field", "On", "Turn on/off high field dependence");
  p->sublist("Mobility ParameterList").set<double>("n1", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("n2", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("an", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("ec", 0., "[V/cm]");
  p->sublist("Mobility ParameterList").set<double>("vsat", 0., "[cm/s]");

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
