
#ifndef CHARON_MOBILITY_LUCENT_IMPL_HPP
#define CHARON_MOBILITY_LUCENT_IMPL_HPP

#include <cmath>
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
The Lucent mobility model is by M. D. Darwish et al., "An improved electron and
hole mobility model for general purpose device simulation," IEEE Transactions on
Electron Devices, vol.44, no.9, pp.1529-1538, 1997.

The implementation here only includes the low- and high-longitudinal field dependences,
while excluding the transverse-field dependence. The transverse-field dependence
requires an accurate calculation of the transverse field, which needs to be investigated
in detail, so leave it for future work.

The low-longitudinal field part is the same as the low-field Philips model (see
charon::Mobility_Philips) except for two things:

(1) Nsceff(e) = Nd + G(Pe)*Na, while it is Nd + G(Pe)*Na + p/F(Pe) in Philips
    Nsceff(h) = Na + G(Ph)*Nd, while it is Na + G(Ph)*Nd + n/F(Ph) in Philips
(2) it does not include clustering, i.e.,
    Nd = Nd0 and Na = Na0, while Nd = Nd0*Zd and Na = Na0*Za in Philips

The high-field dependence uses the expression of
hfMob = 2.*lfMob/(1.+ pow(1.+ 4.*pow(lfMob * hiField /vsat, 2.), 1./2.) );

where lfMob is the low-field mobility, and hiField is computed by one of
the two methods:

(1) when Driving Force = ElectricField (default), hiField = abs(effective electric field) that
includes BGN contribution, i.e., Fn/p,eff = abs(grad(Ei/V0 -/+ 0.5*dEg/V0)) (scaled).

(2) when Driving Force = GradQuasiFermi, hiField = abs(gradient of quasi-fermi
potential) that includes both drift and diffusion contributions, i.e.,
grad_qfp = abs(-grad(n)/n-Fn,eff) for n and = abs(grad(p)/p-Fp,eff) for p (scaled), valid
for isothermal simulation.

For the SUPG-FEM formulation, hiField lives at IPs; while for CVFEM-SG and EFFPG-FEM,
hiField lives at centers of primary edges.

When Driving Force = ElectricField, the calculation here is equivalent to the SGFVM
high field calculation in charon1 for the SG scheme. The reason is that, even though
Charon_SGFVM_CurrentDensityT.h in charon1 uses an approximate current density, because the
current density and electric field used are both parallel to primary edges, the result
is just the electric field parallel to primary edges.

*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Lucent<EvalT, Traits>::
Mobility_Lucent(
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

  // Initialize Lucent mobility parameters, must be called here since
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

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,input_scalar);
  //acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor,input_scalar);
  //donor = MDField<const ScalarT,Cell,Point>(n.field.donor,input_scalar);
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,input_scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,input_scalar);
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,input_scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,input_scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(acceptor);
  this->addDependentField(donor);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);

  if (isEdgedl)  // for CVFEM-SG and EFFPG-FEM
  {
    // Always need these fields when isEdgedl = true
    intrin_fermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi, input_scalar);
    bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap, input_scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, input_scalar);

    this->addDependentField(intrin_fermi);
    this->addDependentField(bandgap);
    this->addDependentField(eff_bandgap);
  }

  else  // for SUPG-FEM
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

  }  // end of outer else block

  std::string name = "Lucent_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Lucent<EvalT, Traits>::
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
Mobility_Lucent<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold low-field Philips mobility
  Kokkos::DynRankView<ScalarT,PHX::Device> philipsMob = Kokkos::createDynRankView(mobility.get_static_view(),"philipsMob",workset.num_cells, num_points);

  // Compute the low-field Philips mobility at IP or BASIS depending on isEdgedl
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Obtain temperature [K]
      ScalarT latt = latt_temp(cell,point)*T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0;

      // Compute vsat in [cm/s]
      if (vsat < 1e-10)  // user did not specify the value
         vsat = 2.4e7/(1.+0.8*exp(latt/600.));  // also in charon1 manual

      // Obtain doping and carrier concentrations
      ScalarT Na = acceptor(cell,point) * C0;
      ScalarT Nd = donor(cell,point) * C0;
      ScalarT eden = edensity(cell,point) * C0;
      ScalarT hden = hdensity(cell,point) * C0;

      // set numerical limits
      if (Na <= 1.0) Na = 1.0;
      if (Nd <= 1.0) Nd = 1.0;
      if (eden <= 1.0) eden = 1.0;
      if (hden <= 1.0) hden = 1.0;

      ScalarT bulkMob = evaluatePhilipsMobility(Na, Nd, eden, hden, latt);
      philipsMob(cell,point) = bulkMob ;  // in [cm2/V.s]
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

       // edgeLatt should be always > 0, but it could become <= 0 due to numerical errors
       // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
       if (Sacado::ScalarValue<ScalarT>::eval(edgeLatt) <= 0.0)  edgeLatt = 300.0;

       // primary edge low-field mobility in [cm2/V.s]
       ScalarT edgelfMob = (philipsMob(cell,node0) + philipsMob(cell,node1)) / 2.0;
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

         // compute the Lucent mobility for a primary edge
         edgehfMob = evalLucentMobForEdgedl(cell, edge, edgelfMob, edgePoints, edgeLatt);
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
       ScalarT lfMob = philipsMob(cell,point);
       ScalarT hfMob = lfMob;

       // enable the high field dependence when hiFieldOn = true
       if (hiFieldOn)
       {
         // compute the Lucent mobility at IP
         hfMob = evalLucentMobForIPdl(cell, point, lfMob);
       }

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
void Mobility_Lucent<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for Lucent electron mobility model
  if (carrType == "Electron")
  {
    string species = "Arsenic";  // by default
    if (mobParamList.isParameter("Dopant Species"))
      species = mobParamList.get<string>("Dopant Species");

    // Retrieve parameters from charon::Material_Properties by default
    if (species == "Arsenic")
    {
      mumax = matProperty.getPropertyValue(matName, "Philips Electron(As) mumax");
      mumin = matProperty.getPropertyValue(matName, "Philips Electron(As) mumin");
      gamma = matProperty.getPropertyValue(matName, "Philips Electron(As) gamma");
      nref = matProperty.getPropertyValue(matName, "Philips Electron(As) nref");
      alpha = matProperty.getPropertyValue(matName, "Philips Electron(As) alpha");
    }
    else if (species == "Phosphorous")
    {
      mumax = matProperty.getPropertyValue(matName, "Philips Electron(P) mumax");
      mumin = matProperty.getPropertyValue(matName, "Philips Electron(P) mumin");
      gamma = matProperty.getPropertyValue(matName, "Philips Electron(P) gamma");
      nref = matProperty.getPropertyValue(matName, "Philips Electron(P) nref");
      alpha = matProperty.getPropertyValue(matName, "Philips Electron(P) alpha");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid Dopant Species ! Must be either Arsenic or Phosphorous !");

    // Set the exponent parameter in the high field part
    beta_hf = 2.0;
    sign = 1.0;   // + for n and - for p
  }

  // Set up parameters for Lucent hole mobility model
  else if (carrType == "Hole")
  {
    // Retrieve parameters from charon::Material_Properties by default
    mumax = matProperty.getPropertyValue(matName, "Philips Hole mumax");
    mumin = matProperty.getPropertyValue(matName, "Philips Hole mumin");
    gamma = matProperty.getPropertyValue(matName, "Philips Hole gamma");
    nref = matProperty.getPropertyValue(matName, "Philips Hole nref");
    alpha = matProperty.getPropertyValue(matName, "Philips Hole alpha");

    // Set the exponent parameter in the high field part
    beta_hf = 2.0;
    sign = -1.0;   // + for n and - for p
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Retrieve carrier-independent parameters
  nref_d = matProperty.getPropertyValue(matName, "Philips Donor nref_d");
  nref_a = matProperty.getPropertyValue(matName, "Philips Acceptor nref_a");
  cref_d = matProperty.getPropertyValue(matName, "Philips Donor cref_d");
  cref_a = matProperty.getPropertyValue(matName, "Philips Acceptor cref_a");

  // Overwrite parameters when specified by users
  if (mobParamList.isParameter("mumax"))
    mumax = mobParamList.get<double>("mumax");
  if (mobParamList.isParameter("mumin"))
    mumin = mobParamList.get<double>("mumin");
  if (mobParamList.isParameter("gamma"))
    gamma = mobParamList.get<double>("gamma");
  if (mobParamList.isParameter("nref"))
    nref = mobParamList.get<double>("nref");
  if (mobParamList.isParameter("alpha"))
    alpha = mobParamList.get<double>("alpha");
  if (mobParamList.isParameter("nref_d"))
    nref_d = mobParamList.get<double>("nref_d");
  if (mobParamList.isParameter("nref_a"))
    nref_a = mobParamList.get<double>("nref_a");
  if (mobParamList.isParameter("cref_d"))
    cref_d = mobParamList.get<double>("cref_d");
  if (mobParamList.isParameter("cref_a"))
    cref_a = mobParamList.get<double>("cref_a");
  if (mobParamList.isParameter("beta_hf"))
    beta_hf = mobParamList.get<double>("beta_hf");

  // Set G Parameters 
  gParaSet = "Klaassen";  // by default
  if (mobParamList.isParameter("G Parameter Set"))
    gParaSet = mobParamList.get<string>("G Parameter Set");

  if (gParaSet == "Klaassen")
  {
    ag = 0.89233;
    bg = 0.41372;
    cg = 0.005978;
    alpha_g = 0.28227;
    alpha_prime_g = 0.72169;
    beta_g = 0.19778;
    gamma_g = 1.80618;
  }
  else if (gParaSet == "Meyer")
  {
    ag = 4.41804;
    bg = 39.9014;
    cg = 0.52896;
    alpha_g = 0.0001;
    alpha_prime_g = 1.595187;
    beta_g = 0.38297;
    gamma_g = 0.25948;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid G Parameter Set ! Must be either Klaassen or Meyer !");

  // Is high field dependence on ?
  hiFieldOn = true;   // high field on by default
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

  // Set saturation velocity when specified by user
  vsat = 0.; // initialization
  if (mobParamList.isParameter("Saturation Velocity"))
    vsat = mobParamList.get<double>("Saturation Velocity");

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
//  evaluatePhilipsMobility()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Lucent<EvalT,Traits>::ScalarT
Mobility_Lucent<EvalT, Traits>::evaluatePhilipsMobility(const ScalarT& Na, const ScalarT& Nd,
   const ScalarT& eden, const ScalarT& hden, const ScalarT& latt)
{
   int maxCount = 500;

   ScalarT refTemp = 300.0;
   ScalarT normTemp = latt / refTemp;

   ScalarT lattMob = mumax * pow(normTemp, gamma);
   ScalarT dopMob = mumax*mumax / (mumax-mumin) *pow(normTemp, (3.0*alpha-1.5) );
   ScalarT carrMob = mumax*mumin / (mumax-mumin) *pow(normTemp, -0.5);

   // Naeff and Ndeff describe clustering effects of dopants at ultrahigh conc.
   // ScalarT Naeff = Na*(1.0 + Na*Na/(cref_a*Na*Na +nref_a*nref_a) );
   // ScalarT Ndeff = Nd*(1.0 + Nd*Nd/(cref_d*Nd*Nd +nref_d*nref_d) );

   // do not include clustering for the Lucent model
   ScalarT Naeff = Na;
   ScalarT Ndeff = Nd;

   // carrier-dependent parameters
   ScalarT Nsc = 1.0;
   ScalarT meff = 1.0;         // in unit of [m0]
   //ScalarT mratio = 1.0;
   if (carrType == "Electron")
   {
      Nsc = Naeff + Ndeff + hden;
      meff = 1.0;
      //mratio = 1.0/1.258;
   }
   else if (carrType == "Hole")
   {
      Nsc = Naeff + Ndeff + eden;
      meff = 1.258;
      //mratio = 1.258/1.0;
   }

   // Psc = screening parameter
   ScalarT tmp1 = 2.459 / (3.97e13 * pow(Nsc,-2./3.) );
   ScalarT tmp2 = 3.828 * (eden + hden) / (1.36e20*meff);
   ScalarT Psc = 1.0/(tmp1 + tmp2) * pow(normTemp, 2.0);

   // Find Pmin for the "Klaassen" G parameters

   // Suzey's code based on Newton-Raphson iteration to solve dG/dP = 0 for Pmin
   ScalarT Pmin = 0.3246;   // initial guess from Charon 1.0
   ScalarT FPmin = 1.0;     // arbitrary values to start the loop
   ScalarT dFPmin = 1.0;
   // converged when |P(k)-P(k-1)| < 1e-5
   int count = 0;
   while ( std::abs(Sacado::ScalarValue<ScalarT>::eval(FPmin/dFPmin)) > 1e-5 &&
           count < maxCount)
   {
      FPmin = cg*gamma_g/(ag*beta_g) * pow(meff/normTemp, (-alpha_prime_g*gamma_g))
            * pow(Pmin, (-gamma_g-1.0)) - pow(normTemp/meff, alpha_g)
            * pow(bg+Pmin*pow(normTemp/meff,alpha_g), (-beta_g-1.0));
      dFPmin = cg*gamma_g*(-gamma_g-1.0)/(ag*beta_g) * pow(meff/normTemp, (-alpha_prime_g*gamma_g))
             * pow(Pmin, (-gamma_g-2.0)) + (beta_g+1.0) * pow(normTemp/meff, 2.0*alpha_g)
             * pow(bg+Pmin*pow(normTemp/meff,alpha_g), (-beta_g-2.0));
      Pmin = Pmin - FPmin/dFPmin;
      count++;
   }
   if (count >= maxCount)
     std::cerr << "WARNING: iteration exceeded maximum allowed in " << __FILE__ << std::endl;

   // Gp describes minority carrier-impurity scattering
   ScalarT Gp = 1.0;
   tmp1 = pow(bg + Psc*pow(normTemp/meff, alpha_g), -beta_g);
   tmp2 = pow(Psc*pow(meff/normTemp, alpha_prime_g), -gamma_g);
   Gp = 1.0 - ag*tmp1 + cg*tmp2;

   if ( (gParaSet == "Klaassen") && (Psc < Pmin) )
   {
      tmp1 = pow(bg + Pmin*pow(normTemp/meff, alpha_g), -beta_g);
      tmp2 = pow(Pmin*pow(meff/normTemp, alpha_prime_g), -gamma_g);
      Gp = 1.0 - ag*tmp1 + cg*tmp2;
   }

   // Fp describes electron-hole scattering effect
   // ScalarT Fp = 1.0;
   // Fp = (0.7643*pow(Psc,0.6478) + 2.2999 + 6.5502*mratio)
   //   / (pow(Psc,0.6478) + 2.3670 - 0.8552*mratio);

   // Nsceff describes the effective carrier-impurity scattering
   ScalarT Nsceff = 1.0;
   if (carrType == "Electron")
   {
      // Nsceff = Ndeff + Gp*Naeff + hden/Fp;
      Nsceff = Ndeff + Gp*Naeff;
   }
   else if (carrType == "Hole")
   {
      // Nsceff = Naeff + Gp*Ndeff + eden/Fp;
      Nsceff = Naeff + Gp*Ndeff;
   }

   // Mobility calculation
   ScalarT muDAeh = dopMob*(Nsc/Nsceff)*pow(nref/Nsc, alpha)
                  + carrMob*(eden + hden)/Nsceff;
   ScalarT bulkMob = 1.0 / (1.0/lattMob + 1.0/muDAeh);

   return bulkMob;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalLucentMobForEdgedl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Lucent<EvalT,Traits>::ScalarT
Mobility_Lucent<EvalT, Traits>::evalLucentMobForEdgedl (const std::size_t& cell, const int& edge,
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

  // check whether the value is zero within machine precision
  if (Sacado::ScalarValue<ScalarT>::eval(hiField) > std::numeric_limits<double>::epsilon())
    edgehfMob = 2.*edgelfMob/(1.+ pow(1.+ 4.*pow(edgelfMob * hiField /vsat, beta_hf), 1./beta_hf) );

  return edgehfMob;

}


///////////////////////////////////////////////////////////////////////////////
//
//  evalLucentMobForIPdl()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Mobility_Lucent<EvalT,Traits>::ScalarT
Mobility_Lucent<EvalT, Traits>::evalLucentMobForIPdl
(const std::size_t& cell, const int& point, const ScalarT& lfMob)
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

  // check whether the value is zero within machine precision
  if (Sacado::ScalarValue<ScalarT>::eval(hiField) > std::numeric_limits<double>::epsilon())
    hfMob = 2.*lfMob/(1.+ pow(1.+ 4.*pow(lfMob * hiField /vsat, beta_hf), 1./beta_hf) );

  return hfMob;

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Lucent<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Carrier Type", "?");
  p->set<std::string>("EquationSet Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Lucent", "Lucent mobility model");

  // for low-field Philips part
  p->sublist("Mobility ParameterList").set<std::string>("Dopant Species", "Arsenic", "Arsenic or Phosphorous as n-type dopant");
  p->sublist("Mobility ParameterList").set<std::string>("G Parameter Set", "Klaassen", "Use Klaassen or Meyer parameter set");
  p->sublist("Mobility ParameterList").set<double>("mumax", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("mumin", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("nref", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("gamma", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("alpha", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("nref_d", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("nref_a", 0., "[cm^-3]");
  p->sublist("Mobility ParameterList").set<double>("cref_d", 0., "[1]");
  p->sublist("Mobility ParameterList").set<double>("cref_a", 0., "[1]");

  // for high-field dependence
  p->sublist("Mobility ParameterList").set<std::string>("Driving Force", "ElectricField", "Different high field calculation methods");
  p->sublist("Mobility ParameterList").set<std::string>("High Field", "On", "Turn on/off high field dependence");
  p->sublist("Mobility ParameterList").set<double>("Saturation Velocity", 0., "[V/cm]");

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
