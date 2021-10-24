
#ifndef CHARON_MOBILITY_IONTEMPDEP_IMPL_HPP
#define CHARON_MOBILITY_IONTEMPDEP_IMPL_HPP

#include <cmath>
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"


/*
The model comes from the supporting document of the paper by Sungho Kim,
ShinHyun Choi, and Wei Lu, ACSNano, Vol.8, No.3, 2369-2376, 2014.

The electron mobility is computed from the electrical conductivity, i.e.,
mun = sigma / (q*n), where mun = electron mobility in [cm^2/(V.s)],
sigma = electrical conductivity in [ohms^(-1).cm^(-1)],
q = electron charge in [C], and n = electron density in [cm^(-3)].

The electrical conductivity sigma is modeled as sigma = sigma0 * exp(-Eac/kbT), where
sigma0 = minSigma0 [ohms^(-1).cm^(-1)] for ionDens <= minIonDens [cm^(-3)],
sigma0 = maxSigma0 for ionDens >= maxIonDens,
sigma0 = (maxSigma0 - minSigma0)/(maxIonDens - minIonDens)*(ionDens - minIonDens) + minSigma0,
for minIonDens < ionDens < maxIonDens, and
Eac = maxActE [eV] for ionDens <= minIonDens,
Eac = minActE for ionDens >= medIonDens,
Eac = (maxActE - minActE)/(minIonDens - medIonDens)*(ionDens - minIonDens) + maxActE.

Because of the exp(-Eac/kbT) dependence, to prevent sigma from becoming unphysically large,
sigma is upper-bounded by Maximum Electrical Conductivity maxSigma.

Specification of the electron mobility model in the input file takes the form of
<ParameterList name="Electron Mobility">
    <Parameter name="Value" type="string" value="IonTempDep" />
    <Parameter name="Maximum Ion Density" type="double" value="1e21" />
    <Parameter name="Minimum Ion Density" type="double" value="0" />
    <Parameter name="Medium Ion Density" type="double" value="5e20" />
    <Parameter name="Maximum Sigma0" type="double" value="940.0" />
    <Parameter name="Minimum Sigma0" type="double" value="10" />
    <Parameter name="Maximum Activation Energy" type="double" value="0.05" />
    <Parameter name="Minimum Activation Energy" type="double" value="-0.006" />
    <Parameter name="Maximum Electrical Conductivity" type="double" value="9400" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_IonTempDep<EvalT, Traits>::
Mobility_IonTempDep(
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

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize the mobility parameters
  initMobilityParams(mobParamList);

  // Evaluated field
  mobility = MDField<ScalarT,Cell,Point>(n.field.elec_mobility, scalar);
  this->addEvaluatedField(mobility);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  ion_density = MDField<const ScalarT,Cell,Point>(n.dof.iondensity,scalar);
  elec_density = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(ion_density);
  this->addDependentField(elec_density);

  std::string name = "IonTempDep_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_IonTempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb and q
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;       // Boltzmann constant in [eV/K]
  double eleQ = cpc.q;      // Electron charge in [C]

  // Compute the mobility at IP or BASIS
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // FieldContainer to temporarily hold the mobility values in a cell
    Kokkos::DynRankView<ScalarT,PHX::Device> tmpMobValue = Kokkos::createDynRankView(mobility.get_static_view(),"tmpMobValue",num_points);
    int numNonzeros = 0;

    for (int point = 0; point < num_points; ++point)
    {
      // Obtain ion density and electron density in [cm^-3]
      ScalarT ionDens = ion_density(cell,point) * C0;
      ScalarT eDens = elec_density(cell,point) * C0;

      // Compute sigma0 which linearly depends on the ion density with a positive slope
      ScalarT sigma0;  // [ohms^(-1).cm^(-1)]
      if (Sacado::ScalarValue<ScalarT>::eval(ionDens) <= minIonDens)
        sigma0 = minSigma0;
      else if (Sacado::ScalarValue<ScalarT>::eval(ionDens) >= maxIonDens)
        sigma0 = maxSigma0;
      else
        sigma0 = slopeSigma0 * (ionDens - minIonDens) + minSigma0;

      // Compute the activation energy which linearly depends on the ion density with a negative slope
      ScalarT actE;  // [eV]
      if (Sacado::ScalarValue<ScalarT>::eval(ionDens) <= minIonDens)
        actE = maxActE;
      else if (Sacado::ScalarValue<ScalarT>::eval(ionDens) >= medIonDens)
        actE = minActE;
      else
        actE = slopeActE * (ionDens - minIonDens) + maxActE;

      // Obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,point) * T0;
      ScalarT sigma;  // [ohms^(-1).cm^(-1)]

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so set sigma = sigma0
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)
        sigma = sigma0;
      else
      {
        ScalarT kbT = kb*lattT;   // [eV]
        sigma = sigma0 * std::exp(-actE / kbT);
      }

      // Reset sigma to maxSigma when sigma becomes unphysically large
      if (Sacado::ScalarValue<ScalarT>::eval(sigma) > maxSigma)  sigma = maxSigma;

      // Note eDens could become negative or a small positive number during simulation,
      // which leads to unphysical value for mobility
      if (Sacado::ScalarValue<ScalarT>::eval(eDens) <= 0.0)
      {
        mobility(cell,point) = 0.0;  // temporarily set to 0.0
        tmpMobValue(point) = 0.0;
      }

      // Compute electron mobility from sigma (electrical conductivity)
      else
      {
        ScalarT mobValue = sigma / (eleQ * eDens);  // [cm^2/(V.s)]
        if (Sacado::ScalarValue<ScalarT>::eval(mobValue) > maxMobValue)
          mobValue = maxMobValue;  // bound to maxMobValue
        mobility(cell,point) = mobValue / Mu0;

        // Save to a tmp vector for averaging
        numNonzeros += 1;
        tmpMobValue(point) = mobValue / Mu0;
      }
    }  // end of the loop over points

    // set entries with 0 value to the average of non-zero positive values in the cell
    for (int point = 0; point < num_points; ++point)
    {
      if (mobility(cell,point) <= 0.0)
      {
        ScalarT average = 0.0;
        for (int pt = 0; pt < num_points; ++pt) average += tmpMobValue(pt);
        average /= numNonzeros;
        mobility(cell,point) = average;
      }
    }

  }  // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_IonTempDep<EvalT, Traits>::initMobilityParams
(const Teuchos::ParameterList& mobParamList)
{
  // Set up parameters for the IonTempDep mobility model
  maxIonDens = mobParamList.get<double>("Maximum Ion Density");
  minIonDens = mobParamList.get<double>("Minimum Ion Density");
  medIonDens = mobParamList.get<double>("Medium Ion Density");
  maxSigma0 = mobParamList.get<double>("Maximum Sigma0");
  minSigma0 = mobParamList.get<double>("Minimum Sigma0");
  maxActE = mobParamList.get<double>("Maximum Activation Energy");
  minActE = mobParamList.get<double>("Minimum Activation Energy");
  maxSigma = mobParamList.get<double>("Maximum Electrical Conductivity");
  maxMobValue = mobParamList.get<double>("Maximum Electron Mobility");

  if (maxIonDens <= minIonDens)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Ion Density must be greater than Minimum Ion Density !");

  if (medIonDens <= minIonDens)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Medium Ion Density must be greater than Minimum Ion Density !");

  if (maxSigma0 <= minSigma0)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Sigma0 must be greater than Minimum Sigma0 !");

  if (maxActE <= minActE)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Activation Energy must be greater than Minimum Activation Energy !");

  if (maxSigma <= maxSigma0)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Electrical Conductivity must be greater than Maximum Sigma0 !");

  // Compute the positive slope of Sigma0 vs. Ion Density (linear relation)
  slopeSigma0 = (maxSigma0 - minSigma0) / (maxIonDens - minIonDens);

  // Compute the negative slope of ActE vs. Ion Density (linear relation)
  slopeActE = (maxActE - minActE) / (minIonDens - medIonDens);

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_IonTempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "IonTempDep", "Ion density and temperature dependent electron mobility");
  p->sublist("Mobility ParameterList").set<double>("Maximum Ion Density", 0., "[cm^(-3)]");
  p->sublist("Mobility ParameterList").set<double>("Minimum Ion Density", 0., "[cm^(-3)]");
  p->sublist("Mobility ParameterList").set<double>("Medium Ion Density", 0., "[cm^(-3)]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Sigma0", 0., "[ohms^(-1).cm^(-1)]");
  p->sublist("Mobility ParameterList").set<double>("Minimum Sigma0", 0., "[ohms^(-1).cm^(-1)]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Activation Energy", 0., "[eV]");
  p->sublist("Mobility ParameterList").set<double>("Minimum Activation Energy", 0., "[eV]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Electrical Conductivity", 0., "[ohms^(-1).cm^(-1)]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Electron Mobility", 0., "[cm^2/(V.s)]");


  return p;
}

}

#endif
