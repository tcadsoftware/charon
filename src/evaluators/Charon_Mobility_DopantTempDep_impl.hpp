
#ifndef CHARON_MOBILITY_DOPANTTEMPDEP_IMPL_HPP
#define CHARON_MOBILITY_DOPANTTEMPDEP_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"


/*
The model comes from the supporting document of the paper by Sungho Kim,
ShinHyun Choi, and Wei Lu, ACSNano, Vol.8, No.3, 2369-2376, 2014.

Since sigma = q*n*mun = sigma0*exp(-Eac/kbT) = c*n*exp(-Eac/kbT), with c = constant,
mun can be modeled as mun = mun0 * exp(-Eac/kbT), with mun0 = constant.

Eac = minActE [eV] for dopDens > maxDens,
Eac = (minActE-maxActE)/(maxDens-minDens) * (dopDens-minDens) + maxActE for minDens <= dopDens <= maxDens,
Eac = maxActE for dopDens < minDens.

Specification of the electron mobility model in the input file takes the form of
<ParameterList name="Electron Mobility">
    <Parameter name="Value" type="string" value="DopantTempDep" />
    <Parameter name="Mobility Multiplier" type="double" value="10" />
    <Parameter name="Maximum Dopant Density" type="double" value="4e21" />
    <Parameter name="Minimum Dopant Density" type="double" value="0" />
    <Parameter name="Maximum Activation Energy" type="double" value="0.1" />
    <Parameter name="Minimum Activation Energy" type="double" value="-0.006" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_DopantTempDep<EvalT, Traits>::
Mobility_DopantTempDep(
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
  dopant_density = MDField<const ScalarT,Cell,Point>(p.get<std::string>("Dopant Name"),scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(dopant_density);

  std::string name = "DopantTempDep_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_DopantTempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb and q
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;       // Boltzmann constant in [eV/K]

  // Compute the mobility at IP or BASIS
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Obtain dopant density in [cm^-3]
      ScalarT dopDens = dopant_density(cell,point) * C0;

      // Compute the activation energy which linearly depends on the dopant density with a negative slope
      ScalarT actE;  // [eV]
      if (Sacado::ScalarValue<ScalarT>::eval(dopDens) < minDopDens)
        actE = maxActE;
      else if (Sacado::ScalarValue<ScalarT>::eval(dopDens) > maxDopDens)
        actE = minActE;
      else
        actE = slopeActE * (dopDens - minDopDens) + maxActE;

      // Obtain lattice temperature in [K]
      ScalarT lattT = latt_temp(cell,point) * T0;

      ScalarT mobValue;  // [cm^2/(V.s)]

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so set mu = mu0
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)
        mobValue = mobVal0;
      else
      {
        ScalarT kbT = kb*lattT;   // [eV]
        mobValue = mobVal0 * std::exp(-actE / kbT);
      }

      mobility(cell,point) = mobValue / Mu0;

    }  // end of the loop over points
  }  // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_DopantTempDep<EvalT, Traits>::initMobilityParams
(const Teuchos::ParameterList& mobParamList)
{
  // Set up parameters for the DopantTempDep mobility model
  mobVal0 = mobParamList.get<double>("Mobility Multiplier");
  maxDopDens = mobParamList.get<double>("Maximum Dopant Density");
  minDopDens = mobParamList.get<double>("Minimum Dopant Density");
  maxActE = mobParamList.get<double>("Maximum Activation Energy");
  minActE = mobParamList.get<double>("Minimum Activation Energy");

  if (maxDopDens < minDopDens)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Dopant Density must be not smaller than Minimum Dopant Density !");

  if (maxActE < minActE)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error: Maximum Activation Energy must be not smaller than Minimum Activation Energy !");

  // Compute the negative slope of ActE vs. dopant density (linear relation)
  slopeActE = (maxActE - minActE) / (minDopDens - maxDopDens);

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_DopantTempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "DopantTempDep", "Ion density and temperature dependent electron mobility");
  p->sublist("Mobility ParameterList").set<double>("Mobility Multiplier", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Dopant Density", 0., "[cm^(-3)]");
  p->sublist("Mobility ParameterList").set<double>("Minimum Dopant Density", 0., "[cm^(-3)]");
  p->sublist("Mobility ParameterList").set<double>("Maximum Activation Energy", 0., "[eV]");
  p->sublist("Mobility ParameterList").set<double>("Minimum Activation Energy", 0., "[eV]");

  p->set<std::string>("Dopant Name","?", "Dopant Name = Donor Concentration or ION_DENSITY");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
