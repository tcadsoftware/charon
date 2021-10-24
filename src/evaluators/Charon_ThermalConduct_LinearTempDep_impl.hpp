
#ifndef CHARON_THERMALCONDUCT_LINEARTEMPDEP_IMPL_HPP
#define CHARON_THERMALCONDUCT_LINEARTEMPDEP_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"


/*
The model comes from the supporting document of the paper by Sungho Kim,
ShinHyun Choi, and Wei Lu, ACSNano, Vol.8, No.3, 2369-2376, 2014.

The thermal conductivity is modeled as kappa = kappa0*[1+lambda*(T-Tref)], where
kappa0 = Thermal Conductivity at Reference Temperature in [W/(cm.K)],
lambda = Linear Thermal Coefficient in [1/K],
Tref = Reference Temperature in [K].
kappa0, lambda, and Tref must be given by user in the input xml file.

Specification of the thermal conductivity model in the input file takes the form of
<ParameterList name="Thermal Conductivity">
  <Parameter name="Value" type="string" value="LinearTempDep" />
  <Parameter name="Thermal Conductivity at Reference Temperature" type="double" value="0.0012" />
  <Parameter name="Linear Thermal Coefficient" type="double" value="0.1" />
  <Parameter name="Reference Temperature" type="double" value="300" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
ThermalConduct_LinearTempDep<EvalT, Traits>::
ThermalConduct_LinearTempDep(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Thermal Conductivity ParameterList
  const ParameterList& plist = p.sublist("Thermal Conductivity ParameterList");

  // Initialize the thermal conductivity model parameters
  initialize(plist);

  // Evaluated fields
  ther_cond = MDField<ScalarT,Cell,Point>(n.field.kappa,scalar);
  this->addEvaluatedField(ther_cond);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  kL0 = scaleParams->scale_params.kL0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,scalar);

  this->addDependentField(latt_temp);

  std::string name = "ThermalConductivity_LinearTempDep";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
ThermalConduct_LinearTempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      const ScalarT& scaledT = latt_temp(cell,point);  // scaled
      ScalarT TL = T0*scaledT;  // in [K]
      ScalarT tc;

      // lattice temperature should be always >= Tref
      if (Sacado::ScalarValue<ScalarT>::eval(TL) >= Tref)
        tc = kappa0 * (1.0 + lambda * (TL - Tref));    // in [W/(K.cm)]
      else
        tc = kappa0;

      ther_cond(cell,point) = tc / kL0;  // scaled
    }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void ThermalConduct_LinearTempDep<EvalT, Traits>::initialize(const Teuchos::ParameterList& plist)
{
  // Set up parameters for the thermal conductivity model
  kappa0 = plist.get<double>("Thermal Conductivity at Reference Temperature");
  lambda = plist.get<double>("Linear Thermal Coefficient");
  Tref = plist.get<double>("Reference Temperature");
  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
ThermalConduct_LinearTempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Thermal Conductivity ParameterList", false, "");
  p->sublist("Thermal Conductivity ParameterList").set<std::string>("Value", "LinearTempDep", "Linear temperature-dependent thermal conductivity");
  p->sublist("Thermal Conductivity ParameterList").set<double>("Thermal Conductivity at Reference Temperature", 0., "[W/(cm.K)]");
  p->sublist("Thermal Conductivity ParameterList").set<double>("Linear Thermal Coefficient", 0., "[1/K]");
  p->sublist("Thermal Conductivity ParameterList").set<double>("Reference Temperature", 0., "[K]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
