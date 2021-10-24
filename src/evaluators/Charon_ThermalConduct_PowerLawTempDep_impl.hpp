
#ifndef CHARON_THERMALCONDUCT_POWERLAWTEMPDEP_IMPL_HPP
#define CHARON_THERMALCONDUCT_POWERLAWTEMPDEP_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
ThermalConduct_PowerLawTempDep<EvalT, Traits>::
ThermalConduct_PowerLawTempDep(
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

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Thermal Conductivity ParameterList
  const ParameterList& plist = p.sublist("Thermal Conductivity ParameterList");

  // Initialize the thermal conductivity model parameters
  initialize(matName, plist);

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

  std::string name = "ThermalConductivity_PowerLawTempDep";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
ThermalConduct_PowerLawTempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT TL = T0 * latt_temp(cell,point);  // in [K]
      ScalarT tc;

      // lattice temperature should be always > 0
      if (Sacado::ScalarValue<ScalarT>::eval(TL) > std::numeric_limits<double>::epsilon())
        tc = kappa300 * std::pow(TL/300., alpha);    // in [W/(K.cm)]
      else
        tc = kappa300;

      ther_cond(cell,point) = tc / kL0;  // scaled
      // std::cout << "cell=" << cell << ", point=" << point << ", ther_cond=" << ther_cond(cell,point) << std::endl;
    }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void ThermalConduct_PowerLawTempDep<EvalT, Traits>::initialize
(const std::string& matName, const Teuchos::ParameterList& plist)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for the thermal conductivity model
  // Retrieve parameters from charon::Material_Properties by default
  if (plist.isParameter("kappa300"))
    kappa300 = plist.get<double>("kappa300");
  else
    kappa300 = matProperty.getPropertyValue(matName, "Thermal Conductivity kappa300");

  if (plist.isParameter("alpha"))
    alpha = plist.get<double>("alpha");
  else
    alpha = matProperty.getPropertyValue(matName, "Thermal Conductivity alpha");

  // std::cout << "Material=" << matName << ", kappa300=" << kappa300 << ", alpha=" << alpha << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
ThermalConduct_PowerLawTempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Thermal Conductivity ParameterList", false, "");
  p->sublist("Thermal Conductivity ParameterList").set<std::string>("Value", "PowerLawTempDep", "Temperature-dependent thermal conductivity");
  p->sublist("Thermal Conductivity ParameterList").set<double>("kappa300", 0., "kappa300:[W/(K.cm)]");
  p->sublist("Thermal Conductivity ParameterList").set<double>("alpha",    0., "alpha:[1]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
