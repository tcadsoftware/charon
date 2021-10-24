
#ifndef CHARON_HEATCAPACITY_POWERLAWTEMPDEP_IMPL_HPP
#define CHARON_HEATCAPACITY_POWERLAWTEMPDEP_IMPL_HPP

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
HeatCapacity_PowerLawTempDep<EvalT, Traits>::
HeatCapacity_PowerLawTempDep(
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

  // Heat Capacity ParameterList
  const ParameterList& plist = p.sublist("Heat Capacity ParameterList");

  // Initialize the heat capacity model parameters
  initialize(matName, plist);

  // Evaluated fields
  heat_cap = MDField<ScalarT,Cell,Point>(n.field.heat_cap,scalar);
  this->addEvaluatedField(heat_cap);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  cL0 = scaleParams->scale_params.cL0;
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,scalar);

  this->addDependentField(latt_temp);

  std::string name = "HeatCapacity_PowerLawTempDep";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
HeatCapacity_PowerLawTempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT TL = T0 * latt_temp(cell,point);  // in [K]
      ScalarT hc;

      // lattice temperature should be always > 0
      if (Sacado::ScalarValue<ScalarT>::eval(TL) > std::numeric_limits<double>::epsilon())
      {
        ScalarT num = std::pow(TL/300., beta) - 1.;
        ScalarT den = std::pow(TL/300., beta) + c1/c300;
        hc = rho * (c300 + c1 * num/den);    // in [J/(K.cm^3)]
      }
      else
        hc = rho * c300;  // in [J/(K.cm^3)]

      heat_cap(cell,point) = hc / cL0;  // scaled
      // std::cout << "cell=" << cell << ", point=" << point << ", heat_cap=" << heat_cap(cell,point) << std::endl;
    }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void HeatCapacity_PowerLawTempDep<EvalT, Traits>::initialize
(const std::string& matName, const Teuchos::ParameterList& plist)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for the heat capacity model
  // Retrieve parameters from charon::Material_Properties by default
  if (plist.isParameter("Mass Density"))
    rho = plist.get<double>("Mass Density");
  else
    rho = matProperty.getPropertyValue(matName, "Mass Density");

  if (plist.isParameter("c300"))
    c300 = plist.get<double>("c300");
  else
    c300 = matProperty.getPropertyValue(matName, "Heat Capacity c300");

  if (plist.isParameter("c1"))
    c1 = plist.get<double>("c1");
  else
    c1 = matProperty.getPropertyValue(matName, "Heat Capacity c1");

  if (plist.isParameter("beta"))
    beta = plist.get<double>("beta");
  else
    beta = matProperty.getPropertyValue(matName, "Heat Capacity beta");

  // std::cout << "rho=" << rho << ", c300=" << c300 << ", c1=" << c1 << ", beta=" << beta << std::endl;

  return;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
HeatCapacity_PowerLawTempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Heat Capacity ParameterList", false, "");
  p->sublist("Heat Capacity ParameterList").set<std::string>("Value", "PowerLawTempDep", "Temperature-dependent heat capacity");
  p->sublist("Heat Capacity ParameterList").set<double>("Mass Density", 0.,  "[g/cm^3]"  );
  p->sublist("Heat Capacity ParameterList").set<double>("c300",         0.,  "[J/(K.g)]" );
  p->sublist("Heat Capacity ParameterList").set<double>("c1",           0.,  "[J/(K.g)]" );
  p->sublist("Heat Capacity ParameterList").set<double>("beta",         0.,  "[1]" );

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
