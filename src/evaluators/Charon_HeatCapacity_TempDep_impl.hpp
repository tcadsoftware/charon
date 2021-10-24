
#ifndef CHARON_HEATCAPACITY_TEMPDEP_IMPL_HPP
#define CHARON_HEATCAPACITY_TEMPDEP_IMPL_HPP

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
HeatCapacity_TempDep<EvalT, Traits>::
HeatCapacity_TempDep(
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

  std::string name = "HeatCapacity_TempDep";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
HeatCapacity_TempDep<EvalT, Traits>::
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
      ScalarT hc;

      // lattice temperature should be always > 0
      if (Sacado::ScalarValue<ScalarT>::eval(TL) > std::numeric_limits<double>::epsilon())
        hc = a + b*TL + c*TL*TL;    // in [J/(K.cm^3)]
      else
        hc = a;

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
void HeatCapacity_TempDep<EvalT, Traits>::initialize
(const std::string& matName, const Teuchos::ParameterList& plist)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for the heat capacity model
  // Retrieve parameters from charon::Material_Properties by default
  if (plist.isParameter("a"))
    a = plist.get<double>("a");
  else
    a = matProperty.getPropertyValue(matName, "Heat Capacity a");

  if (plist.isParameter("b"))
    b = plist.get<double>("b");
  else
    b = matProperty.getPropertyValue(matName, "Heat Capacity b");

  if (plist.isParameter("c"))
    c = plist.get<double>("c");
  else
    c = matProperty.getPropertyValue(matName, "Heat Capacity c");
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
HeatCapacity_TempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Heat Capacity ParameterList", false, "");
  p->sublist("Heat Capacity ParameterList").set<std::string>("Value", "TempDep", "Temperature-dependent heat capacity");
  p->sublist("Heat Capacity ParameterList").set<double>("a", 0., "[J/(K.cm^3)]");
  p->sublist("Heat Capacity ParameterList").set<double>("b", 0., "[J/(K^2.cm^3)]");
  p->sublist("Heat Capacity ParameterList").set<double>("c", 0., "[J/(K^3.cm^3)]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
