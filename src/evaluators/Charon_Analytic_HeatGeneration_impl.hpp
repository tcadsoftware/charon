
#ifndef CHARON_ANALYTIC_HEATGENERATION_IMPL_HPP
#define CHARON_ANALYTIC_HEATGENERATION_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


/*
Analytic expressions for the heat generation. Specification of the analytic
heat generation model in the input file takes the form of

<ParameterList name="Heat Generation">
  <Parameter name="Value" type="string" value="Analytic"/>
  <Parameter name="Type" type="string" value="Constant"/>
  <Parameter name="Constant Value" type="double" value="0"/>
  <!--
  <Parameter name="Type" type="string" value="Linear"/>
  <Parameter name="Linear Factor" type="double" value="1"/>
  !-->
</ParameterList>

According to Type, we can code any analytic expressions as needed. Presently,
we have Type = Constant (H = const) or Linear (H = factor*TL).
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Analytic_HeatGeneration<EvalT, Traits>::
Analytic_HeatGeneration(
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

  // Heat Generation ParameterList
  const ParameterList& plist = p.sublist("Heat Generation ParameterList");
  heatGenType = plist.get<string>("Type");
  if (heatGenType == "Constant")
    constValue = plist.get<double>("Constant Value");
  else if (heatGenType == "Linear")
    linFactor = plist.get<double>("Linear Factor");
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: wrong Type for analytic heat generation !");

  // Evaluated fields
  heat_gen = MDField<ScalarT,Cell,Point>(n.field.heat_gen,scalar);
  this->addEvaluatedField(heat_gen);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  H0 = scaleParams->scale_params.H0;

  // Dependent fields
  if (heatGenType == "Linear")
  {
    T0 = scaleParams->scale_params.T0;
    latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,scalar);
    this->addDependentField(latt_temp);
  }

  std::string name = "Analytic_HeatGeneration";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Analytic_HeatGeneration<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;
  if (heatGenType == "Constant")
  {
    ScalarT scaledValue = constValue/H0;
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int point = 0; point < num_points; ++point)
        heat_gen(cell,point) = scaledValue;  // scaled
  }
  else if (heatGenType == "Linear")
  {
    ScalarT factor = linFactor*T0/H0;
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int point = 0; point < num_points; ++point)
        heat_gen(cell,point) = latt_temp(cell,point)*factor;  // scaled
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Analytic_HeatGeneration<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Heat Generation ParameterList", false, "");
  p->sublist("Heat Generation ParameterList").set<std::string>("Value", "Analytic", "Analytic heat generation");
  p->sublist("Heat Generation ParameterList").set<std::string>("Type", "Constant", "Constant or Linear");
  p->sublist("Heat Generation ParameterList").set<double>("Constant Value", 0., "[W/cm^3]");
  p->sublist("Heat Generation ParameterList").set<double>("Linear Factor", 0., "[W/(cm^3.K)]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
