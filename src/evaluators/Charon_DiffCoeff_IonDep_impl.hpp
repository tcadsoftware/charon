
#ifndef CHARON_DIFFCOEFF_IONDEP_IMPL_HPP
#define CHARON_DIFFCOEFF_IONDEP_IMPL_HPP

#include <cmath>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

/*
The model implements the ion diff. coeff. as
Dion/uion = kbT/q * (1 + c/(1-c)) in physical units, where c = Nion/Nmax.
Nmax = Maximum Ion Density in [cm^(-3)]. In scaled units, Ds/us = Ts * (1+c/(1-c)).

Since c can become <= 0 or >= 1 during simulation, to obtain meaningful value
for Ds under these conditions, we have
Ds = Ts * us for c <= 0,
Ds = Ts * us / ((1-c) + 1/factor ) for 0 < c < 1, (for Reciprocal)
Ds = Ts * us * factor for c >= 1.
factor is given through Maximum Multiply Factor.

If AD Function Type = ReciprocalSqrt, Ds = Ts * us / (sqrt(1-c) + 1/factor)
for 0 < c < 1

An example of using the model is given below:
            <ParameterList name="Ion Diffusion Coefficient">
                <Parameter name="Value" type="string" value="IonDep"/>
                <Parameter name="Maximum Ion Density" type="double" value="5e21" />
                <Parameter name="Maximum Multiply Factor" type="double" value="1e5" />
                <Parameter name="AD Function Type" type="string" value="Reciprocal" />
            </ParameterList>

This model is applicable to ions/vacancies only, neither electrons nor holes.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DiffCoeff_IonDep<EvalT, Traits>::
DiffCoeff_IonDep(
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

  // retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // input ParameterList
  const ParameterList& plist = p.sublist("Diffusion ParameterList");
  maxIonDens = plist.get<double>("Maximum Ion Density");  // [cm^(-3)]
  maxFactor = plist.get<double>("Maximum Multiply Factor");

  funcType = "Reciprocal"; // default
  if (plist.isParameter("AD Function Type"))
    funcType = plist.get<std::string>("AD Function Type");
  TEUCHOS_ASSERT( (funcType == "Reciprocal") || (funcType == "ReciprocalSqrt") )

  // evaluated field
  diffcoeff = MDField<ScalarT,Cell,Point>(n.field.ion_diff_coeff,scalar);
  this->addEvaluatedField(diffcoeff);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  C0 = scaleParams->scale_params.C0;

  // dependent fields
  mobility = MDField<const ScalarT,Cell,Point>(n.field.ion_mobility,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  carr_dens = MDField<const ScalarT,Cell,Point>(n.dof.iondensity,scalar);

  this->addDependentField(mobility);
  this->addDependentField(latt_temp);
  this->addDependentField(carr_dens);

  std::string name = "Diffusion_Coefficient_IonDep";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DiffCoeff_IonDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  //----------------------------------------------------------------------------
  // evaluate diff. ceoff. at IP or BASIS points
  //----------------------------------------------------------------------------

  ScalarT scaledMaxIonDens = maxIonDens / C0;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT mob = mobility(cell,point);
      ScalarT latt = latt_temp(cell,point);  // scaled

      // latt should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0/T0;

      // use the Einstein D/u relation
      ScalarT diff = mob * latt;  // scaled

      // get ion density
      ScalarT iondens = carr_dens(cell,point);
      ScalarT ratio = iondens / scaledMaxIonDens;

      // note that ratio could be negative since iondens can be negative during solving
      if (Sacado::ScalarValue<ScalarT>::eval(ratio) <= 0.0)
        diffcoeff(cell,point) = diff;

      else if (Sacado::ScalarValue<ScalarT>::eval(ratio) >= 1.0)
        diffcoeff(cell,point) = diff * maxFactor;

      else  // 0. < ratio < 1.
      {
        if (funcType == "Reciprocal")
          diffcoeff(cell,point) = diff / (1.0-ratio + 1.0/maxFactor);
        else if (funcType == "ReciprocalSqrt")
          diffcoeff(cell,point) = diff / (std::sqrt(1.0-ratio) + 1.0/maxFactor);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DiffCoeff_IonDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Diffusion ParameterList", false, "");
  p->sublist("Diffusion ParameterList").set<std::string>("Value", "IonDep", "Ion density dependent ion diffusion coefficient");
  p->sublist("Diffusion ParameterList").set<double>("Maximum Ion Density", 0., "[cm^(-3)]");
  p->sublist("Diffusion ParameterList").set<double>("Maximum Multiply Factor", 1., "[unitless]");
  p->sublist("Diffusion ParameterList").set<std::string>("AD Function Type", "Reciprocal", "Reciprocal or ReciprocalSqrt");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
