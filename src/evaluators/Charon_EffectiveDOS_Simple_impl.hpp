
#ifndef CHARON_EFFECTIVEDOS_SIMPLE_IMPL_HPP
#define CHARON_EFFECTIVEDOS_SIMPLE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
The EffectiveDOS_Simple model computes Nc = Nc300*(T/300)^Nc_F, and Nv = Nv300*
(T/300)^Nv_F, where Nc300, Nv300, Nc_F, and Nv_F can be changed by users.
The scaled version is Nc_scaled = Nc/C0 and Nv_scaled = Nv/C0.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
EffectiveDOS_Simple<EvalT, Traits>::
EffectiveDOS_Simple(
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

  // Material name
  const string& materialName = p.get<string>("Material Name");

  // Obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve material parameters
  Nc300 = matProperty.getPropertyValue(materialName, "Electron Effective DOS at 300 K");
  Nv300 = matProperty.getPropertyValue(materialName, "Hole Effective DOS at 300 K");
  Nc_F = matProperty.getPropertyValue(materialName, "Electron Effective DOS Exponent");
  Nv_F = matProperty.getPropertyValue(materialName, "Hole Effective DOS Exponent");

  // Effective DOS ParameterList
  if (p.isSublist("Effective DOS ParameterList"))
  {
    const ParameterList& effDosParamList = p.sublist("Effective DOS ParameterList");
    // Overwrite parameters when specified by users
    if (effDosParamList.isParameter("Nc300"))
      Nc300 = effDosParamList.get<double>("Nc300");
    if (effDosParamList.isParameter("Nv300"))
      Nv300 = effDosParamList.get<double>("Nv300");
    if (effDosParamList.isParameter("Nc_F"))
      Nc_F = effDosParamList.get<double>("Nc_F");
    if (effDosParamList.isParameter("Nv_F"))
      Nv_F = effDosParamList.get<double>("Nv_F");
  }

  // Evaluated fields
  elec_effdos = MDField<ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
  hole_effdos = MDField<ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);

  this->addEvaluatedField(elec_effdos);
  this->addEvaluatedField(hole_effdos);

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);

  this->addDependentField(latt_temp);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  std::string name = "Effective_DOS_Simple";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
EffectiveDOS_Simple<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // loop over the cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature [K]
      ScalarT lattT = latt_temp(cell,point)*T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      // calculate the effective density of states
      ScalarT Nc = Nc300 * pow(lattT/300.0, Nc_F);
      ScalarT Nv = Nv300 * pow(lattT/300.0, Nv_F);

      elec_effdos(cell,point) = Nc / C0;  // scaled
      hole_effdos(cell,point) = Nv / C0;
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
EffectiveDOS_Simple<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Effective DOS ParameterList", false, "");
  p->sublist("Effective DOS ParameterList").set<std::string>("Value", "Simple", "");
  p->sublist("Effective DOS ParameterList").set<double>("Nc300", 0., "[cm^-3]");
  p->sublist("Effective DOS ParameterList").set<double>("Nv300", 0., "[cm^-3]");
  p->sublist("Effective DOS ParameterList").set<double>("Nc_F", 0., "[1]");
  p->sublist("Effective DOS ParameterList").set<double>("Nv_F", 0., "[1]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
