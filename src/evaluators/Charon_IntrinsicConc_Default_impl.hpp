
#ifndef CHARON_INTRINSICCONC_DEFAULT_IMPL_HPP
#define CHARON_INTRINSICCONC_DEFAULT_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IntrinsicConc_Default<EvalT, Traits>::
IntrinsicConc_Default(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Input from closure model
  niValue = p.get<double>("Value");

  // Evaluated fields
  ni = MDField<ScalarT,Cell,Point>(n.field.intrin_conc, scalar);
  effEg = MDField<ScalarT,Cell,Point>(n.field.eff_band_gap, scalar);
  effChi = MDField<ScalarT,Cell,Point>(n.field.eff_affinity, scalar);

  this->addEvaluatedField(ni);
  this->addEvaluatedField(effEg);
  this->addEvaluatedField(effChi);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;

  // Dependent fields
  Eg = MDField<const ScalarT,Cell,Point>(n.field.band_gap, scalar);
  Chi = MDField<const ScalarT,Cell,Point>(n.field.affinity, scalar);

  this->addDependentField(Eg);
  this->addDependentField(Chi);

  std::string name = "Intrinsic_Concentration_Default";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IntrinsicConc_Default<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // scaled intrinsic concentration
  ScalarT niScaled = niValue / C0;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ni(cell,point) = niScaled;
      effEg(cell,point) = Eg(cell,point);  // [eV], Eg in fm is NOT scaled
      effChi(cell,point) = Chi(cell,point); // [eV]
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
IntrinsicConc_Default<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<double>("Value", 0.);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
