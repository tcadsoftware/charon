
#ifndef CHARON_MOBILITY_DEFAULT_IMPL_HPP
#define CHARON_MOBILITY_DEFAULT_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Default<EvalT, Traits>::
Mobility_Default(
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

  // retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);
  num_edges = num_points;   // initialization

  // input from closure model factory
  mobValue = p.get<double>("Value");
  string carrType = p.get<string>("Carrier Type");
  isEdgedl = p.get<bool>("Is Edge Data Layout");

  // retrieve edge data layout if isEdgedl = true
  RCP<DataLayout> output_dl = scalar;
  if (isEdgedl)
  {
    output_dl = p.get< RCP<DataLayout> >("Edge Data Layout");
    num_edges = output_dl->dimension(1);
  }

  // set evaluated field
  if (carrType == "Electron")
    mobility = MDField<ScalarT,Cell,Point>(n.field.elec_mobility,output_dl);
  else if (carrType == "Hole")
    mobility = MDField<ScalarT,Cell,Point>(n.field.hole_mobility,output_dl);
  else if (carrType == "Ion")
    mobility = MDField<ScalarT,Cell,Point>(n.field.ion_mobility,output_dl);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole, or Ion !");

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;

  // evaluated field
  this->addEvaluatedField(mobility);

  std::string name = "Mobility_Default";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Default<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  ScalarT mobScaled = mobValue / Mu0;

  if (isEdgedl)
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int edge = 0; edge < num_edges; ++edge)
        mobility(cell,edge) = mobScaled;
  }
  else
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int point = 0; point < num_points; ++point)
        mobility(cell,point) = mobScaled;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Default<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<double>("Value", 0.0);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);
  p->set("Edge Data Layout", dl);

  p->set("Is Edge Data Layout", false);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
