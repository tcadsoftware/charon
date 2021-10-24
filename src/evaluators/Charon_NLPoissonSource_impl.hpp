
#ifndef CHARON_NLPOISSONSOURCE_IMPL_HPP
#define CHARON_NLPOISSONSOURCE_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"

/*
Evaluate the equilibrium nonlinear Poisson source: nie*exp(Ei/kbT)-nie*exp(-Ei/kbT)+dop
(Ef = 0 at equilibrium), applicable to MB & FD statistics and homogeneous & heterogeneous
devices.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
NLPoissonSource<EvalT, Traits>::
NLPoissonSource(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const RCP<const charon::Names> n = p.get< RCP<const charon::Names> >("Names");

  // Retrieve the data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Define fields
  nlpsrc = MDField<ScalarT,Cell,Point>(p.get<string>("Source Name"),scalar);
  intrin_fermi = MDField<const ScalarT,Cell,Point>(n->field.intrin_fermi,scalar);
  doping = MDField<const ScalarT,Cell,Point>(n->field.doping,scalar);
  intrin_conc = MDField<const ScalarT,Cell,Point>(n->field.intrin_conc,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n->field.latt_temp,scalar);

  // Add fields
  this->addEvaluatedField(nlpsrc);
  this->addDependentField(intrin_fermi);
  this->addDependentField(doping);
  this->addDependentField(intrin_conc);
  this->addDependentField(latt_temp);

  std::string name = "NLPoissonSource";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
NLPoissonSource<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature [K]
      ScalarT lattT = latt_temp(cell,point) * T0;
      ScalarT kbT = kbBoltz*lattT;  // [eV]

      const ScalarT& Ei = intrin_fermi(cell,point);  // [eV]
      const ScalarT& dop = doping(cell,point);       // scaled
      const ScalarT& nie = intrin_conc(cell,point);  // scaled

      // nlpsrc(cell,point) = nie*exp(-phi) - nie*exp(phi) + dop;
      nlpsrc(cell,point) = nie*exp(Ei/kbT) - nie*exp(-Ei/kbT) + dop;
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
NLPoissonSource<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Source Name", "?");

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

