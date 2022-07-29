#ifndef CHARON_QUASIFERMILEVELS_IMPL_HPP
#define CHARON_QUASIFERMILEVELS_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"


const double n_lim = 1e-20;
const double p_lim = 1e-20;

namespace charon {

template<typename EvalT, typename Traits>
QuasiFermiLevels<EvalT, Traits>::
QuasiFermiLevels(
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

  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  C0 = scaleParams->scale_params.C0;

  // Evaluated fields
  eQF = MDField<ScalarT,Cell,Point>(n.field.eQF,scalar);
  hQF = MDField<ScalarT,Cell,Point>(n.field.hQF,scalar);
  this->addEvaluatedField(eQF);
  this->addEvaluatedField(hQF);

  // Dependent fields
  edens = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hdens = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
  lT = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  Ec = MDField<const ScalarT,Cell,Point>(n.field.cond_band,scalar);
  Ev = MDField<const ScalarT,Cell,Point>(n.field.vale_band,scalar);
  Nc = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
  Nv = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);
  e_gamma = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor,scalar);
  h_gamma = MDField<const ScalarT,Cell,Point>(n.field.hole_deg_factor,scalar);
  this->addDependentField(edens);
  this->addDependentField(hdens);
  this->addDependentField(lT);
  this->addDependentField(Ec);
  this->addDependentField(Ev);
  this->addDependentField(Nc);
  this->addDependentField(Nv);
  this->addDependentField(e_gamma);
  this->addDependentField(h_gamma);

  std::string name = "QuasiFermiLevels";
  this->setName(name);

}



template<typename EvalT, typename Traits>
void
QuasiFermiLevels<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain kb in [eV/K]
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb; // Boltzmann constant in [eV/K]

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < num_points; ++point) {
      ScalarT kbT = kbBoltz * lT(cell,point) * T0;  // [eV]
      ScalarT n = edens(cell,point);
      if (n < n_lim/C0) n = n_lim/C0;
      ScalarT p = hdens(cell,point);
      if (p < p_lim/C0) p = p_lim/C0;
      ScalarT log_edens = std::log(e_gamma(cell,point)*Nc(cell,point)/n);
      ScalarT log_hdens = std::log(h_gamma(cell,point)*Nv(cell,point)/p);
      eQF(cell,point) = Ec(cell,point) - kbT*log_edens; // [eV]
      hQF(cell,point) = Ev(cell,point) + kbT*log_hdens; // [eV]
    }
  } // cells

}



template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
QuasiFermiLevels<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

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
