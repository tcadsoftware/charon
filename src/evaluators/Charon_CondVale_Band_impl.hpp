
#ifndef CHARON_CONDVALE_BAND_IMPL_HPP
#define CHARON_CONDVALE_BAND_IMPL_HPP

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
CondVale_Band<EvalT, Traits>::
CondVale_Band(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Evaluated fields
  cond_band = MDField<ScalarT,Cell,Point>(n.field.cond_band,scalar);
  vale_band = MDField<ScalarT,Cell,Point>(n.field.vale_band,scalar);
  this->addEvaluatedField(cond_band);
  this->addEvaluatedField(vale_band);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;

  // Dependent fields
  eff_affinity = MDField<const ScalarT,Cell,Point>(n.field.eff_affinity,scalar);
  eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,scalar);
  potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,scalar);

  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(potential);

  ref_energy = MDField<const ScalarT,Cell,Point>(n.field.ref_energy, scalar);
  this->addDependentField(ref_energy);

  std::string name = "CondVale_Band";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
CondVale_Band<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Reference Energy
  ScalarT Eref = ref_energy(0,0);
  
  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Effective electron affinity and band gap in [eV]
      const ScalarT& Chieff = eff_affinity(cell, point);
      const ScalarT& Egeff = eff_bandgap(cell, point);
      
      // Obtain phi in [V]
      ScalarT phi = potential(cell,point) * V0;  // unscaled in [V]
      
      // Conduction band energy in [eV]
      cond_band(cell,point) = Eref - Chieff - 1.0*phi;  // 1.0 indicates 1.0 [q]
      
      // Valence band energy in [eV]
      vale_band(cell,point) = cond_band(cell,point) - Egeff;
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
CondVale_Band<EvalT, Traits>::getValidParameters() const
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
