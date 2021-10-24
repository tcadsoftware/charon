
#ifndef CHARON_DDLATTICE_VACUUMPOTENTIAL_IMPL_HPP
#define CHARON_DDLATTICE_VACUUMPOTENTIAL_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"


/*
Evaluate the vaccum potential for the DDLattice, DDIon, and DDIonLattice equation sets.
vac_pot = \phi-\theta, where \phi is the electric potential DOF (scaled),
and \theta = (\chi + 0.5*Eg)/kbT0 + 0.5*T*log(Nc/Nv), with all quantities in
scaled units, except that \chi and Eg are in [eV].
\phi corresponds to the intrinsic Fermi potential, and \theta is the so-called
band structure parameter. 
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLattice_VacuumPotential<EvalT, Traits>::
DDLattice_VacuumPotential(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Get data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Evaluated field
  vac_pot = MDField<ScalarT,Cell,Point>(n.field.vac_pot, scalar);
  this->addEvaluatedField(vac_pot);

  // Scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  potential = MDField<const ScalarT,Cell,Point>(n.dof.phi, scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp, scalar);
  eff_affinity = MDField<const ScalarT,Cell,Point>(n.field.eff_affinity, scalar);
  eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos, scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos, scalar);

  this->addDependentField(potential);
  this->addDependentField(latt_temp);
  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);

  std::string name = "DDLattice_VacuumPotential";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLattice_VacuumPotential<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb*T0 in [eV]
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;       // Boltzmann constant in [eV/K]
  ScalarT kbT0 = kbBoltz*T0;  // [eV]

  // compute vacuum potential
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t point = 0; point < num_points; ++point)
    {
      // obtain the effective electron affinity and band gap in [eV]
      const ScalarT& chi = eff_affinity(cell,point);
      const ScalarT& Eg = eff_bandgap(cell,point);

      // obtain the conduction and valence band effective dos [scaled]
      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);

      // obtain the lattice temperature [scaled]
      ScalarT TL = latt_temp(cell,point);

      // TL should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(TL) <= 0.0)  TL = 300.0/T0;

      // compute theta [scaled]
      ScalarT theta = (chi + 0.5*Eg)/kbT0 + 0.5*TL*log(Nc/Nv);

      // compute vacuum potential [scaled]
      vac_pot(cell,point) = potential(cell,point) - theta;
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
DDLattice_VacuumPotential<EvalT, Traits>::getValidParameters() const
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

