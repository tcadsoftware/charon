
#ifndef CHARON_INTRINSIC_FERMIENERGY_IMPL_HPP
#define CHARON_INTRINSIC_FERMIENERGY_IMPL_HPP

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
Intrinsic_FermiEnergy<EvalT, Traits>::
Intrinsic_FermiEnergy(
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
  intrin_fermi = MDField<ScalarT,Cell,Point>(n.field.intrin_fermi,scalar);
  this->addEvaluatedField(intrin_fermi);

  // Dependent fields
  potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,scalar);
  this->addDependentField(potential);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;

    T0 = scaleParams->scale_params.T0;
    eff_affinity = MDField<const ScalarT,Cell,Point>(n.field.eff_affinity,scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,scalar);
    latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);
    elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
    hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);
    ref_energy = MDField<const ScalarT,Cell,Point>(n.field.ref_energy, scalar);
    elec_degfactor = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor,scalar);
    hole_degfactor = MDField<const ScalarT,Cell,Point>(n.field.hole_deg_factor,scalar);

    this->addDependentField(eff_affinity);
    this->addDependentField(eff_bandgap);
    this->addDependentField(latt_temp);
    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
    this->addDependentField(ref_energy);
    this->addDependentField(elec_degfactor);
    this->addDependentField(hole_degfactor);

  std::string name = "Intrinsic_FermiEnergy";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Intrinsic_FermiEnergy<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain kb in [eV/K]
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // Reference Energy
  ScalarT Eref = ref_energy(0,0);

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // Lattice temperature in [K]
      ScalarT kbT = kbBoltz * latt_temp(cell,point) * T0;  // [eV]

      // Effective electron affinity and band gap in [eV]
      const ScalarT& Chieff = eff_affinity(cell, point);
      const ScalarT& Egeff = eff_bandgap(cell, point);

      // Obtain phi in [V]
      ScalarT phi = potential(cell,point) * V0;  // unscaled in [V]

      // Obtain the effective density of states (scaled)
      const ScalarT& Nc = elec_effdos(cell, point);
      const ScalarT& Nv = hole_effdos(cell, point);

      // Obtain the degeneracy factor [unitless]
      const ScalarT& gamma_n = elec_degfactor(cell,point);
      const ScalarT& gamma_p = hole_degfactor(cell,point);

      // Intrinsic fermi energy in [eV]
      intrin_fermi(cell,point) = Eref - Chieff - 1.0*phi - 0.5*Egeff
            - 0.5*kbT*std::log(Nc/Nv) - 0.5*kbT*std::log(gamma_n/gamma_p);
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
Intrinsic_FermiEnergy<EvalT, Traits>::getValidParameters() const
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
