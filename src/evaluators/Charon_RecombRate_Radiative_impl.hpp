
#ifndef CHARON_RECOMBRATE_RADIATIVE_IMPL_HPP
#define CHARON_RECOMBRATE_RADIATIVE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Names.hpp"
#include "Charon_RecombRate_SRH.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_FermiDirac_Integral.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_Radiative<EvalT, Traits>::
RecombRate_Radiative(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Recombination coefficient
  recomb_coeff = p.get<double>("Coefficient");
  bUseFD = p.get<bool>("Fermi Dirac");

  // Evaluated fields
  rad_rate = MDField<ScalarT,Cell,Point>(n.field.rad_recomb,scalar);
  rad_deriv_e = MDField<ScalarT,Cell,Point>(n.field.rad_deriv_e,scalar);
  rad_deriv_h = MDField<ScalarT,Cell,Point>(n.field.rad_deriv_h,scalar);

  this->addEvaluatedField(rad_rate);
  this->addEvaluatedField(rad_deriv_e);
  this->addEvaluatedField(rad_deriv_h);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  R0 = scaleParams->scale_params.R0;

  // Dependent fields
  intrin_conc = MDField<const ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);

  this->addDependentField(intrin_conc);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);

  if (bUseFD)  // need additional fields when bUseFD = true
  {
    elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos, scalar);
    hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos, scalar);
    eff_bandgap =  MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, scalar);
    latt_temp =  MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);
    T0 = scaleParams->scale_params.T0;

    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
    this->addDependentField(eff_bandgap);
    this->addDependentField(latt_temp);
  }

  std::string name = "Radiative_Recombination_Rate";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_Radiative<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using panzer::index_t;

 // Use the Fermi-Dirac statistics
 if (bUseFD)
 {
  // instantiate the FermiDiracIntegral class
  RCP<charon::FermiDiracIntegral<EvalT> > invFDInt =
      rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

  // obtain kb
  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double kbBoltz = phyConst.kb;   // Boltzmann constant in [eV/K]

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      const ScalarT& n = edensity(cell,point);
      const ScalarT& p = hdensity(cell,point);

      if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
          (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
      {
        const ScalarT& ni = intrin_conc(cell,point);
        const ScalarT& Nc = elec_effdos(cell,point);  // scaled
        const ScalarT& Nv = hole_effdos(cell,point);
        const ScalarT& Eg = eff_bandgap(cell,point);  // [eV]
        ScalarT lattT = latt_temp(cell,point) * T0; // [K]

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

        ScalarT n_us = n * C0;  // [cm^-3]
        ScalarT p_us = p * C0;
        ScalarT ni_us = ni * C0;
        ScalarT Nc_us = Nc * C0;
        ScalarT Nv_us = Nv * C0;
        ScalarT kbT = kbBoltz * lattT;  // [eV]

        ScalarT n0p0_us = FermiDiracIntrinsicDensity<EvalT, Traits>::evaluateFDIntrinsicDensity(
                n_us, p_us, ni_us, Nc_us, Nv_us, Eg, kbT, invFDInt);  // [cm^-6]
        ScalarT n0p0 = n0p0_us / C0 / C0;  // scaled

        rad_rate(cell,point) = recomb_coeff *(n*p-n0p0) *C0*C0 /R0;
        rad_deriv_e(cell,point) = recomb_coeff *p *C0*C0 / R0;
        rad_deriv_h(cell,point) = recomb_coeff *n *C0*C0 / R0;
      }
      else
      {
        rad_rate(cell,point) = 0.0;
        rad_deriv_e(cell,point) = 0.0;
        rad_deriv_h(cell,point) = 0.0;
      }
    }
  }
 }  // end of the FD block


 // Use the Maxwell-Boltzmann statistics
 else
 {
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      const ScalarT& n = edensity(cell,point);
      const ScalarT& p = hdensity(cell,point);

      if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
          (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
      {
        const ScalarT& ni = intrin_conc(cell,point);

        rad_rate(cell,point) = recomb_coeff *(n*p-ni*ni) *C0*C0 /R0;
        rad_deriv_e(cell,point) = recomb_coeff *p *C0*C0 / R0;
        rad_deriv_h(cell,point) = recomb_coeff *n *C0*C0 / R0;
      }
      else
      {
        rad_rate(cell,point) = 0.0;
        rad_deriv_e(cell,point) = 0.0;
        rad_deriv_h(cell,point) = 0.0;
      }
    }
  }
 }  // end of the MB block

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
RecombRate_Radiative<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Coefficient", 0.0);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Fermi Dirac", false);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
