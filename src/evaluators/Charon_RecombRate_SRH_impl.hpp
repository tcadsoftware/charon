
#ifndef CHARON_RECOMBRATE_SRH_IMPL_HPP
#define CHARON_RECOMBRATE_SRH_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

namespace charon {


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFDIntrinsicDensity()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename EvalT::ScalarT FermiDiracIntrinsicDensity<EvalT, Traits>::evaluateFDIntrinsicDensity(
         const typename EvalT::ScalarT& n,
         const typename EvalT::ScalarT& p,
         const typename EvalT::ScalarT& niMB,
         const typename EvalT::ScalarT& Nc,
         const typename EvalT::ScalarT& Nv,
         const typename EvalT::ScalarT& Eg,
         const typename EvalT::ScalarT& kbT,
         const Teuchos::RCP<charon::FermiDiracIntegral<EvalT> >& invFDInt)
{
  typedef typename EvalT::ScalarT ScalarT;
  ScalarT n0p0;

  // n-type region
  if (n > p)
  {
    // n-n0=dn=dp=p-p0 = constant for SRH, Auger, and Radiative processes since they require e-h pairs
    ScalarT a = n-p;  // [cm^-3]

    // first solving for n0 by assuming n0*p0 = n0*(n0-a) = niMB*niMB
    ScalarT n0 = (a + std::sqrt(a*a + 4.0*niMB*niMB)) / 2.0;

    // obtain the equilibrium (Ef-Ec)
    ScalarT Ef_minus_Ec = kbT * (*invFDInt)(n0/Nc);  // [eV]

    // compute the equilibrium (Ef-Ev)
    ScalarT Ef_minus_Ev = Ef_minus_Ec + Eg;  // [eV]

    // compute the equilibrium hole concentration
    ScalarT p0 = Nv*exp(-Ef_minus_Ev/kbT);

    // compute the equilibrium electron conc. with only the n-p = n0-p0 condition
    n0 = a + p0;

    n0p0 = n0*p0;  // [cm^-6]
  }

  // p-type region
  else
  {
    // a = p-n = p0-n0 > 0
    ScalarT a = p-n;

    // first solving for p0 by assuming n0*p0 = (p0-a)*p0 = niMB*niMB
    ScalarT p0 = (a + std::sqrt(a*a + 4.0*niMB*niMB)) / 2.0;

    // obtain the equilibrium (Ev-Ef)
    ScalarT Ev_minus_Ef = kbT * (*invFDInt)(p0/Nv);  // [eV]

    // compute the equilibrium (Ec-Ef)
    ScalarT Ec_minus_Ef = Ev_minus_Ef + Eg;  // [eV]

    // compute the equilibrium electron concentration
    ScalarT n0 = Nc*exp(-Ec_minus_Ef/kbT);

    // compute the equilibrium hole conc. with only the p-n = p0-n0 condition
    p0 = a + n0;

    n0p0 = n0*p0;  // [cm^-6]
  }

  return n0p0;
}


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_SRH<EvalT, Traits>::
RecombRate_SRH(
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

  bUseFD = p.get<bool>("Fermi Dirac");

  // Evaluated fields
  srh_rate = MDField<ScalarT,Cell,Point>(n.field.srh_recomb,scalar);
  srh_deriv_e = MDField<ScalarT,Cell,Point>(n.field.srh_deriv_e,scalar);
  srh_deriv_h = MDField<ScalarT,Cell,Point>(n.field.srh_deriv_h,scalar);

  this->addEvaluatedField(srh_rate);
  this->addEvaluatedField(srh_deriv_e);
  this->addEvaluatedField(srh_deriv_h);

  // Dependent fields
  elifetime = MDField<const ScalarT,Cell,Point>(n.field.elec_lifetime,scalar);
  hlifetime = MDField<const ScalarT,Cell,Point>(n.field.hole_lifetime,scalar);
  intrin_conc = MDField<const ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);

  this->addDependentField(elifetime);
  this->addDependentField(hlifetime);
  this->addDependentField(intrin_conc);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");

  if (bUseFD)  // need additional fields when bUseFD = true
  {
    elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos, scalar);
    hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos, scalar);
    eff_bandgap =  MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, scalar);
    latt_temp =  MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);
    C0 = scaleParams->scale_params.C0;
    T0 = scaleParams->scale_params.T0;

    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
    this->addDependentField(eff_bandgap);
    this->addDependentField(latt_temp);
  }

  // unscaled and scaled Rsrh have the same expressions

  std::string name = "SRH_Recombination_Rate";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_SRH<EvalT, Traits>::
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
      const ScalarT& n = edensity(cell,point);      // scaled
      const ScalarT& p = hdensity(cell,point);

      if ((Sacado::ScalarValue<ScalarT>::eval(n) > 0.0) &&
          (Sacado::ScalarValue<ScalarT>::eval(p) > 0.0))  // assure positive densities
      {
        const ScalarT& ni = intrin_conc(cell,point);  // scaled
        const ScalarT& taun = elifetime(cell,point);
        const ScalarT& taup = hlifetime(cell,point);
        const ScalarT& Nc = elec_effdos(cell,point);  // scaled
        const ScalarT& Nv = hole_effdos(cell,point);

        const ScalarT& Eg = eff_bandgap(cell,point);  // [eV]
        ScalarT lattT = latt_temp(cell,point) * T0;  // [K]

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

        ScalarT tmp1 = n*p - n0p0;
        ScalarT tmp2 = taup*(n+std::sqrt(n0p0)) + taun*(p+std::sqrt(n0p0));

        srh_rate(cell,point) = tmp1 / tmp2;
        srh_deriv_e(cell,point) = p/tmp2 - tmp1*taup/(tmp2*tmp2);
        srh_deriv_h(cell,point) = n/tmp2 - tmp1*taun/(tmp2*tmp2);
      }
      else
      {
        srh_rate(cell,point) = 0.0;
        srh_deriv_e(cell,point) = 0.0;
        srh_deriv_h(cell,point) = 0.0;
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
        const ScalarT& taun = elifetime(cell,point);
        const ScalarT& taup = hlifetime(cell,point);

        ScalarT tmp1 = n*p - ni*ni;
        ScalarT tmp2 = taup*(n+ni) + taun*(p+ni);
        srh_rate(cell,point) = tmp1 / tmp2;
        srh_deriv_e(cell,point) = p/tmp2 - tmp1*taup/(tmp2*tmp2);
        srh_deriv_h(cell,point) = n/tmp2 - tmp1*taun/(tmp2*tmp2);
      }
      else
      {
        srh_rate(cell,point) = 0.0;
        srh_deriv_e(cell,point) = 0.0;
        srh_deriv_h(cell,point) = 0.0;
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
RecombRate_SRH<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

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
