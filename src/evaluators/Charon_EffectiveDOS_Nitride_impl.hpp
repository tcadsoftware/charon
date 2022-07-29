#ifndef CHARON_EFFECTIVEDOS_NITRIDE_IMPL_HPP
#define CHARON_EFFECTIVEDOS_NITRIDE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
 * Computes Nitride Density of States as a function of mole fraction
 * Vanheusden et al. Nature Vol386 (1997):587-589
 */

namespace charon {

//**********************************************************************
PHX_EVALUATOR_CTOR(EffectiveDOS_Nitride,p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  auto scalar = p.get< RCP<DataLayout> >("Data Layout"); 
  num_points = scalar->dimension(1);

  // Material name
  materialName = p.get<string>("Material Name");
  
  // Evaluated fields
  elec_effdos = MDField<ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
  hole_effdos = MDField<ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);

  this->addEvaluatedField(elec_effdos);
  this->addEvaluatedField(hole_effdos);
  
  // Dependent fields
  latttemp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  molefrac = MDField<const ScalarT,Cell,Point>(n.field.mole_frac,scalar);

  this->addDependentField(latttemp);
  this->addDependentField(molefrac);
 
  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  std::string name = "Effective_DOS_Nitride";
  this->setName(name);
}


//**********************************************************************
PHX_POST_REGISTRATION_SETUP(EffectiveDOS_Nitride,sd,fm)
{
  this->utils.setFieldData(elec_effdos,fm);
  this->utils.setFieldData(hole_effdos,fm);
  
  this->utils.setFieldData(latttemp,fm);
  this->utils.setFieldData(molefrac,fm);
}


//**********************************************************************
PHX_EVALUATE_FIELDS(EffectiveDOS_Nitride,workset)
{
  using panzer::index_t;
  using std::pow;

  // Obtain physical constants
  auto const& cpc = charon::PhysicalConstants::Instance();
  auto const& kb = cpc.kb;   // [eV/K]
  auto const& m0 = cpc.m0; 
  auto const& q = cpc.q;   
  auto const& pi = cpc.pi;
  auto const& h = cpc.h;

  ScalarT me{0.0}, mh{0.0};

  // loop over the cells 
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature [K]
      ScalarT lattT = latttemp(cell,point)*T0; 
      ScalarT x = molefrac(cell,point); 
      // obtain Density of States Masses
      if (materialName == "AlGaN")
      {
        me = 0.314*x + 0.2*(1-x);
        mh = 0.417*x + 1.0*(1-x);
      }

      if (materialName == "InGaN")
      {
        me = 0.12*x + 0.2*(1-x);
        mh = 0.17*x + 1.0*(1-x);
      }

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;  
      
      // calculate the effective density of states
      ScalarT Nc = pow((2*pi*m0*me*kb*q*lattT)/(pow(h,2)),1.5)*1e-6*2.0;
      ScalarT Nv = pow((2*pi*m0*mh*kb*q*lattT)/(pow(h,2)),1.5)*1e-6*2.0;
  
      elec_effdos(cell,point) = Nc / C0;  // scaled
      hole_effdos(cell,point) = Nv / C0; 
    }
  }  

}


//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
EffectiveDOS_Nitride<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Effective DOS ParameterList", false, "");
  p->sublist("Effective DOS ParameterList").set<std::string>("Value", "Nitride", "");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

//**********************************************************************

}

#endif
