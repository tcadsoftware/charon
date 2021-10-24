
#ifndef CHARON_DDLATTICE_HEATGENERATION_IMPL_HPP
#define CHARON_DDLATTICE_HEATGENERATION_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


/*
The heat generation is modeled in its full form as
H = J0*E0/H0 * (Jn*Fn + Jp*Fp + R*Eg/(kB*T0) + 3*R*T + Jion*Fion), where all quantities
are in scaled units, except Eg and the scaling parameters of T0, J0, E0, and H0
which are in physical units.

When H0 = J0*E0 (the default case), J0*E0/H0 is not needed as it equals to 1.
When H0 is NOT equal to J0*E0, we need the scaling factor of J0*E0/H0.

When solveElectron = true, include the Jn*Fn term.
When solveHole = true, include the Jp*Fp term.
When solveIon = true, include the Jion*Fion term.
When haveSource = true, include the R*Eg/(kB*T0)+3*R*T term.

*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLattice_HeatGeneration<EvalT, Traits>::
DDLattice_HeatGeneration(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  // Determine which carrier equation is solved
  solveElectron = p.get<bool>("Solve Electron");
  solveHole = p.get<bool>("Solve Hole");
  solveIon = p.get<bool>("Solve Ion");

  // Determine if any R/G source is turned on
  haveSource = p.get<bool>("Have Source");

  // Evaluated fields
  heat_gen = MDField<ScalarT,Cell,Point>(n.field.heat_gen,scalar);
  this->addEvaluatedField(heat_gen);

  // Dependent fields
  if (solveElectron)
  {
    elec_curr_dens = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_curr_density,vector);
    elec_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield,vector);
    this->addDependentField(elec_curr_dens);
    this->addDependentField(elec_field);
  }
  if (solveHole)
  {
    hole_curr_dens = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_curr_density,vector);
    hole_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield,vector);
    this->addDependentField(hole_curr_dens);
    this->addDependentField(hole_field);
  }
  if (solveIon)
  {
    ion_curr_dens = MDField<const ScalarT,Cell,Point,Dim>(n.field.ion_curr_density,vector);
    ion_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.ion_efield,vector);
    this->addDependentField(ion_curr_dens);
    this->addDependentField(ion_field);
  }
  if (haveSource)
  {
    latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,scalar);
    total_recomb = MDField<const ScalarT,Cell,Point>(n.field.total_recomb,scalar);
    eff_band_gap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,scalar);
    this->addDependentField(latt_temp);
    this->addDependentField(total_recomb);
    this->addDependentField(eff_band_gap);
  }

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  J0 = scaleParams->scale_params.J0;
  E0 = scaleParams->scale_params.E0;
  H0 = scaleParams->scale_params.H0;
  T0 = scaleParams->scale_params.T0;

  std::string name = "DDLattice_HeatGeneration";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLattice_HeatGeneration<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;

  ScalarT scaling = J0*E0/H0;

  // obtain kb*T0 where T0 is the temperature scaling
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;       // Boltzmann constant in [eV/K]
  ScalarT kbT0 = kbBoltz*T0;  // in [eV]

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT JndotFn = 0.0;
      ScalarT JpdotFp = 0.0;
      ScalarT JidotFi = 0.0;
      ScalarT sourceHG = 0.0;

      if (solveElectron)  // include the Electron Joule heating
      {
        for (int dim = 0; dim < num_dims; ++dim)
          JndotFn += elec_curr_dens(cell,point,dim) * elec_field(cell,point,dim); // [scaled]
      }

      if (solveHole)  // include the Hole Joule heating
      {
        for (int dim = 0; dim < num_dims; ++dim)
          JpdotFp += hole_curr_dens(cell,point,dim) * hole_field(cell,point,dim); // [scaled]
      }

      if (solveIon)  // include the Ion Joule heating
      {
        for (int dim = 0; dim < num_dims; ++dim)
          JidotFi += ion_curr_dens(cell,point,dim) * ion_field(cell,point,dim);  // [scaled]
      }

      if (haveSource) // include the source heating
      {
        const ScalarT& recomb = total_recomb(cell,point); // [scaled]
        const ScalarT& temp = latt_temp(cell,point);  // [scaled]
        const ScalarT& bg = eff_band_gap(cell,point); // [eV]
        sourceHG = recomb*bg/kbT0 + 3.0*recomb*temp;  // [scaled]
      }

      ScalarT hg = JndotFn + JpdotFp + JidotFi + sourceHG;

      if (Sacado::ScalarValue<ScalarT>::eval(hg) > 0.0)
        heat_gen(cell,point) = hg * scaling; // [scaled]
      else
        heat_gen(cell,point) = 0.0;
    }
  }  // end of for loops
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DDLattice_HeatGeneration<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  p->set<bool>("Solve Ion", false);
  p->set<bool>("Solve Electron", false);
  p->set<bool>("Solve Hole", false);
  p->set<bool>("Have Source", false);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
