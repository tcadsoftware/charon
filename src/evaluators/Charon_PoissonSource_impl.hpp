
#ifndef CHARON_POISSONSOURCE_IMPL_HPP
#define CHARON_POISSONSOURCE_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
PoissonSource<EvalT, Traits>::
PoissonSource(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // Retrieve the data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Flags
  solveElectron = p.get<string>("Solve Electron");
  solveHole = p.get<string>("Solve Hole");

  // Additional parameters for DDIon and DDIonLattice equation sets
  solveIon = false;  // default
  ionCharge = 0;
  if (p.isParameter("Solve Ion"))
  {
    solveIon = p.get<bool>("Solve Ion");
    ionCharge = p.get<int>("Ion Charge");
  }

  // Carrier-independent fields
  poissonSource = MDField<ScalarT,Cell,Point>(p.get<string>("Source Name"),scalar);
  //intrin_fermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi,scalar);
  doping = MDField<const ScalarT,Cell,Point>(n.field.doping,scalar);
  //intrin_conc = MDField<const ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);

  this->addEvaluatedField(poissonSource);
  //this->addDependentField(intrin_fermi);
  this->addDependentField(doping);
  //this->addDependentField(intrin_conc);
  this->addDependentField(latt_temp);

  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Carrier-dependent fields
  if (solveElectron == "True")
  {
    edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
    this->addDependentField(edensity);
  }
  if (solveHole == "True")
  {
    hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
    this->addDependentField(hdensity);
  }

  // Add the ion density dependence
  if (solveIon)
  {
    iondensity = MDField<const ScalarT,Cell,Point>(n.dof.iondensity,scalar);
    this->addDependentField(iondensity);
  }

  std::string name = "Poisson Source";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
PoissonSource<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain the Boltzmann constant
  // charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  // double kb = cpc.kb;      // Boltzmann constant in [eV/K]

  // Loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // const ScalarT& Ei = intrin_fermi(cell, point);  // [eV]

      const ScalarT& dop = doping(cell, point);      // scaled
      // const ScalarT& nie = intrin_conc(cell, point); // scaled

      // ScalarT lattT = latt_temp(cell, point)*T0; // [K]

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      // if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      ScalarT eden, hden;
      // ScalarT kbT = kb*lattT;   // [eV]

      if (solveElectron == "True")
      {
        eden = edensity(cell, point);
        if (eden < 0.0)  eden = 0.0;
      }
      else
        // eden = nie*exp(-Ei/kbT);  // assume electron quasi-Fermi remains at 0.
        eden = 0.0;

      if (solveHole == "True")
      {
        hden = hdensity(cell, point);
        if (hden < 0.0) hden = 0.0;
      }
      else
        // hden = nie*exp(Ei/kbT); // assume hole quasi-Fermi remains at 0.
        hden = 0.0;

      // add the ion contribution when solveIon = true
      if (solveIon)
      {
        const ScalarT& ionden = iondensity(cell, point);  // scaled
        poissonSource(cell, point) = hden - eden + dop + ionCharge * ionden;
      }

      // otherwise, just (p-n+dop)
      else
        poissonSource(cell, point) = hden - eden + dop;

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
PoissonSource<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Source Name", "?");

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set("Solve Electron", "?");
  p->set("Solve Hole", "?");

  p->set<bool>("Solve Ion", false, "By default, do not solve the ion continuity equation");
  p->set<int>("Ion Charge", 0, "By default, ion charge = 0");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

