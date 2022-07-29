
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
//  using Teuchos::ParameterList;

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
  doping = MDField<const ScalarT,Cell,Point>(n->field.doping_raw,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n->field.latt_temp,scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n->field.elec_eff_dos,scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n->field.hole_eff_dos,scalar);
  condband = MDField<const ScalarT,Cell,Point>(n->field.cond_band,scalar);
  valeband = MDField<const ScalarT,Cell,Point>(n->field.vale_band,scalar);

  // Determine if Maxwell-Boltzmann or Fermi-Dirac statistics
  UseFD = p.get<string>("Fermi Dirac"); 

  // Add fields
  this->addEvaluatedField(nlpsrc);
  this->addDependentField(doping);
  this->addDependentField(latt_temp);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(condband);
  this->addDependentField(valeband);

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

      const ScalarT& dop = doping(cell,point);       // scaled

      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);
      const ScalarT& Condband = condband(cell,point);  // [eV]
      const ScalarT& Valeband = valeband(cell,point);  // [eV]

      // nlpsrc(cell,point) = nie*exp(-phi) - nie*exp(phi) + dop;
      // nlpsrc(cell,point) = nie*exp(Ei/kbT) - nie*exp(-Ei/kbT) + dop;

      // X. Gao et al., J. Appl. Phys. 114, 164302 (2013)   
      if(UseFD == "False")
        nlpsrc(cell,point) = Nv*exp(Valeband/kbT) - Nc*exp(-Condband/kbT) + dop;
      else
        nlpsrc(cell,point) = Nv*Fhalf(Valeband/kbT) - Nc*Fhalf(-Condband/kbT) + dop;

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

  p->set<std::string>("Fermi Dirac","False");
  
  return p;
}

///////////////////////////////////////////////////////////////////////////////
//
// Fermi-Dirac integral of 1/2 order
//
///////////////////////////////////////////////////////////////////////////////

// X. Gao et al., J. Appl. Phys. 114, 164302 (2013)

template<typename EvalT, typename Traits>
typename NLPoissonSource<EvalT,Traits>::ScalarT
NLPoissonSource<EvalT, Traits>::Fhalf( const ScalarT& x)
{
  if (x > -50.0)
  {
    ScalarT v = x*x*x*x + 50.0 + 33.6*x*(1.0-0.68*exp( -0.17*(x+1.0)*(x+1.0) ) );
    return 1.0/(exp(-x) + 1.329340388*pow(v,-0.375) );
  }
  else
    return exp(x);  
}


}

#endif

