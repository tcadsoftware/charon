
#ifndef CHARON_IC_EQUILIBRIUM_DENSITY_IMPL_HPP
#define CHARON_IC_EQUILIBRIUM_DENSITY_IMPL_HPP

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
IC_Equilibrium_Density<EvalT, Traits>::
IC_Equilibrium_Density(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  m_names = p.get< Teuchos::RCP< const charon::Names> >("Names");

  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> scalar = basis->functional;
  num_basis = scalar->dimension(1);

  dof_name = p.get<string>("DOF Name");
  string edens_name = "ELECTRON_DENSITY";
  string hdens_name = "HOLE_DENSITY";

  if ( (dof_name.length() == edens_name.length()) ||
       (dof_name.length() == hdens_name.length()) )
    haveSuffix = false;
  else
    haveSuffix = true;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters(haveSuffix);
  p.validateParameters(*valid_params);

  // Obtain the parameterlist
  const ParameterList& plist = p.sublist("Equilibrium ParameterList");

  if ((plist.isParameter("Effective Band Gap")) && !(plist.isParameter("Effective Electron Affinity")))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Effective Band Gap and Effective Electron Affinity must be specified together !");

  if (!(plist.isParameter("Effective Band Gap")) && (plist.isParameter("Effective Electron Affinity")))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Effective Band Gap and Effective Electron Affinity must be specified together !");

  if ((plist.isParameter("Effective Band Gap")) && (plist.isParameter("Effective Electron Affinity")))
  {
    effBandGap = plist.get<double>("Effective Band Gap");
    effAffinity = plist.get<double>("Effective Electron Affinity");
  }
  else
  {
    effBandGap = -1.0;
    effAffinity = -1.0;
  }

  // Evaluated field
  carrier_density = MDField<ScalarT,Cell,BASIS>(dof_name, scalar);
  this->addEvaluatedField(carrier_density);

  // Scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  elec_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.elec_eff_dos, scalar);
  hole_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.hole_eff_dos, scalar);
  cond_band = MDField<const ScalarT,Cell,BASIS>(m_names->field.cond_band, scalar);
  vale_band = MDField<const ScalarT,Cell,BASIS>(m_names->field.vale_band, scalar);
  latt_temp = MDField<const ScalarT,Cell,BASIS>(m_names->field.latt_temp, scalar);

  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(cond_band);
  this->addDependentField(vale_band);
  this->addDependentField(latt_temp);

  if (haveSuffix && (effBandGap > 0.0))
  {
    ref_energy = MDField<const ScalarT,Cell,BASIS>(m_names->field.ref_energy, scalar);
    potential = MDField<const ScalarT,Cell,BASIS>(m_names->dof.phi, scalar);
    this->addDependentField(ref_energy);
    this->addDependentField(potential);
  }

  std::string name = "IC_Equilibrium_Density";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Equilibrium_Density<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  // obtain physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;     // Boltzmann constant in [eV/K]
  ScalarT kbT0 = kb * T0;  // [eV[

  for (int cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      ScalarT kbT = latt_temp(cell,basis) * T0 * kb;  // [eV]
      const ScalarT& Nc = elec_effdos(cell, basis);  // scaled
      const ScalarT& Nv = hole_effdos(cell, basis);

      if (dof_name == m_names->dof.edensity)
      {
        ScalarT Ec = 0.0;
        if (!haveSuffix)  // for the case of continuous carrier densities
          Ec = cond_band(cell, basis);
        else   // for the case of discontinuous carrier densities at a heterojunction
        {
          if (effBandGap > 0.0)  // compute Ec
            Ec = ref_energy(0,0) - effAffinity - potential(cell,basis) * kbT0;
          else
            Ec = cond_band(cell, basis);
        }
        carrier_density(cell,basis) = Nc * std::exp(-Ec/kbT); // n0 = Nc*exp(-Ec/kbT), since Ef=0 at equilibrium

        // std::cout <<"Nc=" << Nc << ", Ec=" << Ec << ", dof_name=" << dof_name << ", carr_dens=" << carrier_density(cell,basis) << std::endl;
      }

      else if (dof_name == m_names->dof.hdensity)
      {
        ScalarT Ev = 0.0;
        if (!haveSuffix)  // for the case of continuous carrier densities
          Ev = vale_band(cell, basis);
        else   // for the case of discontinuous carrier densities at a heterojunction
        {
          if (effBandGap > 0.0)  // compute Ev
            Ev = ref_energy(0,0) - effAffinity - potential(cell,basis)*kbT0 - effBandGap;
          else
            Ev = vale_band(cell, basis);
        }
        carrier_density(cell,basis) = Nv * std::exp(Ev/kbT); // p0 = Nv*exp(Ev/kbT), since Ef=0 at equilibrium

        // std::cout <<"Nv=" << Nv << ", Ev=" << Ev << ", dof_name=" << dof_name << ", carr_dens=" << carrier_density(cell,basis) << std::endl;
      }
    }  // loop over basis
  }  // loop over cell
}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
IC_Equilibrium_Density<EvalT, Traits>::getValidParameters(bool haveSuffix) const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("DOF Name", "?");

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->sublist("Equilibrium ParameterList", false, "");
  p->sublist("Equilibrium ParameterList").set<std::string>("Value", "Equilibrium Density", "Compute equilibrium density for initial condition");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  if (haveSuffix)
  {
    p->sublist("Equilibrium ParameterList").set<double>("Effective Band Gap", 0.0, "[eV]");
    p->sublist("Equilibrium ParameterList").set<double>("Effective Electron Affinity", 0.0, "[eV]");
  }

  return p;
}

}

#endif

