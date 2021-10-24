
#ifndef CHARON_REFERENCE_ENERGY_IMPL_HPP
#define CHARON_REFERENCE_ENERGY_IMPL_HPP

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
Reference_Energy<EvalT, Traits>::
Reference_Energy(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  numPoints = scalar->dimension(1);

  // Reference material name
  const string& refMaterial = p.get<string>("Reference Material");

  // Obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve bandgap temp-dep model parameters
  Eg300 = matProperty.getPropertyValue(refMaterial, "Band Gap at 300 K");
  Chi300 = matProperty.getPropertyValue(refMaterial, "Electron Affinity at 300 K");
  alpha = matProperty.getPropertyValue(refMaterial, "Band Gap alpha");
  beta = matProperty.getPropertyValue(refMaterial, "Band Gap beta");

  // Retrieve other parameters
  Nc300 = matProperty.getPropertyValue(refMaterial, "Electron Effective DOS at 300 K");
  Nv300 = matProperty.getPropertyValue(refMaterial, "Hole Effective DOS at 300 K");
  Nc_F = matProperty.getPropertyValue(refMaterial, "Electron Effective DOS Exponent");
  Nv_F = matProperty.getPropertyValue(refMaterial, "Hole Effective DOS Exponent");

  // Default values
  constBg = 0.0;
  constChi = 0.0;
  isBgConst = false;
  isChiConst = false;

  // Constant electron affinity is specified ?
  if (p.isParameter("Constant Electron Affinity"))
  {
    isChiConst = true;
    constChi = p.get<double>("Constant Electron Affinity");
  }

  // Constant band gap is specified ?
  if (p.isParameter("Constant Band Gap"))
  {
    isBgConst = true;
    constBg = p.get<double>("Constant Band Gap");
  }

  // Bandgap parameterlist is specified ?
  if (p.isSublist("Bandgap ParameterList"))
  {
    const ParameterList& bgParamList = p.sublist("Bandgap ParameterList");

    // TempDep band gap is specified ?
    if ( (bgParamList.isType<string>("Value")) &&
         (bgParamList.get<string>("Value") == "TempDep") )
    {
      // Overwrite the parameters when given by users
      if (bgParamList.isParameter("Eg300"))
        Eg300 = bgParamList.get<double>("Eg300");
      if (bgParamList.isParameter("Chi300"))
        Chi300 = bgParamList.get<double>("Chi300");
      if (bgParamList.isParameter("alpha"))
        alpha = bgParamList.get<double>("alpha");
      if (bgParamList.isParameter("beta"))
        beta = bgParamList.get<double>("beta");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid Band Gap model: "
        << bgParamList.get<string>("Value") << "!" << std::endl);
  }

  // Effective DOS parameterlist is given
  if (p.isSublist("Effective DOS ParameterList"))
  {
    const ParameterList& effDosParamList = p.sublist("Effective DOS ParameterList");

    // Effective DOS = Simple is specified
    if ( (effDosParamList.isType<string>("Value")) &&
         (effDosParamList.get<string>("Value") == "Simple") )
    {
      // Overwrite the parameters when given by users
      if (effDosParamList.isParameter("Nc300"))
        Nc300 = effDosParamList.get<double>("Nc300");
      if (effDosParamList.isParameter("Nv300"))
        Nv300 = effDosParamList.get<double>("Nv300");
      if (effDosParamList.isParameter("Nc_F"))
        Nc_F = effDosParamList.get<double>("Nc_F");
      if (effDosParamList.isParameter("Nv_F"))
        Nv_F = effDosParamList.get<double>("Nv_F");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid Effective DOS model: "
        << effDosParamList.get<string>("Value") << "!" << std::endl);
  }

  // Fields
  ref_energy = MDField<ScalarT,Cell,Point>(n.field.ref_energy,scalar);

  // Scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Evaluated fields
  this->addEvaluatedField(ref_energy);

  std::string name = "Constant_Reference_Energy";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Reference_Energy<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Note the Reference Energy saved in the FM has the unit of eV !

  double lattT = T0; // in unit of [K]

  // obtain kb*T
  const charon::PhysicalConstants & cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;          // Boltzmann constant in [eV/K]
  double kbT = kbBoltz*lattT;         // [eV]

  // calculate the effective density of states
  ScalarT Nc = Nc300 * pow(lattT/300.0, Nc_F);
  ScalarT Nv = Nv300 * pow(lattT/300.0, Nv_F);

  // calculate the electron affinity 
  ScalarT affinity = constChi;  // const affinity is given
  if (!isChiConst)              // compute the affinity when it is not given
    affinity = Chi300 - alpha*300.*300. / (2.*(300.+beta))
                      + alpha*lattT*lattT / (2.*(lattT+beta));

  // calculate the band gap
  ScalarT bandgap = constBg;
  if (!isBgConst)
    bandgap = Eg300 + alpha*300.*300./(300.+beta) - alpha*lattT*lattT/(lattT+beta);

  // compute the reference energy in [eV], which is the intrinsic Fermi energy
  // level of the Reference Material from the vacuum level
  ScalarT Eref = affinity + 0.5*bandgap + 0.5*kbT*log(Nc/Nv);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    for (int point = 0; point < numPoints; ++point)
      ref_energy(cell,point) = Eref;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Reference_Energy<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Reference Material", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<double>("Constant Band Gap", 0., "[eV]");
  p->set<double>("Constant Electron Affinity", 0., "[eV]");

  p->sublist("Bandgap ParameterList", false, "");
  p->sublist("Bandgap ParameterList").set<std::string>("Value", "TempDep", "Specify temperature-dep band gap");
  p->sublist("Bandgap ParameterList").set<double>("Eg300", 0., "Band gap at 300 K in [eV]");
  p->sublist("Bandgap ParameterList").set<double>("Chi300", 0., "Electron Affinity at 300 K in [eV]");
  p->sublist("Bandgap ParameterList").set<double>("alpha", 0., "alpha coeff [eV/K] in calculating Eg");
  p->sublist("Bandgap ParameterList").set<double>("beta", 0., "beta coeff [K] in calculating Eg");

  p->sublist("Effective DOS ParameterList", false, "");
  p->sublist("Effective DOS ParameterList").set<std::string>("Value", "Simple", "");
  p->sublist("Effective DOS ParameterList").set<double>("Nc300", 0., "[cm^-3]");
  p->sublist("Effective DOS ParameterList").set<double>("Nv300", 0., "[cm^-3]");
  p->sublist("Effective DOS ParameterList").set<double>("Nc_F", 0., "[1]");
  p->sublist("Effective DOS ParameterList").set<double>("Nv_F", 0., "[1]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
