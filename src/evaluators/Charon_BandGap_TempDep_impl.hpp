
#ifndef CHARON_BANDGAP_TEMPDEP_IMPL_HPP
#define CHARON_BANDGAP_TEMPDEP_IMPL_HPP

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
BandGap_TempDep<EvalT, Traits>::
BandGap_TempDep(
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

  // retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // material name
  const string& materialName = p.get<string>("Material Name");

  // does compute affinity ?
  isCalcAffinity = p.get<bool>("Compute Affinity");

  // bandgap parameterlist
  const ParameterList& bgParamList = p.sublist("Bandgap ParameterList");

  // obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // retrieve bandgap temp-dep model parameters
  Eg300 = matProperty.getPropertyValue(materialName, "Band Gap at 300 K");
  alpha = matProperty.getPropertyValue(materialName, "Band Gap alpha");
  beta = matProperty.getPropertyValue(materialName, "Band Gap beta");
  Chi300 = matProperty.getPropertyValue(materialName, "Electron Affinity at 300 K");

  // overwrite the parameters when specified by users
  if (bgParamList.isParameter("Eg300"))
    Eg300 = bgParamList.get<double>("Eg300");
  if (bgParamList.isParameter("alpha"))
    alpha = bgParamList.get<double>("alpha");
  if (bgParamList.isParameter("beta"))
    beta = bgParamList.get<double>("beta");
  if (bgParamList.isParameter("Chi300"))
    Chi300 = bgParamList.get<double>("Chi300");

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // fields
  band_gap = MDField<ScalarT,Cell,Point>(n.field.band_gap,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);

  // evaluated fields
  this->addEvaluatedField(band_gap);

  // dependent fields
  this->addDependentField(latt_temp);

  // need to compute electron affinity
  if (isCalcAffinity)
  {
    affinity = MDField<ScalarT,Cell,Point>(n.field.affinity,scalar);
    this->addEvaluatedField(affinity);
  }

  std::string name = "BandGap_Temperature_Dependence";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BandGap_TempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Note all energy-related fields saved in the FM are in the units of eV,
  // and all other fields are scaled !

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature in [K]
      ScalarT lattT = latt_temp(cell,point) * T0;

      // compute band gap in [eV]
      band_gap(cell,point) = Eg300 + alpha * 300.*300. / (300.+beta)
                                   - alpha * lattT*lattT / (lattT+beta) ;

      // compute electron affinity in [eV] when not given by user
      if (isCalcAffinity)
        affinity(cell,point) = Chi300 - alpha *300.*300. / (600.+2.*beta)
                             + alpha *lattT*lattT / (2.*lattT + 2.*beta);
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
BandGap_TempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");
  p->set<bool>("Compute Affinity", false);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Bandgap ParameterList", false, "");
  p->sublist("Bandgap ParameterList").set<std::string>("Value", "TempDep", "Specify temperature-dep band gap");
  p->sublist("Bandgap ParameterList").set<double>("Eg300", 0., "Band gap at 300 K in [eV]");
  p->sublist("Bandgap ParameterList").set<double>("Chi300", 0., "Electron Affinity at 300 K in [eV]");
  p->sublist("Bandgap ParameterList").set<double>("alpha", 0., "alpha coeff [eV/K] in calculating Eg");
  p->sublist("Bandgap ParameterList").set<double>("beta", 0., "beta coeff [K] in calculating Eg");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
