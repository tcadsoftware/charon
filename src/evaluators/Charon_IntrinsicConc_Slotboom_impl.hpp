
#ifndef CHARON_INTRINSICCONC_SLOTBOOM_IMPL_HPP
#define CHARON_INTRINSICCONC_SLOTBOOM_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
The BGN model is the new Slotboom model by D. B. M. Klaassen, J. W. Slotboom,
and H. C. De Graaff, "Unified apparent bandgap narrowing in n- and p-type silicon,"
Solid-State Electronics, Vol. 35, No. 2, pp. 125-129, 1992. 

To active this model, let Intrinsic Concentration = Slotboom in the input xml file.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IntrinsicConc_Slotboom<EvalT, Traits>::
IntrinsicConc_Slotboom(
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
  num_points = scalar->dimension(1);

  // Material name
  const string& materialName = p.get<string>("Material Name");

  // BGN flag
  const string bgn = p.get<string>("Band Gap Narrowing");
  if (bgn == "On")
    includeBGN = true;
  else
    includeBGN = false;

  // Obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve material parameters
  V0_BGN = matProperty.getPropertyValue(materialName, "Slotboom BGN V0_BGN");
  N0_BGN = matProperty.getPropertyValue(materialName, "Slotboom BGN N0_BGN");
  CON_BGN = matProperty.getPropertyValue(materialName, "Slotboom BGN CON_BGN");

  // Intrinsic Conc ParameterList
  if (p.isSublist("Intrinsic Conc ParameterList"))
  {
    const ParameterList& niParamList = p.sublist("Intrinsic Conc ParameterList");
    // Overwrite parameters when specified by users
    if (niParamList.isParameter("V0_BGN"))
      V0_BGN = niParamList.get<double>("V0_BGN");
    if (niParamList.isParameter("N0_BGN"))
      N0_BGN = niParamList.get<double>("N0_BGN");
    if (niParamList.isParameter("CON_BGN"))
      CON_BGN = niParamList.get<double>("CON_BGN");
  }

  // Evaluated fields
  nie = MDField<ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  effEg = MDField<ScalarT,Cell,Point>(n.field.eff_band_gap,scalar);
  effChi = MDField<ScalarT,Cell,Point>(n.field.eff_affinity,scalar);

  this->addEvaluatedField(nie);
  this->addEvaluatedField(effEg);
  this->addEvaluatedField(effChi);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  Eg = MDField<const ScalarT,Cell,Point>(n.field.band_gap,scalar);
  Chi = MDField<const ScalarT,Cell,Point>(n.field.affinity,scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);

  this->addDependentField(Eg);
  this->addDependentField(Chi);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(latt_temp);

  // When includeBGN = true, need acceptor and donor concentrations
  if (includeBGN)
  {
    acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,scalar);
    donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,scalar);

    this->addDependentField(acceptor);
    this->addDependentField(donor);
  }

  std::string name = "Intrinsic_Concentration_Slotboom";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IntrinsicConc_Slotboom<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // loop over the cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature [K]
      ScalarT lattT = latt_temp(cell,point)*T0;
      ScalarT kbT = kbBoltz*lattT;  // [eV]

      // obtain band gap without bgn
      const ScalarT& bgNobgn = Eg(cell,point);  // [eV], Eg in fm is NOT scaled.
      effEg(cell,point) = bgNobgn;

      // obtain electron affinity without bgn
      const ScalarT& chiNobgn = Chi(cell,point); // [eV]
      effChi(cell,point) = chiNobgn;

      ScalarT bgn_mod = 1.0;

      // calculate BGN modification when requested by user
      if (includeBGN)
      {
        const ScalarT& Na = acceptor(cell,point); // scaled
        const ScalarT& Nd = donor(cell,point);
        ScalarT Ntot = (Na + Nd) * C0;  // unscaled

        ScalarT deltaEg = 0.0;

        if (Ntot > 1e15)  // compute BGN for doping > 1e15 cm-3
        {
          ScalarT tmp1 = log(Ntot/N0_BGN);
          deltaEg = V0_BGN *(tmp1 + sqrt(tmp1*tmp1 + CON_BGN))*1.0; //[eV] and always > 0
        }

        effEg(cell,point) = bgNobgn - deltaEg;
        bgn_mod = exp(0.5*deltaEg/kbT);

        effChi(cell,point) = chiNobgn + 0.5*deltaEg;
      }

      // obtain the effective density of states (scaled)
      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);

      // compute the effective intrinsic density (scaled)
      nie(cell,point) = sqrt(Nc*Nv) * exp(-0.5*bgNobgn/kbT) * bgn_mod;

    } // end of point loop

  }  // end of cell loop

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
IntrinsicConc_Slotboom<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Band Gap Narrowing", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Intrinsic Conc ParameterList", false, "");
  p->sublist("Intrinsic Conc ParameterList").set<std::string>("Value", "Slotboom", "Use the Slotboom BGN model");
  p->sublist("Intrinsic Conc ParameterList").set<double>("V0_BGN", 0., "[V]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("N0_BGN", 0., "[cm^-3]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("CON_BGN", 0., "[1]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
