
#ifndef CHARON_INTRINSICCONC_PERSSON_IMPL_HPP
#define CHARON_INTRINSICCONC_PERSSON_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
The model is by C. Persson, U. Lindefelt, and
B. E. Sernelius, "Band gap narrowing in n-type and p-type 3C-, 2H-, 4H-, 6H-SiC,
and Si," Journal of Applied Phyiscs 86, 4419 (1999). The default parameters in
charon::Material_Properties are from Table II and III in the reference.
Charon 1.0 also uses the parameters by Persson et al.

This Persson model is similar to the Jain-Roulston model by S. C. Jain and D. J. Roulston,
"A simple expression for band gap narrowing (BGN) in heavily doped Si, Ge, GaAs
and Ge(x)Si(1-x) strained layers," Solid-State Electronics Vol. 34, No. 5,
pp. 453-465, 1991. Note the parameters are somewhat different from the Persson's
paper. 

The model uses general closed-form equations for band gap narrowing for n-type
and p-type semiconductors. The equations are dervied by identifying the exchange
energy shift of the majority band edge, correlation energy shift of the minority
band edge, shift of the majority band edge due to carrier-impurity interactions,
and shift of the minority band edge due to carrier-impurity interactions.

Jain's paper assumes parabolic band, while Persson's work goes beyond the common
parabolic treatments of the ground state energy dispersion by including energy
dispersion and overlap integrals from band structure calculations. The nonparabolic
valence band curvatures influence strongly the energy shifts especially in
p-type materials.

*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IntrinsicConc_Persson<EvalT, Traits>::
IntrinsicConc_Persson(
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
  anc = matProperty.getPropertyValue(materialName, "Persson ANC_BGN");
  bnc = matProperty.getPropertyValue(materialName, "Persson BNC_BGN");
  cnc = matProperty.getPropertyValue(materialName, "Persson CNC_BGN");
  anv = matProperty.getPropertyValue(materialName, "Persson ANV_BGN");
  bnv = matProperty.getPropertyValue(materialName, "Persson BNV_BGN");
  cnv = matProperty.getPropertyValue(materialName, "Persson CNV_BGN");

  apc = matProperty.getPropertyValue(materialName, "Persson APC_BGN");
  bpc = matProperty.getPropertyValue(materialName, "Persson BPC_BGN");
  cpc = matProperty.getPropertyValue(materialName, "Persson CPC_BGN");
  apv = matProperty.getPropertyValue(materialName, "Persson APV_BGN");
  bpv = matProperty.getPropertyValue(materialName, "Persson BPV_BGN");
  cpv = matProperty.getPropertyValue(materialName, "Persson CPV_BGN");

  // Intrinsic Conc ParameterList
  if (p.isSublist("Intrinsic Conc ParameterList"))
  {
    const ParameterList& niParamList = p.sublist("Intrinsic Conc ParameterList");

    // Overwrite parameters when specified by users
    if (niParamList.isParameter("ANC_BGN"))
      anc = niParamList.get<double>("ANC_BGN");
    if (niParamList.isParameter("BNC_BGN"))
      bnc = niParamList.get<double>("BNC_BGN");
    if (niParamList.isParameter("CNC_BGN"))
      cnc = niParamList.get<double>("CNC_BGN");
    if (niParamList.isParameter("ANV_BGN"))
      anv = niParamList.get<double>("ANV_BGN");
    if (niParamList.isParameter("BNV_BGN"))
      bnv = niParamList.get<double>("BNV_BGN");
    if (niParamList.isParameter("CNV_BGN"))
      cnv = niParamList.get<double>("CNV_BGN");

    if (niParamList.isParameter("APC_BGN"))
      apc = niParamList.get<double>("APC_BGN");
    if (niParamList.isParameter("BPC_BGN"))
      bpc = niParamList.get<double>("BPC_BGN");
    if (niParamList.isParameter("CPC_BGN"))
      cpc = niParamList.get<double>("CPC_BGN");
    if (niParamList.isParameter("APV_BGN"))
      apv = niParamList.get<double>("APV_BGN");
    if (niParamList.isParameter("BPV_BGN"))
      bpv = niParamList.get<double>("BPV_BGN");
    if (niParamList.isParameter("CPV_BGN"))
      cpv = niParamList.get<double>("CPV_BGN");
  }

  // Evaluated fields
  nie = MDField<ScalarT,Cell,Point>(n.field.intrin_conc, scalar);
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
  Eg = MDField<const ScalarT,Cell,Point>(n.field.band_gap, scalar);
  Chi = MDField<const ScalarT,Cell,Point>(n.field.affinity, scalar);
  elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos, scalar);
  hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos, scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);

  this->addDependentField(Eg);
  this->addDependentField(Chi);
  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(latt_temp);

  // When includeBGN = true, need acceptor and donor concentrations
  if (includeBGN)
  {
    acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw, scalar);
    donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw, scalar);

    this->addDependentField(acceptor);
    this->addDependentField(donor);
  }

  std::string name = "Intrinsic_Concentration_Persson";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IntrinsicConc_Persson<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& phyConst = charon::PhysicalConstants::Instance();
  double kbBoltz = phyConst.kb;       // Boltzmann constant in [eV/K]

  // loop over the cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature [K]
      ScalarT lattT = latt_temp(cell,point)*T0;
      ScalarT kbT = kbBoltz*lattT;  // [eV]

      // obtain band gap without bgn
      const ScalarT& bgNobgn = Eg(cell,point); //[eV],Eg in the fm is NOT scaled.
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
        ScalarT Nnet = (Nd - Na) * C0;         // unscaled
        ScalarT tmp1 = std::abs(Nnet/1.0e18);

        ScalarT deltaEc = 0.0;
        ScalarT deltaEv = 0.0;

        // use the net doping
        if (Nnet > 0.0)   // n-type
        {
          deltaEc = anc*pow(tmp1, 1.0/3.0) + bnc*pow(tmp1, 0.25) + cnc*pow(tmp1, 0.5);
          deltaEv = anv*pow(tmp1, 1.0/3.0) + bnv*pow(tmp1, 0.25) + cnv*pow(tmp1, 0.5);
        }
        else   // p-type
        {
          deltaEc = apc*pow(tmp1, 1.0/3.0) + bpc*pow(tmp1, 0.25) + cpc*pow(tmp1, 0.5);
          deltaEv = apv*pow(tmp1, 1.0/3.0) + bpv*pow(tmp1, 0.25) + cpv*pow(tmp1, 0.5);
        }

        ScalarT deltaEg = deltaEv - deltaEc;  // [eV] and always > 0
        effEg(cell,point) = bgNobgn - deltaEg;
        bgn_mod = exp(0.5*deltaEg/kbT);

        effChi(cell,point) = chiNobgn - deltaEc;   // [eV], deltaEc < 0 given in Persson's paper
      }

      // obtain the effective density of states (scaled)
      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);

      // compute the intrinsic density (scaled)
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
IntrinsicConc_Persson<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Material Name", "?");
  p->set<std::string>("Band Gap Narrowing", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Intrinsic Conc ParameterList", false, "");
  p->sublist("Intrinsic Conc ParameterList").set<std::string>("Value", "Persson", "Use the Persson model");

  p->sublist("Intrinsic Conc ParameterList").set<double>("ANC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("BNC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("CNC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("ANV_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("BNV_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("CNV_BGN", 0., "[eV]");

  p->sublist("Intrinsic Conc ParameterList").set<double>("APC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("BPC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("CPC_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("APV_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("BPV_BGN", 0., "[eV]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("CPV_BGN", 0., "[eV]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
