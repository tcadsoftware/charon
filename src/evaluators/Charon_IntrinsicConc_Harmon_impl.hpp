
#ifndef CHARON_INTRINSICCONC_HARMON_IMPL_HPP
#define CHARON_INTRINSICCONC_HARMON_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
The BGN model is the Harmon model by E. S. Harmon, M. R. Melloch, and M. S. Lundstrom,
"Effective band-gap shrinkage in GaAs," Appl. Phys. Lett. 64, 502 (1994).
To active this model, let Intrinsic Concentration = Harmon in the input xml file.

When BGN File is given, read in the doping-dependent dEc and dEv values and use them
as BGN values, and neglect the BGN due to the Harmon model.

By default, the Fermi-Dirac correction terms are not included. When
Enable Fermi Correction = true, then include the FD correction terms.

<ParameterList name="Intrinsic Concentration">
    <Parameter name="Value" type="string" value="Harmon" />
    <!-- <Parameter name="BGN File" type="string" value="BGN_File.txt" /> !-->
    <Parameter name="Enable Fermi Correction" type="bool" value="false" />
    <Parameter name="An" type="double" value="3.23e-8" />
    <Parameter name="Ap" type="double" value="2.55e-8" />
</ParameterList>

BGN_File.txt looks like the following:
dop1  dEc1  dEv1
dop2  dEc2  dEv2
...

*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IntrinsicConc_Harmon<EvalT, Traits>::
IntrinsicConc_Harmon(
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

  // BGN flag
  const string bgn = p.get<string>("Band Gap Narrowing");
  if (bgn == "On")
    includeBGN = true;
  else
    includeBGN = false;

  // Get the parameterlist
  const ParameterList& niParamList = p.sublist("Intrinsic Conc ParameterList");

  // Initialize the parameters
  initialize(niParamList);


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
    doping = MDField<const ScalarT,Cell,Point>(n.field.doping_raw,scalar);
    this->addDependentField(doping);
  }

  std::string name = "Intrinsic_Concentration_Harmon";
  this->setName(name);

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral =
    Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IntrinsicConc_Harmon<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // compute effEg, effChi, and nie without bgn
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // set effective band gap and electron affinity without bgn
      const ScalarT& bgNobgn = Eg(cell,point);
      effEg(cell,point) = bgNobgn;   // [eV]
      effChi(cell,point) = Chi(cell,point); // [eV]

      // obtain temperature [K]
      ScalarT lattT = latt_temp(cell,point)*T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

      ScalarT kbT = kbBoltz*lattT;  // [eV]

      // obtain the effective density of states (scaled)
      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);

      // compute the effective intrinsic density (scaled)
      nie(cell,point) = std::sqrt(Nc*Nv) * std::exp(-0.5*bgNobgn/kbT);
    }
  }

  // include bgn when requested by user and read bgn values from a file
  if (includeBGN && bgnFromFile)
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        // obtain temperature [K]
        ScalarT lattT = latt_temp(cell,point)*T0;

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

        ScalarT kbT = kbBoltz*lattT;  // [eV]

        // obtain net doping in [cm^-3]
        const ScalarT& dop = std::fabs(doping(cell,point)*C0);  // doping < 0 for p-type

        // obtain bgn values from a table read from a file
        ScalarT dEc = 0.0;
        ScalarT dEv = 0.0;
        evaluateBGNFromFile(dop, dEc, dEv);

        // include bgn in effEg, effChi, and nie
        ScalarT dEg = dEc + dEv;
        effEg(cell,point) -= dEg;
        effChi(cell,point) += dEc;
        nie(cell,point) *= std::exp(0.5*dEg/kbT);
      }
    }
  }

  // include bgn when requested by user and use the Harmon model when BGN File is not given
  if (includeBGN && (!bgnFromFile))
  {
    // include the doping-induced bgn
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        // obtain temperature [K]
        ScalarT lattT = latt_temp(cell,point)*T0;

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

        ScalarT kbT = kbBoltz*lattT;  // [eV]

        // obtain net doping in [cm^-3]
        const ScalarT& dop = doping(cell,point)*C0;
        ScalarT deltaEg = 0.0;

        if (dop > 0.0)  // n-type
          deltaEg = An * std::pow(dop, 1.0/3.0);  // [eV]
        else  // p-type
          deltaEg = Ap * std::pow(std::fabs(dop), 1.0/3.0);

        // include bgn in effEg, effChi, and nie
        effEg(cell,point) -= deltaEg;
        effChi(cell,point) += 0.5*deltaEg;
        nie(cell,point) *= std::exp(0.5*deltaEg/kbT);
      }
    }

    // include the Fermi-Dirac correction
    if (enableFD)
    {
      for (index_t cell = 0; cell < workset.num_cells; ++cell)
      {
        for (int point = 0; point < num_points; ++point)
        {
          // obtain temperature in [K]
          ScalarT lattT = latt_temp(cell,point)*T0;

          // lattT should be always > 0, but it could become <= 0 due to numerical errors
          // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
          if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

          ScalarT kbT = kbBoltz*lattT;  // [eV]

          // obtain net doping in [cm^-3]
          const ScalarT& dop = doping(cell,point)*C0;

          // obtain the effective density of states in [cm^-3]
          const ScalarT& Nc = elec_effdos(cell,point)*C0;
          const ScalarT& Nv = hole_effdos(cell,point)*C0;

          ScalarT ratio = 0.0;
          if (dop > 0.0)  // n-type
            ratio = dop / Nc;
          else  // p-type
            ratio = std::fabs(dop) / Nv;

          ScalarT fdCorrection = 0.0;
          if (ratio > 1e-4)  // FD correction is negligible for ratio <= 1e-4
            fdCorrection = kbT * (std::log(ratio) - (*inverseFermiIntegral)(ratio));

          // include FD correction in effEg, effChi, and nie
          effEg(cell,point) -= fdCorrection;
          effChi(cell,point) += 0.5*fdCorrection;
          nie(cell,point) *= std::exp(0.5*fdCorrection/kbT);
        }
      }
    }  // end if (enableFD)

  }  // end if (includeBGN)

}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void IntrinsicConc_Harmon<EvalT, Traits>::initialize(const Teuchos::ParameterList& plist)
{
  using std::string;
  using std::ifstream;

  // Set the default FD flag
  enableFD = false;

  // Set the default A coefficients
  An = 3.23e-8;  // in unit of [eV.cm], for n-GaAs
  Ap = 2.55e-8;  // in unit of [eV.cm], for p-GaAs

  // Overwrite parameters when specified by users
  if (plist.isParameter("Enable Fermi Correction"))
    enableFD = plist.get<bool>("Enable Fermi Correction");
  if (plist.isParameter("An"))
    An = plist.get<double>("An");
  if (plist.isParameter("Ap"))
    Ap = plist.get<double>("Ap");

  bgnFromFile = false;  // by default
  string bgnFile = "";  // Empty bgn file name by default
  if (plist.isParameter("BGN File"))
  {
    bgnFromFile = true;

    bgnFile = plist.get<string>("BGN File");
    if (bgnFile == "")
       TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "Error ! BGN File name cannot be empty !" << std::endl);

    ifstream inf(bgnFile.c_str());
    if (!inf)
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
            << "Error ! Cannot read the BGN File '" << bgnFile << "'" << std::endl);

    dopBGNStruct dbs;
    while (inf)
    {
      inf >> dbs.dop >> dbs.dEc >> dbs.dEv;
      dopBGNVec.push_back(dbs);
    }

    // sort the data according to doping values in ascending order
    std::sort(dopBGNVec.begin(), dopBGNVec.end());

    // eliminate duplicate (dop,dEc,dEv) entry
    dopBGNVec.resize( unique (dopBGNVec.begin(), dopBGNVec.end()) - dopBGNVec.begin() );

    // BGN data after sorting and eliminating duplicate entry
    // std::cout << "BGN data after sorting and eliminating duplicate entry ..." << std::endl;
    // for (int i = 0; i < dopBGNVec.size(); i++)
    //   std::cout << dopBGNVec[i].dop << ", " << dopBGNVec[i].dEc << ", " << dopBGNVec[i].dEv << std::endl;

  }  // end of if (BGN File)
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateBGNFromFile()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void IntrinsicConc_Harmon<EvalT, Traits>::evaluateBGNFromFile(const ScalarT& dop, ScalarT& dEc, ScalarT& dEv)
{
  int index = binarySearch(dop);
  int vectorSize = dopBGNVec.size();

  // index is located at the last ending point
  if (index >= vectorSize-1)
  {
    dEc = dopBGNVec[index].dEc;
    dEv = dopBGNVec[index].dEv;
  }
  else  // use linear interpolation
  {
    // obtain dEc
    ScalarT slope = (dopBGNVec[index+1].dEc - dopBGNVec[index].dEc)
                  / (dopBGNVec[index+1].dop - dopBGNVec[index].dop);
    dEc = slope * (dop - dopBGNVec[index].dop) + dopBGNVec[index].dEc;

    // obtain dEv
    slope = (dopBGNVec[index+1].dEv - dopBGNVec[index].dEv)
          / (dopBGNVec[index+1].dop - dopBGNVec[index].dop);
    dEv = slope * (dop - dopBGNVec[index].dop) + dopBGNVec[index].dEv;
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  binarySearch()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
int IntrinsicConc_Harmon<EvalT, Traits>::binarySearch(const ScalarT& dop)
{
  int vectorSize = dopBGNVec.size();
  int xmin, xmax;
  bool forward = true;

  if (dopBGNVec[0].dop < dopBGNVec[vectorSize-1].dop)  // ascendingly sorted
  {
    xmin = 0;
    xmax = vectorSize-1;
  }
  else  // descendingly sorted
  {
    xmin = vectorSize-1;
    xmax = 0;
    forward=false;
  }

  //if dop is not in the vector, return the closest endpoint
  if (dop <= dopBGNVec[xmin].dop)  return xmin;
  if (dop >= dopBGNVec[xmax].dop)  return xmax;

  int imin = 0;
  int imax = vectorSize-1;

  while (imax > imin + 1)
  {
    int imid = (imin+imax)/2;  // equivalent to imin+(imax-imin)/2 ? Note the integer division
    if (forward)
    {
      if (dop < dopBGNVec[imid].dop)
        imax = imid;
      else
        imin = imid;
    }
    else
    {
      if (dop > dopBGNVec[imid].dop)
        imax = imid;
      else
        imin = imid;
    }
  }

  return imin;
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
IntrinsicConc_Harmon<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Band Gap Narrowing", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Intrinsic Conc ParameterList", false, "");
  p->sublist("Intrinsic Conc ParameterList").set<std::string>("Value", "Harmon", "Use the Harmon BGN model");
  p->sublist("Intrinsic Conc ParameterList").set<double>("An", 0., "The A coefficient for n-GaAs in [eV.cm]");
  p->sublist("Intrinsic Conc ParameterList").set<double>("Ap", 0., "The A coefficient for p-GaAs in [eV.cm]");
  p->sublist("Intrinsic Conc ParameterList").set<bool>("Enable Fermi Correction", false, "Include the FD correction in the bgn calculation");
  p->sublist("Intrinsic Conc ParameterList").set<std::string>("BGN File", "", "Get doping-dependent BGN values from an external file");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
