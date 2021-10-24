
#ifndef CHARON_THERMODIFFCOEFF_CUSTOM_IMPL_HPP
#define CHARON_THERMODIFFCOEFF_CUSTOM_IMPL_HPP

#include <cmath>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

/*
The model implements the ion thermodiffusion coefficient as
D_{T} = D_{T0}/(kb*T^2)*exp(-Ea/kbT)*sign, and D_{Ts} = D_{T}*T0/D0, where T0 and D0
are scaling factors, and
Ea = Emax for T < Tmin,
Ea = Emin for T > Tmax,
Ea = slope*(T-Tmin) + Emax for Tmin <= T <= Tmax, with
slope = (Emax - Emin) / (Tmin - Tmax) in unit of [eV/K].

This model is applicable to ions/vacancies only, neither electrons nor holes.

*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
ThermodiffCoeff_Custom<EvalT, Traits>::
ThermodiffCoeff_Custom(
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
  num_edges = num_points;  // initialization

  // input from closure model factory
  isEdgedl = p.get<bool>("Is Edge Data Layout");

  // input from ParameterList
  const ParameterList& plist = p.sublist("Thermodiffusion Coefficient ParameterList");

  // Initialize the parameters
  initialize(plist);

  // retrieve edge data layout if isEdgedl = true
  RCP<DataLayout> output_dl = scalar;
  if (isEdgedl)
  {
    RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
    RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
    output_dl = cellTopoInfo->edge_scalar;
    num_edges = output_dl->dimension(1);
    cellType = cellTopoInfo->getCellTopology();
  }

  // evaluated field
  thermodiff_coeff = MDField<ScalarT,Cell,Point>(n.field.ion_thermodiff_coeff,output_dl);
  this->addEvaluatedField(thermodiff_coeff);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  D0 = scaleParams->scale_params.D0;

  // dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);

  this->addDependentField(latt_temp);

  std::string name = "Thermodiffusion_Coefficient_Custom";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
ThermodiffCoeff_Custom<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;       // Boltzmann constant in [eV/K]

  if (isEdgedl) // evaluate thermodiff. coeff. at edge midpoints
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        const ScalarT& latt = ((latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0) * T0;   // [K]

        ScalarT tdc = sign * multiplier / (kb * latt * latt);  // [cm^2/(s.K)]

        ScalarT actE = 0.0;  // [eV]
        if (Sacado::ScalarValue<ScalarT>::eval(latt) < minTemp)
          actE = maxActE;
        else if (Sacado::ScalarValue<ScalarT>::eval(latt) > maxTemp)
          actE = minActE;
        else
          actE = slope * (latt - minTemp) + maxActE;

        tdc *= std::exp(-actE / (kb*latt) );

        thermodiff_coeff(cell,edge) = tdc * T0 / D0;
      }
    }
  }

  else  // evaluate thermodiff. coeff. at IP or BASIS points
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        const ScalarT& latt = latt_temp(cell,point) * T0;   // [K]

        ScalarT tdc = sign * multiplier / (kb * latt * latt);  // [cm^2/(s.K)]

        ScalarT actE = 0.0;  // [eV]
        if (Sacado::ScalarValue<ScalarT>::eval(latt) < minTemp)
          actE = maxActE;
        else if (Sacado::ScalarValue<ScalarT>::eval(latt) > maxTemp)
          actE = minActE;
        else
          actE = slope * (latt - minTemp) + maxActE;

        tdc *= std::exp(-actE / (kb*latt) );

        thermodiff_coeff(cell,point) = tdc * T0 / D0;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  initialize()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void ThermodiffCoeff_Custom<EvalT, Traits>::initialize
(const Teuchos::ParameterList& plist)
{
  using std::string;

  // Set up parameters for the model
  multiplier = plist.get<double>("Thermodiffusion Multiplier");

  // Thermodiffusion Coefficient Sign
  sign = 1.0;   // positive thermodiffusion coefficient by default

  if (plist.isParameter("Thermodiffusion Coefficient Sign"))
  {
    string tdcSign = plist.get<string>("Thermodiffusion Coefficient Sign");
    if (tdcSign == "Positive")
      sign = 1.0;
    else if (tdcSign == "Negative")
      sign = -1.0;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
          "Error ! Thermodiffusion Coefficient Sign must be either Positive or Negative !" << std::endl);
  }

  minTemp = plist.get<double>("Minimum Temperature");       // [K]
  maxTemp = plist.get<double>("Maximum Temperature");       // [K]
  minActE = plist.get<double>("Minimum Activation Energy"); // [eV]
  maxActE = plist.get<double>("Maximum Activation Energy"); // [eV]

  slope = (maxActE - minActE) / (minTemp - maxTemp);  // [eV/K]
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
ThermodiffCoeff_Custom<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->sublist("Thermodiffusion Coefficient ParameterList", false, "");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<std::string>("Value", "Custom", "Customized model for the thermodiffusion coefficient");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<std::string>("Thermodiffusion Coefficient Sign", "Positive", "Thermodiffusion coefficient can be either Positive or Negative");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<double>("Thermodiffusion Multiplier", 1.0, "[cm^2.eV/s]");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<double>("Minimum Temperature", 300.0, "[K]");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<double>("Maximum Temperature", 300.0, "[K]");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<double>("Minimum Activation Energy", 0.0, "[eV]");
  p->sublist("Thermodiffusion Coefficient ParameterList").set<double>("Maximum Activation Energy", 0.0, "[eV]");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
