
#ifndef CHARON_THERMALCONDUCT_TEMPDEP_IMPL_HPP
#define CHARON_THERMALCONDUCT_TEMPDEP_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"


/*
The ion/vacancy soret coefficient is modeled as S = Ua/(kB*T^2), or -Ua/(kB*T^2) ?
to be tested and verified ?
where Ua is the potential energy barrier in [eV], and T is the lattice temperature in [K].
S is in the unit of [1/K], so the scaling parameter for S is 1/T0.

The Ua parameter can be changed by user in the input file or come from
charon::Material_Properties when not specified/changed.

Specification of the soret coefficient model in the input file takes the form of

<ParameterList name="Ion Soret Coefficient">
  <Parameter name="Value" type="string" value="TempDep" />
  <Parameter name="Soret Energy Barrier" type="double" value="1.2" />
  <Parameter name="Soret Coefficient Sign" type="string" value="Positive" />
</ParameterList>
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SoretCoeff_TempDep<EvalT, Traits>::
SoretCoeff_TempDep(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using Teuchos::ParameterList;
  using panzer::BasisIRLayout;

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

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Soret Coefficient ParameterList
  const ParameterList& plist = p.sublist("Soret Coefficient ParameterList");

  // Initialize the soret coefficient model parameters
  initialize(matName, plist);

  // Evaluated fields
  soret_coeff = MDField<ScalarT,Cell,Point>(n.field.ion_soret_coeff,output_dl);
  this->addEvaluatedField(soret_coeff);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.dof.latt_temp,scalar);
  this->addDependentField(latt_temp);

  std::string name = "Soret_Coefficient_TempDep";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SoretCoeff_TempDep<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;

  // obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;       // Boltzmann constant in [eV/K]

  if (isEdgedl) // calculate at edge midpoints
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        const ScalarT& scaledT = (latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0;   // scaled
        ScalarT lattT = T0*scaledT;  // in [K]

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

        ScalarT sc = sign * Ua / (kb*lattT*lattT);  // in [1/K]
        soret_coeff(cell,edge) = sc * T0;  // scaled
      }
    }
  }
  else // calculate at IP
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        const ScalarT& scaledT = latt_temp(cell,point);  // scaled
        ScalarT lattT = T0*scaledT;  // in [K]

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(lattT) <= 0.0)  lattT = 300.0;

        ScalarT sc = sign * Ua / (kb*lattT*lattT);  // in [1/K]
        soret_coeff(cell,point) = sc * T0;  // scaled
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
void SoretCoeff_TempDep<EvalT, Traits>::initialize
(const std::string& matName, const Teuchos::ParameterList& plist)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Retrieve parameters from xml
  if (plist.isParameter("Soret Energy Barrier"))
    Ua = plist.get<double>("Soret Energy Barrier");

  else  // otherwise, retrieve parameters from charon::Material_Properties
    Ua = matProperty.getPropertyValue(matName, "Soret Energy Barrier");

  // Soret Coefficient Sign
  sign = 1.0;   // positive Soret coefficient by default

  if (plist.isParameter("Soret Coefficient Sign"))
  {
    string scSign = plist.get<string>("Soret Coefficient Sign");
    if (scSign == "Positive")
      sign = 1.0;
    else if (scSign == "Negative")
      sign = -1.0;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl <<
          "Error ! Soret Coefficient Sign must be either Positive or Negative !" << std::endl);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SoretCoeff_TempDep<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->sublist("Soret Coefficient ParameterList", false, "");
  p->sublist("Soret Coefficient ParameterList").set<std::string>("Value", "TempDep", "Temperature-dependent Soret coefficient");
  p->sublist("Soret Coefficient ParameterList").set<double>("Soret Energy Barrier", 1.2, "[eV]");
  p->sublist("Soret Coefficient ParameterList").set<std::string>("Soret Coefficient Sign", "Positive", "Soret coefficient can be either Positive or Negative");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
