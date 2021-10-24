
#ifndef CHARON_THERMODIFFCOEFF_DEFAULT_IMPL_HPP
#define CHARON_THERMODIFFCOEFF_DEFAULT_IMPL_HPP

#include <cmath>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

/*
The model implements the ion thermodiffusion coefficient as
D_{Ts} = u_{s} * T_{s} * S_{s} with all quantities in scaled units.

This model is used when "Ion Thermodiffusion Coefficient" is not specified in the
input xml file for DDIonLattice simulations. It is applicable to ions/vacancies only,
neither electrons nor holes.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
ThermodiffCoeff_Default<EvalT, Traits>::
ThermodiffCoeff_Default(
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

  // dependent fields
  mobility = MDField<const ScalarT,Cell,Point>(n.field.ion_mobility,output_dl);
  soret_coeff = MDField<const ScalarT,Cell,Point>(n.field.ion_soret_coeff,output_dl);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);

  this->addDependentField(mobility);
  this->addDependentField(soret_coeff);
  this->addDependentField(latt_temp);

  std::string name = "Thermodiffusion_Coefficient_Default";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
ThermodiffCoeff_Default<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  if (isEdgedl) // evaluate thermodiff. coeff. at edge midpoints
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        const ScalarT& latt = (latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0;   // scaled
        const ScalarT& mob = mobility(cell,edge);   // scaled
        const ScalarT& soret = soret_coeff(cell,edge);
        thermodiff_coeff(cell,edge) = mob * latt * soret;
      }
    }
  }

  else  // evaluate thermodiff. coeff. at IP or BASIS points
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        const ScalarT& mob = mobility(cell,point);   // scaled
        const ScalarT& soret = soret_coeff(cell,point);
        const ScalarT& latt = latt_temp(cell,point);
        thermodiff_coeff(cell,point) = mob * latt * soret;
      }
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
ThermodiffCoeff_Default<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  return p;
}

}

#endif
