
#ifndef CHARON_DIFFCOEFF_DEFAULT_IMPL_HPP
#define CHARON_DIFFCOEFF_DEFAULT_IMPL_HPP

#include <cmath>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"

/*
The default diffusion coefficient is determined by the Einstein relation, i.e.,
D/u = kbT/q (in physical units). After introducing scaling factors, the scaled
relaltion becomes Ds/us = Ts, where s subscript denotes scaled. The scaled relation
makes use of the fact that V0 = kbT0/q and u0 = D0/V0 with 0 denoting scaling
parameters.

Note that the Einstein relation always holds for D. Schroeder's DD formulation,
independent of the carrier statistics.

When bUseFD = true and FDFormula = Diffusion, D/u = kbT/q * F1/2(eta)/F(-1/2)(eta).
The scaled version is Ds/us = Ts * F1/2(eta)/F(-1/2)(eta).

When Carrier Type = Ion, always use the D/u = kbT/q Einstein relation.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DiffCoeff_Default<EvalT, Traits>::
DiffCoeff_Default(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
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
  carrType = p.get<string>("Carrier Type");
  isEdgedl = p.get<bool>("Is Edge Data Layout");
  bUseFD = p.get<bool>("Fermi Dirac");
  fdFormula = p.get<string>("FD Formula");

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
  if (carrType == "Electron")
    diffcoeff = MDField<ScalarT,Cell,Point>(n.field.elec_diff_coeff,output_dl);
  else if (carrType == "Hole")
    diffcoeff = MDField<ScalarT,Cell,Point>(n.field.hole_diff_coeff,output_dl);
  else if (carrType == "Ion")
    diffcoeff = MDField<ScalarT,Cell,Point>(n.field.ion_diff_coeff,output_dl);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  this->addEvaluatedField(diffcoeff);

  // dependent fields
  if (carrType == "Electron")
    mobility = MDField<const ScalarT,Cell,Point>(n.field.elec_mobility,output_dl);
  else if (carrType == "Hole")
    mobility = MDField<const ScalarT,Cell,Point>(n.field.hole_mobility,output_dl);
  else if (carrType == "Ion")
    mobility = MDField<const ScalarT,Cell,Point>(n.field.ion_mobility,output_dl);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;

  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  this->addDependentField(mobility);
  this->addDependentField(latt_temp);

  // dependent fields for bUseFD = true and fdFormula = Diffusion
  if ( bUseFD && (fdFormula == "Diffusion") )
  {
    if (carrType == "Electron")
    {
      carr_dens = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
      eff_dos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
    }
    else if (carrType == "Hole")
    {
      carr_dens = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
      eff_dos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);
    }
    this->addDependentField(carr_dens);
    this->addDependentField(eff_dos);
  }

  std::string name = "Diffusion_Coefficient_Default";
  this->setName(name);

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral = Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));
  forwardFermiIntegral = Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::forward_MinusOneHalf));

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DiffCoeff_Default<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  //----------------------------------------------------------------------------
  // evaluate diff. coeff. at the center of a primary edge
  //----------------------------------------------------------------------------
  if (isEdgedl)
  {
   // use the modified D/u relation
   if ( bUseFD && (fdFormula == "Diffusion") )
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        const ScalarT& mob = mobility(cell,edge);  // mobility has values at primary edge

        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        // obtain values at the center of a primary edge
        ScalarT latt = (latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0;

        // latt should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0/T0;

        ScalarT dens = (carr_dens(cell,node0) + carr_dens(cell,node1)) / 2.0;
        ScalarT dos = (eff_dos(cell,node0) + eff_dos(cell,node1)) / 2.0;
        ScalarT ratio = dens / dos;  // C0 is cancelled out

        // avoid to pass negative argument to the inverse function
        if (ratio <= 1e-4)  // similar treatment as in charon::Degeneracy_Factor
          diffcoeff(cell,edge) = mob * latt;
        else
        {
          ScalarT eta = (*inverseFermiIntegral)(ratio);    // F^(-1)_1/2(x)
          ScalarT fdint = (*forwardFermiIntegral)(eta);    // F_{-1/2}(x)
          diffcoeff(cell,edge) = mob * latt * ratio / fdint;

        }
      }
    }
   }   // end of if (bUseFD && fdFormula) block

   // use the Einstein D/u relation
   else
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        const ScalarT& mob = mobility(cell,edge);  // mobility has values at primary edge

        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        // latt_temp has values at IP or BASIS points, not at primary edge
        ScalarT latt = (latt_temp(cell,node0) + latt_temp(cell,node1)) / 2.0;

        // latt should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0/T0;

        diffcoeff(cell,edge) = mob * latt;

      }
    }
   }  // end of else (bUseFD && fdFormula) block

  }  // end of if (isEdgedl) block


  //----------------------------------------------------------------------------
  // evaluate diff. ceoff. at IP or BASIS points
  //----------------------------------------------------------------------------
  else
  {
   // use the modified D/u relation
   if ( bUseFD && (fdFormula == "Diffusion") )
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        ScalarT mob = mobility(cell,point);
        ScalarT latt = latt_temp(cell,point);
        ScalarT dos = eff_dos(cell,point);
        ScalarT dens = carr_dens(cell,point);
        ScalarT ratio = dens/dos;

        // latt should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0/T0;

        // avoid to pass negative argument to the inverse function
        if (ratio <= 1e-4)  // similar treatment as in charon::Degeneracy_Factor
          diffcoeff(cell,point) = mob * latt;  // nondegenerate
        else
        {
          ScalarT eta = (*inverseFermiIntegral)(ratio);    // F^(-1)_1/2(x)
          ScalarT fdint = (*forwardFermiIntegral)(eta);    // F_{-1/2}(x)
          diffcoeff(cell,point) = mob * latt * ratio / fdint;
        }

      }
    }
   }  // end of if (bUseFD && fdFormula) block

   // use the Einstein D/u relation
   else
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        ScalarT mob = mobility(cell,point);
        ScalarT latt = latt_temp(cell,point);

        // latt should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latt) <= 0.0)  latt = 300.0/T0;

        diffcoeff(cell,point) = mob * latt;  // scaled

        //if (carrType == "Ion")
        //  std::cout << "cell=" << cell << ", point=" << point << ", latt=" << latt << ", mob=" << mob << ", diff=" << diffcoeff(cell,point) << std::endl;
      }
    }
   }  // end of else (bUseFD && fdFormula) block

  }  // end of else (isEdgedl) block

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DiffCoeff_Default<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Is Edge Data Layout", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<bool>("Fermi Dirac", false);

  p->set<std::string>("FD Formula", "?");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
