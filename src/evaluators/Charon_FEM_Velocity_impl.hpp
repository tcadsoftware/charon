
#ifndef CHARON_FEM_VELOCITY_IMPL_HPP
#define CHARON_FEM_VELOCITY_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Charon_Names.hpp"

/*
Compute carrier velocity at IPs including the BGN effect when BGN = On, i.e.,
vn = -mun*Fn for electrons and vp = mup*Fp for holes. The carrier velocity is
used in computing the FEM stab. parameter 'tau'.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
FEM_Velocity<EvalT, Traits>::
FEM_Velocity(
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

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  carrType = p.get<string>("Carrier Type");
  isEdgedl = false;  // default
  if (p.isParameter("Is Edge Data Layout"))
    isEdgedl = p.get<bool>("Is Edge Data Layout");

  // Carrier dependent fields
  if (carrType == "Electron")
  {
    velocity = MDField<ScalarT,Cell,Point,Dim>(n.field.elec_velocity,vector);
    efield = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield,vector);
    mobility = MDField<const ScalarT,Cell,Point>(n.field.elec_mobility,scalar);
    sign = -1.0;
    this->addEvaluatedField(velocity);
    this->addDependentField(efield);
    this->addDependentField(mobility);
  }

  else if (carrType == "Hole")
  {
    velocity = MDField<ScalarT,Cell,Point,Dim>(n.field.hole_velocity,vector);
    efield = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield,vector);
    mobility = MDField<const ScalarT,Cell,Point>(n.field.hole_mobility,scalar);
    sign = 1.0;
    this->addEvaluatedField(velocity);
    this->addDependentField(efield);
    this->addDependentField(mobility);
  }

  else if (carrType == "Ion")
  {
    if (isEdgedl)
    {
      RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
      basis_name = basis->name();
      RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
      RCP<DataLayout> basis_scalar = basis->functional;
      RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
      num_edges = edge_scalar->dimension(1);
      cellType = cellTopoInfo->getCellTopology();

      edge_velocity = MDField<ScalarT,Cell,Point>(n.field.ion_velocity,edge_scalar);
      mobility = MDField<const ScalarT,Cell,Point>(n.field.ion_mobility,edge_scalar);
      potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,basis_scalar);
      this->addEvaluatedField(edge_velocity);
      this->addDependentField(potential);
      this->addDependentField(mobility);
    }
    else
    {
      velocity = MDField<ScalarT,Cell,Point,Dim>(n.field.ion_velocity,vector);
      efield = MDField<const ScalarT,Cell,Point,Dim>(n.field.ion_efield,vector);
      mobility = MDField<const ScalarT,Cell,Point>(n.field.ion_mobility,scalar);
      this->addEvaluatedField(velocity);
      this->addDependentField(efield);
      this->addDependentField(mobility);
    }

    int ion_charge = p.get<int>("Ion Charge"); // retrieve the integer ion charge
    if (ion_charge > 0)  // positive charge, similar to hole
      sign = 1.0;
    else if (ion_charge < 0) // negative charge, similar to electron
      sign = -1.0;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Ion Charge cannot be zero !");
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole or Ion !");

  std::string name = "FEM_Carrier_Velocity";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_Velocity<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if ((carrType == "Ion") && isEdgedl)
    basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_Velocity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  if ((carrType == "Ion") && isEdgedl)  // compute ion velocity at primary edges
  {
    double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

        // get local node coordinates
        x0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,0);
        x1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,0);
        if (num_dims > 1)  // 2D or 3D
        {
          y0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,1);
          y1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,1);
        }
        if (num_dims > 2)  // 3D
        {
          z0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,2);
          z1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,2);
        }

        // compute the primary cell edge length [scaled]
        // double edgeLen = std::sqrt(std::pow((x0-x1),2.0) + std::pow((y0-y1),2.0) + std::pow((z0-z1),2.0) );
        double edgeLen = std::sqrt(pow((x0-x1),2.0) + pow((y0-y1),2.0) + pow((z0-z1),2.0) );

        // compute the field [scaled] at the primary edge
        ScalarT edgeField = (potential(cell,node0) - potential(cell,node1)) / edgeLen;

        // obtain the ion mobility at the primary edge
        ScalarT edgeMob = mobility(cell,edge);

        // compute the ion velocity at the primary edge
        edge_velocity(cell,edge) = sign * edgeMob * edgeField;
      }
    }
  }  // end of if ((carrType == "Ion") && isEdgedl)

  else  // compute velocity at IPs
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t point = 0; point < num_points; ++point)
        for (std::size_t dim = 0; dim < num_dims; ++dim)
          velocity(cell,point,dim) = sign * mobility(cell,point) * efield(cell,point,dim);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
FEM_Velocity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<int>("Ion Charge", 1, "The integer number of ion charge");
  p->set<bool>("Is Edge Data Layout", false);

  return p;
}

}

#endif
