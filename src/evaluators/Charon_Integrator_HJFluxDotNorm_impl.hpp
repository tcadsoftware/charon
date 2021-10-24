
#ifndef CVFEM_INTEGRATOR_HJFLUXDOTNORM_IMPL_HPP
#define CVFEM_INTEGRATOR_HJFLUXDOTNORM_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Teuchos_TestForException.hpp"
#include "Charon_Names.hpp"

namespace charon {

/**
 * @brief This evaluator integrates the normal heterojunction (HJ) current density
 * J_{HJ}, where J_{HJ} = J_{TE} * (1 + tunnel), over a HJ, and adds the integrated
 * value to the corresponding dof residual for the CVFEM-SG scheme.
*/

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Integrator_HJFluxDotNorm<EvalT, Traits>::
Integrator_HJFluxDotNorm(
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

  // Obtain the BASIS layout
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> node_scalar = basis->functional;
  basis_name = basis->name();

  // Integration rule
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_scalar->dimension(1);

  // Obtain the Edge data
  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
  RCP<DataLayout> edge_scalar = cellTopoInfo->edge_scalar;
  RCP<DataLayout> edge_vector = cellTopoInfo->edge_vector;
  num_dims = edge_vector->dimension(2);

  // Get the primary cell topology
  cellType = cellTopoInfo->getCellTopology();

  // Evaluated field
  residual = MDField<ScalarT,Cell,BASIS>(p.get<string>("Residual Name"),node_scalar);
  this->addEvaluatedField(residual);

  // Dependent field
  hj_flux = MDField<const ScalarT,Cell,IP>(p.get<string>("Flux Name"),ip_scalar);
  this->addDependentField(hj_flux);

  // get multiplier
  multiplier = p.get<double>("Multiplier");

  // optional field multiplier
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"))
  {
    const std::vector<std::string>& field_multiplier_names =
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin();
      name != field_multiplier_names.end(); ++name)
    {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  // get string names
  residual_name = p.get<string>("Residual Name");
  flux_name = p.get<string>("Flux Name");

  std::string name = "Integrator_HJFluxDotNorm";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_HJFluxDotNorm<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0], this->wda);
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0], this->wda);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_HJFluxDotNorm<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  residual.deep_copy(ScalarT(0.0));

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
     int subcell_ind = this->wda(workset).subcell_index;
     unsigned nodes_per_subcell = cellType->getVertexCount(num_dims-1,subcell_ind);

     // std::cout << "cell=" << cell << ", subcell_ind=" << subcell_ind << ", nodes_per_subcell=" << nodes_per_subcell;
     // std::cout << ", residual_name = " << residual_name << ", flux_name = " << flux_name << std::endl;

     for (unsigned inode = 0; inode < nodes_per_subcell; ++inode)
     {
        // cell node number for node on subcell edge/face associated with heterojunction
        int cell_node = cellType->getNodeMap(num_dims-1,subcell_ind,inode);

        const panzer::IntegrationValues2<double> &iv = *this->wda(workset).int_rules[int_rule_index];

        ScalarT tmpVar = 1.0;
        for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
             field != field_multipliers.end(); ++field)
           tmpVar = tmpVar * (*field)(cell,inode);

        // For now this is assuming lowest order CVFEM with one integration point per node
        // on the boundary face.
        residual(cell,cell_node) = multiplier * tmpVar * hj_flux(cell,inode) * iv.weighted_measure(cell,inode);

        // std::cout << std::setprecision(20) << "cell=" << cell << ", cell_node=" << cell_node << ", inode=" << inode << ", residual=" << residual(cell,cell_node) << std::endl;

     }

  } // end cell loop
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Integrator_HJFluxDotNorm<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Flux Name", "?");
  p->set("Multiplier", 1.0);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);

  return p;
}

}

#endif
