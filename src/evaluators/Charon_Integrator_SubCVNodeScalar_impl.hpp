
#ifndef CHARON_INTEGRATOR_SUBCVNODESCALAR_T_HPP
#define CHARON_INTEGRATOR_SUBCVNODESCALAR_T_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Shards_CellTopology.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Integrator_SubCVNodeScalar<EvalT, Traits>::
Integrator_SubCVNodeScalar(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // Obtain the BASIS layout
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> node_scalar = basis->functional;
  basis_name = basis->name();
  num_nodes = node_scalar->dimension(1);

  // Obtain integration rule
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_scalar->dimension(1);

  // Get input parameters
  multiplier = p.get<double>("Multiplier");
  WithInterpolation = p.get<bool>("WithInterpolation");

  // Evaluated field
  residual = MDField<ScalarT,Cell,Point>(p.get<string>("Residual Name"), node_scalar);
  this->addEvaluatedField(residual);

  // Dependent field
  if (WithInterpolation) // dependent field is located at BASIS
    value = MDField<const ScalarT,Cell,Point>(p.get<string>("Value Name"), node_scalar);
  else // dependent field is located at IPs
    value = MDField<const ScalarT,Cell,Point>(p.get<string>("Value Name"), ip_scalar);
  this->addDependentField(value);

  std::string name = "Integrator_SubCV_NodeScalar: " + residual.fieldTag().name();
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_SubCVNodeScalar<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_SubCVNodeScalar<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  residual.deep_copy(ScalarT(0.0));

  // loop over cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    Kokkos::DynRankView<ScalarT,PHX::Device> cptsValues;

    // interpolate values at nodes to centroids of subcv
    if (WithInterpolation)  
    {
      cptsValues = createDynRankView(residual.get_static_view(),"cptsValues",num_ips);
      Kokkos::deep_copy(cptsValues,ScalarT(0.0));

      for (int inode = 0; inode < num_nodes; ++inode)
        for (int ip = 0; ip < num_ips; ++ip)
          cptsValues(ip) += (workset.bases[basis_index])->basis_scalar(cell, inode, ip) * value(cell, inode);
    }

    // integrate a scalar over a subcv and scatter to the nodal residual
    for (int node = 0; node < num_nodes; ++node)
    {
        ScalarT measure = (workset.int_rules[int_rule_index])->weighted_measure(cell,node);
        if (WithInterpolation)
          residual(cell, node) = multiplier * measure * cptsValues(node);
        else
          residual(cell, node) = multiplier * measure * value(cell, node);
    }

  } // end of cell loop
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Integrator_SubCVNodeScalar<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Residual Name", "??");
  p->set<std::string>("Value Name", "??");
  p->set("Multiplier", 1.0);
  p->set("WithInterpolation", true);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Names> c_names = Teuchos::rcp(new charon::Names(1,"","",""));
  p->set<Teuchos::RCP<const charon::Names> >("Names", c_names);

  return p;
}

}

#endif
