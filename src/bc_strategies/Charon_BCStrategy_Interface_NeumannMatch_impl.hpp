

#ifndef CHARON_BC_STRATEGY_INTERFACE_NEUMANNMATCH_IMPL_HPP
#define CHARON_BC_STRATEGY_INTERFACE_NEUMANNMATCH_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BasisIRLayout.hpp"

// Evaluators
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"

#include "Charon_Names.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Interface_NeumannMatch<EvalT>::
BCStrategy_Interface_NeumannMatch(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Interface_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Interface Neumann Match" );
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_NeumannMatch<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using std::string;
  using Teuchos::RCP;

  // are we setting up the left or right side of the interface?
  const int di = this->getDetailsIndex();

  // obtain the dof names
  const string dof_name = (di == 0) ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  // unique residual name
  const string residual_name = "Residual_" + this->m_bc.equationSetName();
  const string flux_name = "Other_Flux";

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);

  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);

}

// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_NeumannMatch<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
                           const Teuchos::ParameterList& /* models */,
                           const Teuchos::ParameterList& /* user_data */) const
{

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
        Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  string residual_name = std::get<0>(data[0]);
  string dof_name = std::get<1>(data[0]);
  string flux_name = std::get<2>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name);

  // are we setting up the left or right side of the interface?
  const int di = this->getDetailsIndex();

  if( di == 0 ) {
    { // add contribution to the residual
      using panzer::EvaluatorStyle;
      using panzer::Integrator_BasisTimesScalar;
      using panzer::Traits;
      using PHX::Evaluator;
      const RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
        residual_name, flux_name, *basis, *ir, -1));
      fm.template registerEvaluator<EvalT>(op);
    }
  } else {
    const std::string dof_grad_name = dof_name + "_gradient";
    { // Compute side 2 normal.
      ParameterList p("Side Normal");
      p.set("Name", "Other_Side_Normal");
      p.set("Side ID", pb.cellData().side());
      p.set("IR", ir);
      p.set("Normalize", true);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
    { // Compute dof side 2 gradient.
      ParameterList p("Other DOF gradient");
      p.set("Name", dof_name);
      p.set("Gradient Name", dof_grad_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
    { // Compute dot(dof side 2 gradient, side 2 normal).
      ParameterList p("dot(Other DOF gradient, other normal)");
      p.set("Result Name", flux_name);
      p.set("Vector A Name", dof_grad_name);
      p.set("Vector B Name", "Other_Side_Normal");
      p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
  }
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_NeumannMatch<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::PhysicsBlock& pb,
                                               const panzer::LinearObjFactory<panzer::Traits> & lof,
                                               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // Gather
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm,lof,user_data);

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_NeumannMatch<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Interface_NeumannMatch<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
