

#ifndef CHARON_BC_STRATEGY_INTERFACE_SIMPLE_IMPL_HPP
#define CHARON_BC_STRATEGY_INTERFACE_SIMPLE_IMPL_HPP

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
#include "Panzer_DOF.hpp"
#include "Panzer_Sum.hpp"
//#include "Panzer_Product.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_FieldSpy.hpp"

#include "Charon_Names.hpp"

// Physics Block Suffixes
// +-------+-------+
// +   1   |   2   |
// +-------+-------+
//

template <typename EvalT>
charon::BCStrategy_Interface_Simple<EvalT>::
BCStrategy_Interface_Simple(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data)
  : panzer::BCStrategy_Interface_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Interface Simple" );
}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using std::string;
  using Teuchos::RCP;

  // Get the Data parameter list.
  {
    const RCP<Teuchos::ParameterList> dataPList = this->m_bc.nonconstParams();
    TEUCHOS_ASSERT(Teuchos::nonnull(dataPList));
    const char* coeff_names[] = {"a", "b", "c", "d"};
    {
      Teuchos::ParameterList vp;
      vp.set<string>("Coupling DOF Name", "ELECTRIC_POTENTIAL", "Field used for coupling at interface");
      vp.set<bool>("Field Spy", false, "Turn on field spy debugging?");
      for (int i = 0; i < 2; ++i)
        vp.set<double>(coeff_names[i], 1.0, "Coefficient");
      for (int i = 2; i < 4; ++i)
        vp.set<double>(coeff_names[i], 0.0, "Coefficient");
      dataPList->validateParametersAndSetDefaults(vp);
    }
    for (int i = 0; i < 4; ++i) {
      coeffs_[i] = dataPList->get<double>(coeff_names[i]);
    }
    coupling_dof_name_ = dataPList->get<string>("Coupling DOF Name");
    field_spy_ = dataPList->get<bool>("Field Spy");
  }

  // Are we setting up my side or the other side of the interface?
  const int di = this->getDetailsIndex();

  // my side's DOF name
  dof_name_ = di == 0 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();
  // other side's DOF name
  other_dof_name_ = di == 1 ? this->m_bc.equationSetName() : this->m_bc.equationSetName2();

  // get the physics block parameter list
  RCP<const Teuchos::ParameterList> pbParamList = side_pb.getParameterList();

  // get the equation set parameter list
  const Teuchos::ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any suffix parameters
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // unique name for the coupling DOF on this side (needed for the extra gather)
  my_coupling_dof_name_ = coupling_dof_name_ + discsuffix;

  // unique residual name
  const string residual_name = "Residual_" + this->m_bc.equationSetName();

  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);
  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name, dof_name_, "", integration_order, side_pb);
}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
setCombineValues(Teuchos::ParameterList& p,
                 const std::string value_name1, const double scalar1,
                 const std::string value_name2, const double scalar2,
                 const std::string value_name3, const double scalar3,
                 const std::string value_name4, const double scalar4)
{
  std::vector<std::string> values_names(2);
  values_names[0] = value_name1;
  values_names[1] = value_name2;
  if (value_name3 != "") values_names.push_back(value_name3);
  if (value_name4 != "") values_names.push_back(value_name4);
  p.set< Teuchos::RCP<std::vector<std::string> > >(
    "Values Names", Teuchos::rcp(new std::vector<std::string>(values_names)));
  std::vector<double> scalars(2);
  scalars[0] = scalar1;
  scalars[1] = scalar2;
  if (values_names.size() > 2) scalars.push_back(scalar3);
  if (values_names.size() > 3) scalars.push_back(scalar4);
  p.set< Teuchos::RCP<const std::vector<double> > >(
    "Scalars", Teuchos::rcp(new std::vector<double>(scalars)));
}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
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

  const string residual_name = std::get<0>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name_);

  // Are we setting up my side or the other side of the interface?
  const int di = this->getDetailsIndex();

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // build and register all closure models
  // pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // Implement the interface condition by substituting the rhs of the condition
  // into the flux term for f. The flux term is the natural BC on this side of
  // the element block, so it doesn't appear in the following evaluators, just
  // like in the implementation of the constant Neumann boundary condition.

  const string
    normal_name = (0 == di) ? "my_normal" : "other_normal",
    coupling_dof_grad_name = (0 == di)? "my_" + coupling_dof_name_ + "_gradient" : "other_" + coupling_dof_name_ + "_gradient",
    normal_dot_coupling_grad_name = (0 == di)? "my_normal_dot_coupling_grad" : "other_normal_dot_coupling_grad",
    other_normal_dot_coupling_grad_name = (0 == di)? "other_normal_dot_coupling_grad" : "my_normal_dot_coupling_grad",
    sum_contributions_name = (0 == di) ? "sum_contributions1" : "sum_contributions2";

  // Either get My DOF or Other DOF
  {
    const string plist_name = (0 == di) ? "My DOF" : "Other DOF";
    ParameterList p(plist_name);
    p.set("Name", dof_name_);
    p.set("Basis", basis);
    p.set("IR", ir);
    const RCP< PHX::Evaluator<panzer::Traits> >
      op = rcp(new panzer::DOF<EvalT,panzer::Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Either get My Side Normal or Other Side Normal
  {
    const string plist_name = (0 == di) ? "My Side Normal" : "Other Side Normal";
    ParameterList p(plist_name);
    p.set("Name", normal_name);
    p.set("Side ID", pb.cellData().side());
    p.set("IR", ir);
    p.set("Normalize", true);
    const RCP< PHX::Evaluator<panzer::Traits> >
      op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
   this->template registerEvaluator<EvalT>(fm, op);
  }

  // Either get My Side Grad Phi or Other Side Grad Phi
  {
    const string plist_name = (0 == di) ? "My Side Grad Phi" : "Other Side Grad Phi";
    ParameterList p(plist_name);
    p.set("Name", my_coupling_dof_name_);
    p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
    p.set("Gradient Name", coupling_dof_grad_name);
    p.set("Basis", basis);
    p.set("IR", ir);
    const RCP< PHX::Evaluator<panzer::Traits> >
      op = rcp(new panzer::DOFGradient<EvalT,panzer::Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Print the Electric Potential
  if(field_spy_) {
    const RCP<panzer::FieldSpy<EvalT,panzer::Traits> > opEval
      = rcp(new panzer::FieldSpy<EvalT,panzer::Traits>(my_coupling_dof_name_,basis->functional));
    const RCP< PHX::Evaluator<panzer::Traits> > op = opEval;
    this->template registerEvaluator<EvalT>(fm, op);

    fm.requireField<EvalT>(opEval->getRequiredFieldTag());
  }

  // Dot product
  {
    const string plist_name = (0 == di) ? "My Grad Phi Dot Normal" : "Other Grad Phi Dot Normal";
    ParameterList p(plist_name);
    p.set("Result Name", normal_dot_coupling_grad_name);
    p.set("Vector A Name", coupling_dof_grad_name);
    p.set("Vector B Name", normal_name);
    p.set("Point Rule", Teuchos::rcp_dynamic_cast<const panzer::PointRule>(ir));
    const RCP< PHX::Evaluator<panzer::Traits> >
      op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
    this->template registerEvaluator<EvalT>(fm, op);
  }

  // Print the Dot Product
  if(field_spy_) {
    const RCP<panzer::FieldSpy<EvalT,panzer::Traits> > opEval
      = rcp(new panzer::FieldSpy<EvalT,panzer::Traits>(normal_dot_coupling_grad_name,ir->dl_scalar));
    const RCP< PHX::Evaluator<panzer::Traits> > op = opEval;
    this->template registerEvaluator<EvalT>(fm, op);

    fm.requireField<EvalT>(opEval->getRequiredFieldTag());
  }

  // The interface condition is
  //   c dot(normal1, grad phi1) + d dot(normal2, grad phi2) + a dof1 + b dof2
  //
  // When di == 0, dof_name = dof1, other_dof_name = dof2
  // When di == 1, dof_name = dof2, other_dof_name = dof1
  //

  if (0 == di) {
    { // add contribution to the residual
      using panzer::EvaluatorStyle;
      using panzer::Integrator_BasisTimesScalar;
      using panzer::Traits;
      using PHX::Evaluator;
      const RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
        residual_name, sum_contributions_name, *basis, *ir));
      this->template registerEvaluator<EvalT>(fm, op);
    }
    { // compute the linear combination
      ParameterList p("a dof_me + b dof_other");
      p.set("Sum Name", sum_contributions_name);
      setCombineValues(p,
                       dof_name_, coeffs_[0],
                       other_dof_name_, coeffs_[1],
                       normal_dot_coupling_grad_name, coeffs_[2],
                       other_normal_dot_coupling_grad_name, coeffs_[3]);
      p.set("Data Layout", ir->dl_scalar);
      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::Sum<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }
  }

}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                               const panzer::PhysicsBlock& pb,
                                               const panzer::LinearObjFactory<panzer::Traits> & lof,
                                               const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // First do the standard gather of all the DOFs within this physics block
  pb.buildAndRegisterGatherAndOrientationEvaluators(fm, lof, user_data);

  // Get tangent fields
  const std::vector<panzer::StrPureBasisPair> tangent_fields = pb.getTangentFields();

  // We don't support tangent fields yet
  TEUCHOS_ASSERT( 0 == tangent_fields.size() );

  // Second do an extra gather for any continuous fields needed on both sides of the interface

  // Get the vector of residual contributions data
  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
    Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  // Charon assumes that there is only one residual contribution
  TEUCHOS_ASSERT( 1 == data.size() );

  // Get the basis from the tuple
  const RCP<const panzer::PureBasis> basis = std::get<4>(data[0]);

  // Register the extra gather evaluator
  {
    ParameterList p("Extra Gather");
    p.set("Basis", basis);
    RCP<std::vector<string> > dof_names = rcp(new std::vector<string>);
    RCP<std::vector<string> > indexer_names = rcp(new std::vector<string>);
    dof_names->push_back(my_coupling_dof_name_);
    indexer_names->push_back(coupling_dof_name_);
    p.set("DOF Names", dof_names);
    p.set("Indexer Names", indexer_names);
    p.set("Sensitivities Name", "");
    p.set("First Sensitivities Available", true);

    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);

    this->template registerEvaluator<EvalT>(fm, op);
  }

}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{
}

template <typename EvalT>
void charon::BCStrategy_Interface_Simple<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{
}
#endif
