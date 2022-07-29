// Input arguments
#include "Teuchos_ParameterList.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_Factory.H"
#include "NOX_Utils.H"
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_STK_NOXObserverFactory.hpp"
#include "Charon_CurrentConstraintList.hpp"
#include "Charon_NOXObserverFactory.hpp"

// Householder Constraint Solver headers
#include "NOX_TpetraTypedefs.hpp"
#include "LOCA_Tpetra_Factory.hpp"
#include "LOCA_Thyra_Group.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"
#include "LOCA_Tpetra_ConstraintModelEvaluator.hpp"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_Stepper.H"

namespace charon {

void solve_steadystate_constraint_problem(const Teuchos::RCP<Teuchos::ParameterList>& input_parameters,
                                          const Teuchos::RCP<Thyra::ModelEvaluator<double>>& thyra_me,
                                          charon::CurrentConstraintList& current_constraints,
                                          const std::vector<std::pair<std::string, std::string>>& contact_sides, // for response names
                                          const Teuchos::RCP<std::ostream>& serial_ostream,
                                          const Teuchos::RCP<std::ostream>& parallel_ostream,
                                          const panzer_stk::ModelEvaluatorFactory<double>& me_factory,
                                          const panzer_stk::NOXObserverFactory& nox_observer_factory)
{
  const bool print_debug = std::getenv("CHARON_PRINT_DEBUG");
  auto& out = *serial_ostream;
  auto& pout = *parallel_ostream;

  if (print_debug) {
    out << "***********\nBEGIN DEBUG: Householder solver\n***********\n";
    input_parameters->print(std::cout);
    out << "***********\n";
    out << current_constraints << std::endl;
    out << "***********\n";
    for (int i=0; i < thyra_me->Np(); ++i) {
      auto p_names = thyra_me->get_p_names(i);
      for (int j=0; j < p_names->size(); ++j)
        pout << "ME Parameter(" << i << ")=" << (*p_names)[j] << std::endl; 
    }
    for (int i=0; i < thyra_me->Ng(); ++i) {
      auto g_names = thyra_me->get_g_names(i);
      for (int j=0; j < g_names.size(); ++j)
        pout << "ME Response(" << i << ")=" << (g_names)[j] << std::endl; 
    }
    out << "***********\nEND DEBUG: Householder solver\n***********\n";
  }

  auto initial_guess = thyra_me->getNominalValues().get_x();

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess,
                                       thyra_me,
                                       Teuchos::null, Teuchos::null, Teuchos::null,
                                       false));
  
  Teuchos::RCP<LOCA::Abstract::Factory> tpetra_factory = Teuchos::rcp(new LOCA::Tpetra::Factory);
  
  auto top_params = Teuchos::sublist(input_parameters,"Solution Control",true);
  Teuchos::RCP<LOCA::GlobalData> loca_global_data = LOCA::createGlobalData(top_params, tpetra_factory);
  
  Teuchos::RCP<LOCA::ParameterVector> p_vec = Teuchos::rcp(new LOCA::ParameterVector);
  Teuchos::RCP<std::vector<std::string>> constraint_p_names = Teuchos::rcp(new std::vector<std::string>);
  for (int i=0; i < current_constraints.size(); ++i) {
    p_vec->addParameter(current_constraints[i]->parameterName(),
                        current_constraints[i]->initialVoltage());
    constraint_p_names->push_back(current_constraints[i]->parameterName());
  }

  // If doing a parameter sweep, add the loca continuation parameter
  const auto solver_type = top_params->get<std::string>("Piro Solver");
  if (solver_type == "LOCA-Constrained") {
    const auto sweep_param_name = top_params->sublist("LOCA").sublist("Stepper").get<std::string>("Continuation Parameter");
    p_vec->addParameter(sweep_param_name,0.0);
  }

  // Map Constraint params to model evaluator parameter indices
  std::vector<int> me_p_indices;
  for (int p=0; p < p_vec->length(); ++p) {
    bool found = false;
    for (int l=0; l < thyra_me->Np(); ++l) {
      auto names_vec = thyra_me->get_p_names(l);
      Teuchos::Array<std::string>::const_iterator it = std::find(names_vec->begin(),names_vec->end(),p_vec->getLabel(p));
      if (it != names_vec->end()) {
        me_p_indices.push_back(l);
        found = true;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!found,std::runtime_error,"ERROR: Constraint parameter \"" << p_vec->getLabel(p) <<  "\" does not exist in the ModelEvaluator!");
  }

  if (print_debug) {
    p_vec->print(out);
    for (int p=0; p < p_vec->length(); ++p) {
      out << "Charon_Main: me_p_indices(" << p << ")=" << me_p_indices[p] << std::endl;
    }
  }

  Teuchos::RCP<LOCA::Thyra::Group> loca_group = Teuchos::rcp(new LOCA::Thyra::Group(loca_global_data,
                                                                                    *nox_group,
                                                                                    *p_vec,
                                                                                    me_p_indices));

  auto g_names = Teuchos::rcp(new std::vector<std::string>);
  for (int i=0; i < current_constraints.size(); ++i) {
    g_names->push_back(current_constraints[i]->responseName());
  }
  auto x_thyra = ::Thyra::createMember(thyra_me->get_x_space(),"x");
  NOX::Thyra::Vector x(x_thyra);
  auto constraints = Teuchos::rcp(new LOCA::MultiContinuation::ConstraintModelEvaluator(thyra_me,*p_vec,*g_names,x));
  // Set initial parameter conditions
  constraints->setX(x);
  for (int i=0; i < current_constraints.size(); ++i)
    constraints->setParam(i,current_constraints[i]->initialVoltage());

  // Create the constraints list
  auto& locaParamsList = top_params->sublist("LOCA");
  auto& constraint_list = locaParamsList.sublist("Constraints");
  constraint_list.set("Bordered Solver Method", "Householder");
  constraint_list.set<Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>>("Constraint Object", constraints);
  constraint_list.set("Constraint Parameter Names", constraint_p_names);

  auto loca_parser = Teuchos::rcp(new LOCA::Parameter::SublistParser(loca_global_data));
  loca_parser->parseSublists(top_params);

  NOX::StatusTest::Factory stf;
  Teuchos::RCP<NOX::Utils> utils = Teuchos::rcp(new NOX::Utils(top_params->sublist("NOX").sublist("Printing")));
  auto status_tests = stf.buildStatusTests(top_params->sublist("NOX").sublist("Status Tests"),*utils);

  auto cnof = dynamic_cast<const charon::NOXObserverFactory&>(nox_observer_factory);
  const_cast<charon::NOXObserverFactory&>(cnof).setModelEvaluator(thyra_me);

  Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo = nox_observer_factory.buildNOXObserver(me_factory.getMesh(),
                                                                                           me_factory.getGlobalIndexer(),
                                                                                           me_factory.getLinearObjFactory());

  top_params->sublist("NOX").sublist("Solver Options").set("User Defined Pre/Post Operator", ppo);

  if (solver_type == "NOX-Constrained") {

    // Index into p_vec for constraint params. Constraints only!
    std::vector<int> param_ids(constraint_p_names->size());
    for (size_t i=0; i < constraint_p_names->size(); ++i)
      param_ids[i] = p_vec->getIndex((*constraint_p_names)[i]);
    auto constraint_list_ptr = Teuchos::rcpFromRef(constraint_list);

    Teuchos::RCP<LOCA::MultiContinuation::ConstrainedGroup> loca_constrained_group =
      Teuchos::rcp(new LOCA::MultiContinuation::ConstrainedGroup(loca_global_data,
                                                                 loca_parser,
                                                                 constraint_list_ptr,
                                                                 loca_group,
                                                                 constraints,
                                                                 param_ids,
                                                                 false));

    auto solver = NOX::Solver::buildSolver(loca_constrained_group, status_tests, Teuchos::sublist(top_params,"NOX"));
    NOX::StatusTest::StatusType solve_status = solver->solve();
    TEUCHOS_TEST_FOR_EXCEPTION(solve_status != NOX::StatusTest::Converged,std::runtime_error,
                               "ERROR: Nonlinear solver with Householder Constraints Failed to Converge!");
  }
  else {
    // Hard code these for continuation with Householder constraints
    // in case arc-length is used. This will flatten the linear solve
    // blocks for arc-length with the constraints.
    auto stepper_list = top_params->sublist("LOCA").sublist("Stepper");
    stepper_list.set("Bordered Solver Method", "Nested");
    stepper_list.sublist("Nested Bordered Solver").set("Bordered Solver Method", "Householder");
    
    LOCA::Stepper stepper(loca_global_data,loca_group,status_tests,top_params);
    auto solve_status = stepper.run();
    TEUCHOS_TEST_FOR_EXCEPTION(solve_status != LOCA::Abstract::Iterator::Finished,std::runtime_error,
                               "ERROR: Parameter sweep with Householder Constraints Failed to Converge!");
  }

  LOCA::destroyGlobalData(loca_global_data); // Breaks RCP circular dependency

}

} // namespace charon
