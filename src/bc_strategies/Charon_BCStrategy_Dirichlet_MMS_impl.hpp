
#ifndef CHARON_BCSTRATEGY_DIRICHLET_MMS_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_MMS_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PureBasis.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

// Evaluators
#include "Charon_MMS_AnalyticSolutions.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_MMS<EvalT>::
BCStrategy_Dirichlet_MMS(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "MMS");
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_MMS<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::vector;
  using std::string;
  using std::pair;

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = side_pb.getParameterList();
  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // get MMS type
  if ( dataPList->isParameter("MMS Type") )
  {
    paramName.push_back("MMS Type");
    value = dataPList->template get<string>("MMS Type");
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: You must provide MMS Type for BC, "
        << "either MMS_RDH_1 or MMS_RDH_2 !" );
  }

  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  const charon::Names& n = *names;

  // set dirichlet conditions on all dofs (potential and carrier density)
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
   dofs.begin(); dof_it != dofs.end(); ++dof_it)
  {
    // for the isothermal DD formulations
    if ( (dof_it->first == n.dof.phi) || (dof_it->first == n.dof.edensity) ||
         (dof_it->first == n.dof.hdensity) )
    {
      // need the dof values to form the residual
      this->required_dof_names.push_back(dof_it->first);

      // unique residual name
      std::string residual_name = "Residual_" + dof_it->first;

      // map residual to dof
      this->residual_to_dof_names_map[residual_name] = dof_it->first;

      // map residual to target field
      this->residual_to_target_field_map[residual_name] = "Target_" + dof_it->first;

      // For now assume that potential and carrier density use the same basis
      this->basis = dof_it->second;
    }

  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
                             "Error: the name \"" << this->m_bc.equationSetName()
                             << "\" is not a valid DOF for the boundary condition:\n"
                             << this->m_bc << "\n");

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_MMS<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
                           const Teuchos::ParameterList& /* models */,
                           const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // build and register all closure models
  //pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // get the element block id and physics id for the given physics block
  string ebID_pb = pb.elementBlockID();
  string pbID = pb.physicsBlockID();

  // get the element block ID for the given boundary condition
  string ebID_bc = this->m_bc.elementBlockID();

  if (ebID_pb != ebID_bc)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: " << pbID << " corresponds to "
      << ebID_pb << ", while the BC corresponds to " << ebID_bc << "! \n");


  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // allow only one equation set per physics block
  if (pbParamList->numParams() > 1)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
       "The physics block " << pbParamList->name() << " has more than one equation sets ! ");

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));

  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

 // ParameterList empty_list;
  ParameterList input;

  if (value == "MMS_RDH_1")
  {
  // build the BC_MMS evaluator for MMS_RDH_1
  input.set("Scaling Parameters", scaleParams);
  RCP< PHX::Evaluator<panzer::Traits> > op =
    rcp(new charon::MMS_DD_RDH_1_AnalyticSolution<EvalT,panzer::Traits>("Target_",*names,pb.getFieldLibraryBase(),Teuchos::null,input));
  fm.template registerEvaluator<EvalT>(op);
  }
  else if (value == "MMS_RDH_2")
  {
  // build the BC_MMS evaluator for MMS_RDH_2
  input.set("Scaling Parameters", scaleParams);
  RCP< PHX::Evaluator<panzer::Traits> > op =
    rcp(new charon::MMS_DD_RDH_2_AnalyticSolution<EvalT,panzer::Traits>("Target_",*names,pb.getFieldLibraryBase(),Teuchos::null,input));
  fm.template registerEvaluator<EvalT>(op);
  }

}

#endif
