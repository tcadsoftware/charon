
#ifndef CHARON_BCSTRATEGY_DIRICHLET_SCHOTTKYCONTACT_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_SCHOTTKYCONTACT_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_GlobalData.hpp"

#include "Charon_Names.hpp"
#include "Charon_BCStrategy_Dirichlet_SchottkyContact.hpp"

#include "Charon_Scaling_Parameters.hpp"

// evaluators
#include "Charon_BC_DirichletSchottkyContact.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_SchottkyContact<EvalT>::
BCStrategy_Dirichlet_SchottkyContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Dirichlet Schottky Contact" );
}

template <typename EvalT>
void charon::BCStrategy_Dirichlet_SchottkyContact<EvalT>::
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
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? 
    eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? 
    eqSetPList.get<string>("Discontinuous Suffix") : "";

  // get the name of electrostatic potential
  m_names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  string phi_name = m_names->dof.phi;

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // check if Equation Set Name = ALL_DOFS or ELECTRIC_POTENTIAL
  string eqsn = this->m_bc.equationSetName();
  if (eqsn != phi_name)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Equation Set Name must be either ALL_DOFS or "
        << phi_name << ". But you entered \"" << eqsn << "\" ! \n");

  // get dofs
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
    dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    if (dof_it->first == phi_name) { // Dirichlet, NLP, ELECTRIC_POTENTIAL 
      // need the dof name to form the residual
      this->required_dof_names.push_back(dof_it->first);

      // unique residual name
      string residual_name = "Residual_" + dof_it->first;

      // map residual to dof
      this->residual_to_dof_names_map[residual_name] = dof_it->first;

      // map residual to target field
      this->residual_to_target_field_map[residual_name] = "Target_" + dof_it->first;

      // find the basis for this dof
      this->basis = dof_it->second;
    } 
  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
      "Error: the name \"" << this->m_bc.equationSetName()
      << "\" is not a valid DOF for the boundary condition:\n"
      << this->m_bc << "\n");
  

  
}


template <typename EvalT>
void charon::BCStrategy_Dirichlet_SchottkyContact<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& side_pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                           const Teuchos::ParameterList& models,
                           const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // build and register all closure models
  side_pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // provide a target value to map into residual
  {
    RCP<const charon::Names> names = this->m_names;
    Teuchos::RCP<charon::Scaling_Parameters> scaleParams = 
      user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

    // build the BC_SchottkyContact evaluator 
    ParameterList p("BC Dirichlet Schottky Contact");
    p.set("Prefix", "Target_");
    p.set("Field Library", side_pb.getFieldLibraryBase());
    p.set("Names", names);
    p.set("Scaling Parameters", scaleParams);
    if (this->m_bc.params()->isParameter("Voltage")) {
      p.setEntry("Voltage", this->m_bc.params()->getEntry("Voltage"));
    } else if (this->m_bc.params()->isParameter("Varying Voltage")) {
      p.setEntry("Varying Voltage",
		 this->m_bc.params()->getEntry("Varying Voltage"));
    } else {
      p.set<double>("Voltage", 0);
    }
    if (this->m_bc.params()->isParameter("Work Function")) {
      p.setEntry("Work Function", this->m_bc.params()->getEntry("Work Function"));
    } 

    p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::BC_DirichletSchottkyContact<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }
}




#endif
