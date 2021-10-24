
#ifndef CHARON_BCSTRATEGY_DIRICHLET_THERMALCONTACT_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_THERMALCONTACT_IMPL_HPP

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
#include "Charon_BC_ThermalContact.hpp"

#include "Charon_Scaling_Parameters.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_ThermalContact<EvalT>::
BCStrategy_Dirichlet_ThermalContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Thermal Contact" );
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_ThermalContact<EvalT>::
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

  // get the name of lattice temperature
  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  string temp_name = names->dof.latt_temp;

  // check if Equation Set Name = ALL_DOFS or Lattice Temperature
  string eqsn = this->m_bc.equationSetName();
  if ( (eqsn != "ALL_DOFS") && (eqsn != temp_name) )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Equation Set Name must be either ALL_DOFS or "
        << temp_name << ". But you entered \"" << eqsn << "\" ! \n");

  // get dofs
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
    dofs.begin(); dof_it != dofs.end(); ++dof_it)
  {
    if (dof_it->first == temp_name)
    {
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


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_ThermalContact<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                           const Teuchos::ParameterList& models,
                           const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // build and register all closure models
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  double value = this->m_bc.params()->template get<double>("Temperature");

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // build the BC_ThermalContact evaluator
  ParameterList p("BC Dirichlet Thermal Contact");
  p.set("Prefix", "Target_");
  p.set("Field Library", pb.getFieldLibraryBase());
  p.set("Names", names);
  p.set("Temperature", value);
  p.set("Scaling Parameters", scaleParams);

  RCP< PHX::Evaluator<panzer::Traits> > op =
    rcp(new charon::BC_ThermalContact<EvalT,panzer::Traits>(p));

  fm.template registerEvaluator<EvalT>(op);
}

#endif
