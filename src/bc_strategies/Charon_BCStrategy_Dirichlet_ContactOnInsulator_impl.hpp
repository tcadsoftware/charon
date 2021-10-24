
#ifndef CHARON_BCSTRATEGY_DIRICHLET_CONTACTONINSULATOR_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_CONTACTONINSULATOR_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

// Evaluators
#include "Charon_BC_ContactOnInsulator.hpp"

// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>::
BCStrategy_Dirichlet_ContactOnInsulator(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data,
                                        Teuchos::RCP<Teuchos::ParameterList> input_pl) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Contact On Insulator");

  this->m_names = (input_pl->isParameter("Names") ? input_pl->get<Teuchos::RCP<charon::Names> >("Names") : Teuchos::rcp(new charon::Names(1,"","","")));
  this->basis = (input_pl->isParameter("Names") ? input_pl->get<Teuchos::RCP<panzer::PureBasis> >("Basis") : Teuchos::RCP<panzer::PureBasis>());
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  // this setup method is hit for time domain simulations, and avoided for frequency domain simulations
  this->isFreqDom = false;

  // set dirichlet conditions on all dofs (potential and carrier density)
  const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
   dofs.begin(); dof_it != dofs.end(); ++dof_it)
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

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
                             "Error: the name \"" << this->m_bc.equationSetName()
                             << "\" is not a valid DOF for the boundary condition:\n"
                             << this->m_bc << "\n");

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                           const panzer::PhysicsBlock& pb,
                           const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                           const Teuchos::ParameterList& models,
                           const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // build and register all closure models
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  std::string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<std::string>("Prefix") : "";
  std::string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<std::string>("Discontinuous Fields") : "";
  std::string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<std::string>("Discontinuous Suffix") : "";

  RCP<const charon::Names> names = (isFreqDom ? rcp(new charon::Names(1,prefix,discfields,discsuffix,this->m_names->FDsuffix())) : 
                                                rcp(new charon::Names(1,prefix,discfields,discsuffix)) );

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // build the BC_ContactOnInsulator evaluator
  {
    ParameterList p("BC Dirichlet Contact On Insulator");
    p.set("Prefix", "Target_");
    // p.set("Basis", basis);
    p.set("Field Library", pb.getFieldLibraryBase());
    p.set("Names", names);
    p.set("Scaling Parameters", scaleParams);
    p.set("Frequency Domain", isFreqDom);

    //Loca stuff
    p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);

    if(this->m_bc.params()->template isType<std::string>("Voltage"))
      p.set("Voltage", this->m_bc.params()->template get<std::string>("Voltage"));
    else if (this->m_bc.params()->isParameter("Varying Voltage"))
      p.setEntry("Varying Voltage",
                 this->m_bc.params()->getEntry("Varying Voltage"));
    else
      p.set("Voltage", this->m_bc.params()->template get<double>("Voltage"));
    p.set("Work Function", this->m_bc.params()->template get<double>("Work Function"));

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::BC_ContactOnInsulator<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

#endif
