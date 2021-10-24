

#ifndef CHARON_BC_STRATEGY_NEUMANN_THERMALCONTACT_IMPL_HPP
#define CHARON_BC_STRATEGY_NEUMANN_THERMALCONTACT_IMPL_HPP

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
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Constant.hpp"

#include "Charon_Names.hpp"
#include "Charon_NeumannBC_ThermalContact.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Neumann_ThermalContact<EvalT>::
BCStrategy_Neumann_ThermalContact(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Neumann_DefaultImpl<EvalT>(bc,global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Neumann Thermal Contact" );
}

// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_ThermalContact<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::string;

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

  // Power is given in [W/cm^2]
  if ( dataPList->isParameter("Power") )
  {
    paramName.push_back("Power");
    value = dataPList->template get<double>("Power");
    if ( dataPList->isParameter("Temperature") )
      temp = dataPList->template get<double>("Temperature");
    else
      temp = 300.0;
  }

  // Surface Resistance is given in [K.cm^2/W]
  if ( dataPList->isParameter("Surface Resistance") )
  {
    paramName.push_back("Surface Resistance");
    value = dataPList->template get<double>("Surface Resistance");
    temp = dataPList->template get<double>("Temperature");  // Temperature must be given in [K]
  }

  // Surface Conductance is given in [W/(K.cm^2)]
  if ( dataPList->isParameter("Surface Conductance") )
  {
    paramName.push_back("Surface Conductance");
    value = dataPList->template get<double>("Surface Conductance");
    temp = dataPList->template get<double>("Temperature");  // Temperature must be given in [K]
  }

  // allow only one of {Power, Surface Resistance, Surface Conductance} is given
  if ( (paramName.size() > 1) || (paramName.size() < 1) )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: You must and can only provide either Power, or Surface Resistance, "
        << "or Surface Conductance for the Neumann Thermal Contact !" );

  // get the name of lattice temperature
  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));
  string dof_name = names->dof.latt_temp;

  // check if Equation Set Name = ALL_DOFS or dof.latt_temp
  string eqsn = this->m_bc.equationSetName();
  if ( (eqsn != "ALL_DOFS") && (eqsn != dof_name) )
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
        << "Error: Equation Set Name must be either ALL_DOFS or "
        << dof_name << ". But you entered \"" << eqsn << "\" ! \n");

  const string residual_name = "Residual_" + dof_name;
  const string flux_name = "Heat_Flux";

  // For now, assume that only one IR is used for everything.
  const std::map<int,RCP< panzer::IntegrationRule > >& ir = side_pb.getIntegrationRules();
  TEUCHOS_ASSERT(ir.size() == 1);  // Doesn't support mixed ir yet!

  const int integration_order = ir.begin()->second->order();

  this->addResidualContribution(residual_name,dof_name,flux_name,integration_order,side_pb);

}

// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_ThermalContact<EvalT>::
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

  const std::vector<std::tuple<string,string,string,int,Teuchos::RCP<panzer::PureBasis>,
        Teuchos::RCP<panzer::IntegrationRule> > > data = this->getResidualContributionData();

  string residual_name = std::get<0>(data[0]);
  string dof_name = std::get<1>(data[0]);
  string flux_name = std::get<2>(data[0]);

  RCP<panzer::IntegrationRule> ir = std::get<5>(data[0]);
  RCP<const panzer::FieldLayoutLibrary> fll = pb.getFieldLibrary()->buildFieldLayoutLibrary(*ir);
  RCP<panzer::BasisIRLayout> basis = fll->lookupLayout(dof_name);

  // build and register all closure models
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);
  pb.buildAndRegisterDOFProjectionsToIPEvaluators(fm,Teuchos::null,user_data);

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // obtain the heat flux at thermal Neumann BC
  {
    ParameterList p("NeumannBC Heat Flux");
    p.set("Names", names);
    p.set("Data Layout", ir->dl_scalar);
    p.set("Flux Name", flux_name);
    p.set("DOF Name", dof_name);
    p.set("Parameter Name", paramName[0]);
    p.set("Value", value);
    p.set("Temperature", temp);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::NeumannBC_ThermalContact<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // add contribution to the residual
  {
    using panzer::EvaluatorStyle;
    using panzer::Integrator_BasisTimesScalar;
    using panzer::Traits;
    using PHX::Evaluator;
    RCP<Evaluator<Traits>> op = rcp(new 
      Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
      residual_name, flux_name, *basis, *ir, -1));
    fm.template registerEvaluator<EvalT>(op);
  }
}

// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_ThermalContact<EvalT>::
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
void charon::BCStrategy_Neumann_ThermalContact<EvalT>::
postRegistrationSetup(typename panzer::Traits::SetupData /* d */,
                      PHX::FieldManager<panzer::Traits>& /* vm */)
{

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Neumann_ThermalContact<EvalT>::
evaluateFields(typename panzer::Traits::EvalData /* d */)
{

}

#endif
