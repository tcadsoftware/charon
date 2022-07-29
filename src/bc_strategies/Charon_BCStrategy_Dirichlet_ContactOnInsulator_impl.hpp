
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

  Teuchos::RCP<const Teuchos::ParameterList> bc_params = bc.params();

  this->m_names = (input_pl->isParameter("Names") ? input_pl->get<Teuchos::RCP<charon::Names> >("Names") : Teuchos::rcp(new charon::Names(1,"","","")));
  this->basis = (input_pl->isParameter("Names") ? input_pl->get<Teuchos::RCP<panzer::PureBasis> >("Basis") : Teuchos::RCP<panzer::PureBasis>());
  this->small_signal_perturbation = (bc_params->isParameter("Small Signal Perturbation") ? bc_params->get<double>("Small Signal Perturbation") : 0.0 );
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp; 
  using std::vector;
  using std::string;
  using std::pair;

  // get the Data parameter list
  RCP<const ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  // validate the Data parameter list
  RCP<ParameterList> valid_params = this->getValidParameters();
  dataPList->validateParameters(*valid_params);

  // linear ramp voltage parameter list
  bLinRamp = false; 
  if (dataPList->isSublist("Linear Ramp"))
  {
    bLinRamp = true; 
    const ParameterList& tmpPL = dataPList->sublist("Linear Ramp"); 
    linPList = rcp(new ParameterList(tmpPL));
  }

  // trapezoid pulse voltage parameter list
  bTrapezoid = false; 
  if (dataPList->isSublist("Trapezoid Pulse"))
  {
    bTrapezoid = true; 
    const ParameterList& tmpPL = dataPList->sublist("Trapezoid Pulse"); 
    trapzPList = rcp(new ParameterList(tmpPL));
  }

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
    p.set("Sideset ID",this->m_bc.sidesetID());
    p.set("Frequency Domain", isFreqDom);
    
    //Loca stuff
    p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);

    if(this->m_bc.params()->template isType<std::string>("Voltage"))
      p.set("Voltage", this->m_bc.params()->template get<std::string>("Voltage"));

    else if(this->m_bc.params()->template isType<double>("Voltage"))
    {
      p.set("Voltage", this->m_bc.params()->template get<double>("Voltage"));
    }

    else if (this->m_bc.params()->isParameter("Varying Voltage"))
    {
      p.setEntry("Varying Voltage", this->m_bc.params()->getEntry("Varying Voltage"));
      if (this->m_bc.params()->isParameter("Initial Voltage"))
        p.set("Initial Voltage", this->m_bc.params()->template get<double>("Initial Voltage"));
    }

    else if (bLinRamp)  // time-dependent linear ramp voltage source
    {
      p.set<bool>("Enable Linear Ramp", bLinRamp); 
      p.set<RCP<ParameterList>>("Linear Ramp ParameterList", linPList); 
    }

    else if (bTrapezoid)  // time-dependent linear ramp voltage source
    {
      p.set<bool>("Enable Trapezoid Pulse", bTrapezoid); 
      p.set<RCP<ParameterList>>("Trapezoid Pulse ParameterList", trapzPList); 
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error!\
        Valid parameters/lists are Voltage, Varying Voltage, Linear Ramp or Trapezoid Pulse!" << std::endl);
       
  if(this->isFreqDom)
  {
    // if using a small signal analysis, pass the 
    p.set<double>("Small Signal Perturbation", this->small_signal_perturbation);
  }

    p.set("Work Function", this->m_bc.params()->template get<double>("Work Function"));

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::BC_ContactOnInsulator<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<Teuchos::ParameterList>
charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::string;

  RCP<ParameterList> p = Teuchos::rcp(new ParameterList);
  p->set<double>("Work Function", 0.0, "Metal work function in (eV)");
  p->set<double>("Voltage", 0.0, "Apply a single DC voltage in (V)");
  p->set<string>("Varying Voltage", "", "Apply sweeping voltages in (V)");
  p->set<double>("Initial Voltage", 0.0, "Initial voltage for a voltage sweep in (V)"); 

  ParameterList& linPL = p->sublist("Linear Ramp", false, "Sublist defining Linear Ramp voltage source");
  linPL.set<double>("Initial Time", 0.0, "Initial time in (s)");
  linPL.set<double>("Final Time", 0.0, "Final time in (s)");
  linPL.set<double>("Initial Voltage", 0.0, "Initial voltage in (V)");
  linPL.set<double>("Final Voltage", 0.0, "Final voltage in (V)"); 

  ParameterList& trapzPL = p->sublist("Trapezoid Pulse", false, "Sublist defining Trapezoid Pulse voltage source");
  trapzPL.set<double>("DC Offset",0.0);
  trapzPL.set<double>("Amplitude",0.0);
  trapzPL.set<double>("Period",0.0);
  trapzPL.set<double>("Rise Time",0.0);
  trapzPL.set<double>("Fall Time",0.0);
  trapzPL.set<double>("Delay",0.0);
  trapzPL.set<double>("Duty Cycle",1.0);
  trapzPL.set<int>("Number Pulses",1);

  return p;
}

#endif
