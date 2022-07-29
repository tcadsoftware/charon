
#ifndef CHARON_BCSTRATEGY_DIRICHLET_TRAPEZOID_IMPL_HPP
#define CHARON_BCSTRATEGY_DIRICHLET_TRAPEZOID_IMPL_HPP

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
#include "Panzer_GlobalData.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

// Evaluators
#include "Charon_BC_Trapezoid.hpp"
#include "Charon_DDLatticeBC_Trapezoid.hpp"


// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_Dirichlet_Trapezoid<EvalT>::
BCStrategy_Dirichlet_Trapezoid(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "Trapezoid");

  // default
  isLatTDof = false;
  isIonDof = false;
  isFermiPin = false;
  ionDens = 0.0; 
}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_Trapezoid<EvalT>::
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

  // check if Fermi Level Pinning is turned on
  if ( dataPList->isParameter("Fermi Level Pinning") )
    isFermiPin = dataPList->template get<bool>("Fermi Level Pinning");

  if (isFermiPin) // Ion density is pinned at a given value on a contact
    ionDens = dataPList->template get<double>("Contact Ion Density");

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

    // for the DD+Lattice and DD+Lattice+Ion formulations, need to gather T (lattice temperature), since (phi,n,p) depend on T
    if (dof_it->first == n.dof.latt_temp)
    {
      isLatTDof = true;
      this->required_dof_names.push_back(dof_it->first);
    }

    // for the DD+Ion and DD+Lattice+Ion formulations, need to gather ION_DENSITY
    if (dof_it->first == n.dof.iondensity)
    {
      isIonDof = true;
      this->required_dof_names.push_back(dof_it->first);
    }

  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
                             "Error: the name \"" << this->m_bc.equationSetName()
                             << "\" is not a valid DOF for the boundary condition:\n"
                             << this->m_bc << "\n");

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_Dirichlet_Trapezoid<EvalT>::
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

  // get the element block id and physics id for the given physics block
  string ebID_pb = pb.elementBlockID();
  string pbID = pb.physicsBlockID();

  // get the element block ID for the given boundary condition
  string ebID_bc = this->m_bc.elementBlockID();

  if (ebID_pb != ebID_bc)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: " << pbID << " corresponds to "
      << ebID_pb << ", while the BC corresponds to " << ebID_bc << "! \n");

  // flag to determine if Fermi Dirac is turned on
  bool bUseFD = false;  // default

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();

  // allow only one equation set per physics block
  if (pbParamList->numParams() > 1)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
       "The physics block " << pbParamList->name() << " has more than one equation sets ! ");

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  const ParameterList& options = eqSetPList.sublist("Options");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  // determine if Fermi Dirac is turned on
  if (options.isParameter("Fermi Dirac"))
  {
    string fermiDirac = options.get<string>("Fermi Dirac");
    if (fermiDirac == "True")  bUseFD = true;
  }

  ParameterList incmpl_ioniz;
  incmpl_ioniz.sublist("Acceptor"); incmpl_ioniz.sublist("Donor");
  // determine if Acceptor Incomplete Ionization is turned on
  if(options.isParameter("Acceptor Incomplete Ionization")) {
    string AccIncmplIoniz = options.get<string>("Acceptor Incomplete Ionization");
    if(AccIncmplIoniz == "On") {
      assert(eqSetPList.isParameter("Model ID"));
      string model_name = eqSetPList.get<string>("Model ID");
      assert(models.isSublist(model_name));
      const ParameterList& model_param = models.sublist(model_name);
      if(model_param.isSublist("Incomplete Ionized Acceptor")) {
        incmpl_ioniz.sublist("Acceptor") =
          model_param.sublist("Incomplete Ionized Acceptor").sublist("Model");
      }
    }
  }

  // determine if Donor Incomplete Ionization is turned on
  if(options.isParameter("Donor Incomplete Ionization")) {
    string DonIncmplIoniz = options.get<string>("Donor Incomplete Ionization");
    if(DonIncmplIoniz == "On") {
      assert(eqSetPList.isParameter("Model ID"));
      string model_name = eqSetPList.get<string>("Model ID");
      assert(models.isSublist(model_name));
      const ParameterList& model_param = models.sublist(model_name);
      if(model_param.isSublist("Incomplete Ionized Donor")) {
        incmpl_ioniz.sublist("Donor") =
          model_param.sublist("Incomplete Ionized Donor").sublist("Model");
      }
    }
  }

  // determine if Ion Charge is given
  int ionCharge = 1;
  if (options.isParameter("Ion Charge"))
    ionCharge = options.get<int>("Ion Charge");

  RCP<const charon::Names> names = rcp(new charon::Names(1,prefix,discfields,discsuffix));

  // get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // get the Data parameter list
  RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
  TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));

  ParameterList p("BC Dirichlet Trapezoidal");
  p.set("Prefix", "Target_");
  p.set("Field Library", pb.getFieldLibraryBase());
  p.set("Names", names);
  //p.set("DC Offset", this->m_bc.params()->template get<double>("DC Offset"));
  p.set("Amplitude", this->m_bc.params()->template get<double>("Amplitude"));
  p.set("Period", this->m_bc.params()->template get<double>("Period"));
  p.set("Rise Time", this->m_bc.params()->template get<double>("Rise Time"));
  p.set("Fall Time", this->m_bc.params()->template get<double>("Fall Time"));
  //p.set("Delay", this->m_bc.params()->template get<double>("Delay"));
  //p.set("Duty Cycle", this->m_bc.params()->template get<double>("Duty Cycle"));
  p.set("Sideset ID",this->m_bc.sidesetID());
 
  //Insert the paramaeter library
  p.set<RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);

  // It should be optional to specify these values
  if (dataPList->isParameter("DC Offset"))
    p.set("DC Offset", this->m_bc.params()->template get<double>("DC Offset"));
  else
    p.set("DC Offset", 0.0);
  if (dataPList->isParameter("Delay"))
    p.set("Delay", this->m_bc.params()->template get<double>("Delay"));
  else
    p.set("Delay", 0.0);
  if (dataPList->isParameter("Duty Cycle"))  
    p.set("Duty Cycle", this->m_bc.params()->template get<double>("Duty Cycle"));
  else
    p.set("Duty Cycle", 1.0);
  if (dataPList->isParameter("Number Pulses"))  
    p.set("Number Pulses", this->m_bc.params()->template get<int>("Number Pulses"));
  else
    p.set("Number Pulses", 1);

  p.set<bool>("Fermi Dirac", bUseFD);
  p.set("Scaling Parameters", scaleParams);
  p.sublist("Incomplete Ionization") = incmpl_ioniz;

  // build DDLatticeBC_Trapezoid for DD+Lattice, DD+Ion, and DD+Lattice+Ion formulations
  if (isLatTDof || isIonDof)
  {
    p.set<bool>("Solve Ion", isIonDof);
    p.set<int>("Ion Charge", ionCharge);
    p.set<bool>("Fermi Level Pinning", isFermiPin);
    p.set<double>("Contact Ion Density", ionDens); 

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::DDLatticeBC_Trapezoid<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // build the BC_Trapezoid evaluator for isothermal DD formulations
  else
  {
    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::BC_Trapezoid<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

#endif
