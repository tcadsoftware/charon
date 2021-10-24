
#ifndef   Charon_BCStrategy_Dirichlet_CurrentConstraint_impl_hpp
#define   Charon_BCStrategy_Dirichlet_CurrentConstraint_impl_hpp

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Charon
#include "Charon_BC_CurrentConstraint.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

// Panzer
#include "Panzer_GlobalData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_ScalarParameterEntry.hpp"

// Phalanx
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

// Teuchos
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Default Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
charon::BCStrategy_Dirichlet_CurrentConstraint<EvalT>::
BCStrategy_Dirichlet_CurrentConstraint(
  const panzer::BC&                       bc,
  const Teuchos::RCP<panzer::GlobalData>& globalData)
  :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, globalData)
{
  using panzer::accessScalarParameter;
  using std::cout;
  using std::endl;
  using std::string;
  string strategy(this->m_bc.strategy());
  TEUCHOS_ASSERT((strategy == "Constant Current")  or
                 (strategy == "Resistor Contact"))
  string ssID(this->m_bc.sidesetID()), controlName("");
  if (strategy == "Constant Current")
    controlName = ssID + "ConstantCurrentVoltage";
  else // if (strategy == "Resistor Contact")
    controlName = ssID + "ResistorContactVoltage";
  cout << "ssID = "         << ssID        << ", "
       << "control name = " << controlName << endl;
  voltageParameter_ =
    accessScalarParameter<EvalT>(controlName, *globalData->pl);
} // end of Default Constructor

///////////////////////////////////////////////////////////////////////////////
//
//  setup()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void
charon::BCStrategy_Dirichlet_CurrentConstraint<EvalT>::
setup(
  const panzer::PhysicsBlock&   sidePB,
  const Teuchos::ParameterList& /* userData */)
{
  using charon::Names;
  using std::runtime_error;
  using std::string;
  using Teuchos::is_null;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef std::vector<std::pair<std::string, Teuchos::RCP<panzer::PureBasis>>>
    DofVec;
  typedef
    std::vector<std::pair<std::string, Teuchos::RCP<panzer::PureBasis>>>::
    const_iterator DofVecIter;

  // get the physics block parameter list
  RCP<const ParameterList> pbPL = sidePB.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPL = pbPL->sublist("child0");

  // get any prefix or suffix parameters
  string prefix(""), discFields(""), discSuffix("");
  if (eqSetPL.isParameter("Prefix"))
    prefix = eqSetPL.get<string>("Prefix");
  if (eqSetPL.isParameter("Discontinuous Fields"))
    discFields = eqSetPL.get<string>("Discontinuous Fields");
  if (eqSetPL.isParameter("Discontinuous Suffix"))
    discSuffix = eqSetPL.get<string>("Discontinuous Suffix");

  // get charon::Names
  RCP<const Names> names = rcp(new Names(3, prefix, discFields, discSuffix));

  string carrierDofName("");
  bjt1DBaseContact_ = this->m_bc.params()->isParameter("BJT1D Base Doping Type");
  if (bjt1DBaseContact_)
  {
    string baseDopingType(
      this->m_bc.params()->template get<string>("BJT1D Base Doping Type"));
    if (baseDopingType == "Acceptor")
      carrierDofName = names->dof.hdensity;
    else if (baseDopingType == "Donor")
      carrierDofName = names->dof.edensity;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, runtime_error, "Error:  Invalid \""    \
        "Base Doping Type\" (\"" << baseDopingType << "\") for BJT 1-D Base " \
        "Contact; must be either \"Acceptor\" or \"Donor\".")
  } // end if (bjt1DBaseContact)

  // Set the dirichlet conditions on all of the degrees of freedom (potential
  // and carrier density).
  const DofVec& dofs = sidePB.getProvidedDOFs();
  for (DofVecIter dofIt = dofs.begin(); dofIt != dofs.end(); ++dofIt)
  {
    // If we're not specifying a BJT1D Base Contact, then we want to do the
    // following on all the degrees of freedom.  In the BJT1D Base Contact case,
    // we want to do the following only for the electric potential and either
    // the electron or hole density, depeding on the "BJT1D Base Doping Type".
    // When solving the temperature equation, we need to gather "Lattice Temperature"
    // here, not to set a value for the temperature.

    if ((not bjt1DBaseContact_         )    or
        (dofIt->first == names->dof.phi)    or
        (dofIt->first == carrierDofName)    or
        (dofIt->first == names->dof.latt_temp))
    {
      // for the DDLattice and DDIonLattice formulations, need to gather
      // T (Lattice Temperature), since (phi,n,p) depend on T
      if (dofIt->first == names->dof.latt_temp)
        this->required_dof_names.push_back(dofIt->first);

      else  // do the following for potential and carrier densities
      {
        // We need the degree of freedom values to form the residual.
        this->required_dof_names.push_back(dofIt->first);

        // Create the unique residual name.
        string residualName("Residual_" + dofIt->first);

        // Map the residual to the degree of freedom and to the target field.
        this->residual_to_dof_names_map[residualName]    = dofIt->first;
        this->residual_to_target_field_map[residualName] = "Target_" + dofIt->first;

        // For now assume all dofs use the same basis.
        this->basis_ = dofIt->second;
      }
    }  // end if we're doing a BJT 1-D Base Contact case
  } // end loop over the degrees of freedom

  TEUCHOS_TEST_FOR_EXCEPTION(is_null(this->basis_), runtime_error, "Error:  " \
    "The name \"" << this->m_bc.equationSetName() << "\" is not a valid DOF " \
    "for the boundary condition:\n" << this->m_bc << "\n");

} // end of setup()

///////////////////////////////////////////////////////////////////////////////
//
//  buildAndRegisterEvaluators()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void
charon::BCStrategy_Dirichlet_CurrentConstraint<EvalT>::
buildAndRegisterEvaluators(
  PHX::FieldManager<panzer::Traits>&                                 fm,
  const panzer::PhysicsBlock&                                        pb,
  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
  const Teuchos::ParameterList&                                      models,
  const Teuchos::ParameterList&                                      userData)
  const
{
  using charon::BC_CurrentConstraint;
  using charon::Names;
  using charon::Scaling_Parameters;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::invalid_argument;
  using std::runtime_error;
  using std::string;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Build and register all of the closure models.
  pb.buildAndRegisterClosureModelEvaluators(fm, factory, models, userData);

  // Get the element block ID and physics ID for the given physics block, along
  // with the element block ID for the given boundary condition, and ensure
  // that the element block IDs match.
  string ebIdPB(pb.elementBlockID()), pbId(pb.physicsBlockID()),
    ebIdBC(this->m_bc.elementBlockID());
  TEUCHOS_TEST_FOR_EXCEPTION(ebIdPB != ebIdBC, runtime_error, "Error:  "
    << pbId << " corresponds to " << ebIdPB << ", while the BC corresponds to "
    << ebIdBC << "!\n");

  // Get the physics block parameter list, and ensure there's only one equation
  // set per physics block.
  RCP<const ParameterList> pbParamList = pb.getParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(pbParamList->numParams() > 1, runtime_error,
    "The physics block " << pbParamList->name() << " has more than one "    \
    "equation set!");

  // Get the equation set parameter list, along with any prefix or suffix
  // parameters.
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  string prefix(""), discFields(""), discSuffix("");
  if (eqSetPList.isParameter("Prefix"              ))
    prefix     = eqSetPList.get<string>("Prefix"              );
  if (eqSetPList.isParameter("Discontinuous Fields"))
    discFields = eqSetPList.get<string>("Discontinuous Fields");
  if (eqSetPList.isParameter("Discontinuous Suffix"))
    discSuffix = eqSetPList.get<string>("Discontinuous Suffix");

  // Determine if Fermi-Dirac is turned on.
  bool bUseFD(false);
  const ParameterList& options = eqSetPList.sublist("Options");
  if (options.isParameter("Fermi Dirac"))
  {
    string fermiDirac(options.get<string>("Fermi Dirac"));
    if (fermiDirac == "True")
      bUseFD = true;
  } // end if (options.isParameter("Fermi Dirac"))

  // Determine if bUseRefE = false
  const string eqnSetType = eqSetPList.get<string>("Type");
  bool bUseRefE = true;  // default

  // Comment out the following two lines, because we want to always use the Reference Energy for
  // BC calculation after the Poisson equation reformulation for the equation sets listed below.
  // if ((eqnSetType == "DDLattice") || (eqnSetType == "DDIon") || (eqnSetType == "DDIonLattice") )
  //  bUseRefE = false;

  // Determine if either Donor or Acceptor Incomplete Ionization is turned on.
  bool donIncIonOn(false), accIncIonOn(false);
  if ((options.isParameter("Donor Incomplete Ionization")           )  and
      (options.get<string>("Donor Incomplete Ionization")    == "On"))
    donIncIonOn = true;
  if ((options.isParameter("Acceptor Incomplete Ionization")        )  and
      (options.get<string>("Acceptor Incomplete Ionization") == "On"))
    accIncIonOn = true;
  ParameterList modelParam;
  if ((donIncIonOn) or (accIncIonOn))
  {
    TEUCHOS_TEST_FOR_EXCEPTION(not eqSetPList.isParameter("Model ID"),
      invalid_argument, "Error:  Donor or Acceptor Incomplete Ionization is " \
      "On, but there is no \"Model ID\".");
    string modelName(eqSetPList.get<string>("Model ID"));
    TEUCHOS_TEST_FOR_EXCEPTION(not models.isSublist(modelName),
      invalid_argument, "Error:  Donor or Acceptor Incomplete Ionization is " \
      "On with \"Model ID\" = \"" + modelName + "\", but there is no \"" +
      modelName + "\" sublist.");
    modelParam = models.sublist(modelName);
  }
  ParameterList incompleteIonizationPL;
  incompleteIonizationPL.sublist("Donor");
  incompleteIonizationPL.sublist("Acceptor");
  if (donIncIonOn and (modelParam.isSublist("Incomplete Ionized Donor"   )))
    incompleteIonizationPL.sublist("Donor"   ) =
      modelParam.sublist("Incomplete Ionized Donor"   ).sublist("Model");
  if (accIncIonOn and (modelParam.isSublist("Incomplete Ionized Acceptor")))
    incompleteIonizationPL.sublist("Acceptor") =
      modelParam.sublist("Incomplete Ionized Acceptor").sublist("Model");

  // Build the BC_CurrentConstraint evaluator.
  RCP<const Names> names = rcp(new Names(1, prefix, discFields, discSuffix));
  RCP<Scaling_Parameters> scaleParams =
    userData.get<RCP<Scaling_Parameters>>("Scaling Parameter Object");
  ParameterList p("BC Dirichlet Current Constraint");

  // Get the data holder for the empirical damage model
  auto damage_data =
    userData.get<Teuchos::RCP<charon::EmpiricalDamage_Data> >("empirical damage data");
  p.set("empirical damage data", damage_data);

  p.set("Prefix",             "Target_"               );
  p.set("Field Library",      pb.getFieldLibraryBase());
  p.set("Names",              names                   );
  p.set("Voltage Control",    voltageParameter_       );
  p.set("Scaling Parameters", scaleParams             );
  p.set<bool>("Fermi Dirac",        bUseFD           );
  p.set<bool>("BJT1D Base Contact", bjt1DBaseContact_);
  p.set<bool>("Use Reference Energy", bUseRefE);
  p.sublist("Incomplete Ionization") = incompleteIonizationPL;
  RCP<Evaluator<Traits>> op = rcp(new BC_CurrentConstraint<EvalT, Traits>(p));
  fm.template registerEvaluator<EvalT>(op);
} // end of buildAndRegisterEvaluators()

#endif // Charon_BCStrategy_Dirichlet_CurrentConstraint_impl_hpp
