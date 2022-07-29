
#ifndef   CHARON_CLOSURE_MODEL_FACTORY_IMPL_HPP
#define   CHARON_CLOSURE_MODEL_FACTORY_IMPL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// Boost
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

// Charon
// Non-Evaluator files
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_EmpiricalDamage_Data.hpp"

// Scaling parameters and constant temperature
#include "Charon_LatticeTemp_Constant.hpp"

// Doping profiles
#include "Charon_Doping_StepJunction.hpp"
#include "Charon_Doping_Function.hpp"
#include "Charon_DopingRaw_Function.hpp"
#include "Charon_Space_Charge.hpp"

// Bandgap and other energy-related evaluators
#include "Charon_Reference_Energy.hpp"
#include "Charon_BandGap_TempDep.hpp"
#include "Charon_Intrinsic_FermiEnergy.hpp"
#include "Charon_CondVale_Band.hpp"
#include "Charon_BandGap_Nitride.hpp"
#include "Charon_Permittivity_Nitride.hpp"
#include "Charon_Relative_Permittivity.hpp"

// Intrinsic conc models
#include "Charon_IntrinsicConc_Default.hpp"
#include "Charon_IntrinsicConc_OldSlotboom.hpp"
#include "Charon_IntrinsicConc_Persson.hpp"
#include "Charon_IntrinsicConc_Slotboom.hpp"
#include "Charon_IntrinsicConc_Harmon.hpp"
#include "Charon_EffectiveDOS_Simple.hpp"
#include "Charon_EffectiveDOS_Nitride.hpp"

// Mobility models
#include "Charon_Mobility_Default.hpp"
#include "Charon_Mobility_Analytic.hpp"
#include "Charon_Mobility_Arora.hpp"
#include "Charon_Mobility_UniBo.hpp"
#include "Charon_Mobility_Masetti.hpp"
#include "Charon_Mobility_PhilipsThomas.hpp"
#include "Charon_Mobility_Lucent.hpp"
#include "Charon_Mobility_RigidPointIon.hpp"
#include "Charon_Mobility_GaAs.hpp"
#include "Charon_Mobility_Albrecht.hpp"
#include "Charon_Mobility_Farahmand.hpp"
#include "Charon_Mobility_DopantTempDep.hpp"
#include "Charon_DiffCoeff_Default.hpp"
#include "Charon_DiffCoeff_IonDep.hpp"
#include "Charon_Mobility_MOSFET.hpp"
#include "Charon_Mobility_Shirahata.hpp"

// Recombination models
#include "Charon_SRHLifetime_Constant.hpp"
#include "Charon_SRHLifetime_Function.hpp"
#include "Charon_RecombRate_SRH.hpp"
#include "Charon_RecombRate_TrapSRH.hpp"
#include "Charon_RecombRate_DynamicTraps.hpp"
#include "Charon_QuasiFermiLevels.hpp"
#include "Charon_RecombRate_Defect_Cluster.hpp"
#include "Charon_RecombRate_Empirical_Defect.hpp"
#include "Charon_RecombRate_Radiative.hpp"
#include "Charon_RecombRate_Auger.hpp"
#include "Charon_RecombRate_Total.hpp"

// Generation models
#include "Charon_Avalanche_vanOverstraeten.hpp"
#include "Charon_Avalanche_Okuto.hpp"
#include "Charon_Avalanche_Lackner.hpp"
#include "Charon_Avalanche_UniBo.hpp"
#include "Charon_Avalanche_UniBoNew.hpp"
#include "Charon_Avalanche_Selberherr.hpp"
#include "Charon_Avalanche_CrowellSze.hpp"
#include "Charon_Ionization_ParticleStrike.hpp"
#include "Charon_KimptonTID.hpp"
#include "Charon_Band2Band_Tunneling_Local.hpp"

// Initial conditions
#include "Charon_IC_Equilibrium_Density.hpp"
#include "Charon_IC_Equilibrium_Potential.hpp"
#include "Charon_IC_Gauss.hpp"
#include "Charon_IC_Function.hpp"
#include "Charon_IC_Remap.hpp"

// Allows the initial conditions read from an Exodus file to be scaled
// prior to use as an initial guess
#include "Charon_GatherScaledFields_decl.hpp"

// Additional evaluators for the Fermi-Dirac statistics
#include "Charon_Degeneracy_Factor.hpp"

//Electric field
#include "Charon_FEM_ElectricField.hpp"

// Heat capacity and thermal conductivity models
#include "Charon_HeatCapacity_TempDep.hpp"
#include "Charon_HeatCapacity_PowerLawTempDep.hpp"
#include "Charon_ThermalConduct_TempDep.hpp"
#include "Charon_ThermalConduct_PowerLawTempDep.hpp"
#include "Charon_ThermalConduct_LinearTempDep.hpp"
#include "Charon_ThermalConduct_LinearIonDep.hpp"
#include "Charon_Analytic_HeatGeneration.hpp"

// Comparison of computed to analytic results
#include "Charon_AnalyticComparison.hpp"
#include "Charon_MMS_AnalyticSolutions.hpp"
#include "Charon_NormCalculation.hpp" //NEW ADDITION
#include "Charon_ManufacturedSolution.hpp" //NEW ADDITION

// Miscellaneous
#include "Charon_FEM_GradNegPotential.hpp"
#include "Charon_FEM_Velocity.hpp"
#include "Charon_SoretCoeff_TempDep.hpp"
#include "Charon_OptGen_Function.hpp"
#include "Charon_MoleFraction_Function.hpp"
#include "Charon_BulkFixCharge_Function.hpp"
#include "Charon_ThermodiffCoeff_Default.hpp"
#include "Charon_ThermodiffCoeff_Custom.hpp"
#include "Charon_Vector.hpp"
#include "Charon_Initial_PotentialGrad.hpp"

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_GlobalStatistics.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_STK_GatherFields.hpp"

// Teuchos
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TypeNameTraits.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Using statements used in the various routines below.
//
///////////////////////////////////////////////////////////////////////////////
#define CHARON_USINGS_ETC(input) CHARON_USINGS_ETC_##input
#define CHARON_USINGS_ETC_DEFINE_NONE                                         \
  using charon::Names;                                                        \
  using panzer::Traits;                                                       \
  using PHX::Evaluator;                                                       \
  using std::size_t;                                                          \
  using std::string;                                                          \
  using std::stringstream;                                                    \
  using std::vector;                                                          \
  using Teuchos::ParameterList;                                               \
  using Teuchos::RCP;                                                         \
  using Teuchos::rcp;
#define CHARON_USINGS_ETC_DEFINE_BIRL                                         \
  CHARON_USINGS_ETC_DEFINE_NONE                                               \
  typedef panzer::BasisIRLayout BIRL;
#define CHARON_USINGS_ETC_DEFINE_IR                                           \
  CHARON_USINGS_ETC_DEFINE_NONE                                               \
  typedef panzer::IntegrationRule IR;
#define CHARON_USINGS_ETC_DEFINE_BIRL_IR                                      \
  CHARON_USINGS_ETC_DEFINE_BIRL                                               \
  typedef panzer::IntegrationRule IR;
#define CHARON_USINGS_ETC_DEFINE_NAMES                                        \
  CHARON_USINGS_ETC_DEFINE_NONE                                               \
  const RCP<const Names>& names = defaults.get<RCP<const Names>>("Names");    \
  const RCP<const Names>& inp_names =                                         \
    Teuchos::rcp(new const charon::Names(1,names->prefix(),                   \
                 names->discfields(),names->discsuffix()));
#define CHARON_USINGS_ETC_DEFINE_BIRL_NAMES                                   \
  CHARON_USINGS_ETC_DEFINE_BIRL                                               \
  const RCP<const Names>& names = defaults.get<RCP<const Names>>("Names");    \
  const RCP<const Names>& inp_names =                                         \
    Teuchos::rcp(new const charon::Names(1,names->prefix(),                   \
                 names->discfields(),names->discsuffix()));
#define CHARON_USINGS_ETC_DEFINE_IR_NAMES                                     \
  CHARON_USINGS_ETC_DEFINE_IR                                                 \
  const RCP<const Names>& names = defaults.get<RCP<const Names>>("Names");    \
  const RCP<const Names>& inp_names =                                         \
    Teuchos::rcp(new const charon::Names(1,names->prefix(),                   \
                 names->discfields(),names->discsuffix()));
#define CHARON_USINGS_ETC_DEFINE_BIRL_IR_NAMES                                \
  CHARON_USINGS_ETC_DEFINE_BIRL_IR                                            \
  const RCP<const Names>& names = defaults.get<RCP<const Names>>("Names");    \
  const RCP<const Names>& inp_names =                                         \
    Teuchos::rcp(new const charon::Names(1,names->prefix(),                   \
                 names->discfields(),names->discsuffix()));                   \

///////////////////////////////////////////////////////////////////////////////
//
//  buildClosureModels()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
Teuchos::RCP<std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
charon::ClosureModelFactory<EvalT>::buildClosureModels(
  const std::string&                           modelId,
  const Teuchos::ParameterList&                models,
  const panzer::FieldLayoutLibrary&            fl,
  const Teuchos::RCP<panzer::IntegrationRule>& ir,
  const Teuchos::ParameterList&                defaults,
  const Teuchos::ParameterList&                userData,
  const Teuchos::RCP<panzer::GlobalData>&      globalData,
  PHX::FieldManager<panzer::Traits>&           fm) const
{
  CHARON_USINGS_ETC(DEFINE_NONE);
  using charon::Material_Properties;
  EvaluatorVector evaluators = rcp(new vector<RCP<Evaluator<Traits>>>);
  if (not models.isSublist(modelId))
  {
    stringstream msg;
    msg << "Failed to find requested model, \"" << modelId
        << "\" for equation set:\n" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(not models.isSublist(modelId),
      std::logic_error, msg.str());
  }
  const string& modelName(models.name());
  const ParameterList& myModels(models.sublist(modelId));
 

  // **************************************************************************
  // Build and Register Charon Parameters->Closure Models.
  // **************************************************************************
  if (modelName == m_model_name_prefix + "Closure Models")
  {
    // Mandatory to specify material name for Closure Models, except for global
    // sections.
    if (not myModels.isParameter("Material Name"))
    {
      if (m_throw_if_model_not_found)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(not myModels.isParameter("Material Name"),
          std::logic_error, "Material Name must be specified!");
      }
      else
        return evaluators;
    }

    // Get the Material Name and validate it.
    Material_Properties& matProperty = Material_Properties::getInstance();
    const string& matName(myModels.get<string>("Material Name"));
    matProperty.validateMaterialName(matName);
   
    // setup mole fraction materials
    setupMoleFraction(myModels);

    // Make sure charon::Names is available.
    if (not defaults.isType<RCP<const charon::Names>>("Names"))
    {
      if (m_throw_if_model_not_found)
      {
        stringstream msg;
        msg << "charon::ClosureModelFactory failed to build evaluator "
            << "\n in model \"" << modelId
            << "\".  No charon::Names object." << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION(false, std::logic_error, msg.str());
      }
      else
        return evaluators;
    }
    const RCP<const Names>& names = defaults.get<RCP<const Names>>("Names");
    const RCP<const Names>& inp_names = 
      Teuchos::rcp(new const charon::Names(1,names->prefix(),names->discfields(),names->discsuffix()));
    // Note for harmonic balance support: 
    // when referring to input deck naming, use inp_names
    // when referring to evaluated field names, use names; these fields will have the _TPi_ suffix required to perform HB

    // Obtain information on the input equation set
    const ParameterList& iesParams(defaults.sublist("Options"));
    string eqnSetType(defaults.get<string>("Type"));

    // ACH: TESTING
    string fd_suffix = m_fd_suffix;
    // NOTE: We don't do the following:
    // string fd_suffix = (eqnSetType == "Frequency Domain" ? iesParams.get<std::string>("Frequency Domain Suffix") : "");
    // because that string is reserved for suffixing the DOFs, etc. 
    // We handle the fd_suffix for the Closure Model through the Closure Model constructor, using a different charon::Names object.
    //std::cout << "ClosureModelFactory has fd_suffix value: " << fd_suffix << std::endl;
    string timeDomEqnSet = eqnSetType;
    eqnSetType = (eqnSetType == "Frequency Domain" ? iesParams.get<std::string>("Time Domain Equation Set") : timeDomEqnSet);

    // When BGN = On, Intrinsic Concentration MUST be either Old Slotboom, Persson,
    // Slotboom, or Harmon!
    string bgn("Off");
    if (iesParams.isParameter("Band Gap Narrowing"))
      bgn = iesParams.get<string>("Band Gap Narrowing");
    if (bgn == "On")
    {
      if (not myModels.isSublist(inp_names->field.intrin_conc))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Intrinsic Concentration MUST be specified when Band Gap "   \
          "Narrowing = On");
      const ParameterList& niParamList(
        myModels.sublist(inp_names->field.intrin_conc));
      if (niParamList.isType<double>("Value"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Intrinsic Concentration Value MUST be string when Band "    \
          "Gap Narrowing = On");
      const string& niModel(niParamList.get<string>("Value"));
      if ((niModel != "Old Slotboom"  )  and
          (niModel != "Persson"  )  and
          (niModel != "Slotboom")  and
          (niModel != "Harmon"  ))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Intrinsic Concentration Value MUST be Old Slotboom, Persson, "     \
          "Slotboom, or Harmon when Band Gap Narrowing = On");
    } // end if (bgn == "On")

    // Get the Incomplete Ionization information.
    string withAccIncompleteIoniz("Off"), withDonIncompleteIoniz("Off");
    bool withIonizAcc(false), withIonizDon(false);
    if (iesParams.isParameter("Acceptor Incomplete Ionization"))
      withAccIncompleteIoniz =
        iesParams.get<string>("Acceptor Incomplete Ionization");
    if (iesParams.isParameter("Donor Incomplete Ionization"))
      withDonIncompleteIoniz =
        iesParams.get<string>("Donor Incomplete Ionization");
    if ((withAccIncompleteIoniz == "On"                       )  and
        (not myModels.isSublist("Incomplete Ionized Acceptor")))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error! An Acceptor Incomplete Ionization model MUST be specified "   \
        "when Acceptor Incomplete Ionization = On\n");
    if ((withDonIncompleteIoniz == "On"                    )  and
        (not myModels.isSublist("Incomplete Ionized Donor")))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error! A Donor Incomplete Ionization model MUST be specified when "  \
        "Donor Incomplete Ionization = On\n");

    // Check the Incomplete Ionization Acceptor sublist.
    if ((withAccIncompleteIoniz == "On"                   )  and
        (myModels.isSublist("Incomplete Ionized Acceptor")))
    {
      const ParameterList& incmpl_param =
        myModels.sublist("Incomplete Ionized Acceptor");
      if (not incmpl_param.isType<string>("Value"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! "          \
          "Incomplete Ionization Acceptor MUST specify a string Value");
      if (incmpl_param.get<string>("Value") != "Model")
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Value MUST be Model");
      if (not incmpl_param.isSublist("Model"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor MUST define a Model");
      const ParameterList& incmpl_param_def = incmpl_param.sublist("Model");
      if (not incmpl_param_def.isType<double>("Critical Doping Value"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Model MUST specify a "       \
          "Critical Doping Value of double type");
      if (not incmpl_param_def.isType<double>("Degeneracy Factor"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Model MUST specify a "       \
          "Degeneracy Factor of double type");
      if ((not incmpl_param_def.isType<double>("Ionization Energy")  )  and
          (not incmpl_param_def.isType<string>("AccIncmplIoniz File")))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Model MUST specify an"       \
          "Ionization Energy of double type or a two-column table in a file");
      if (not incmpl_param_def.isType<string>("Approximation"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Model MUST specify an "      \
          "Approximation of string type");
      withIonizAcc = true;
    } // end of "Check the Incomplete Ionization Acceptor sublist".

    // Check the Incomplete Ionization Donor sublist.
    if ((withDonIncompleteIoniz == "On"                )  and
        (myModels.isSublist("Incomplete Ionized Donor")))
    {
      const ParameterList& incmpl_param =
        myModels.sublist("Incomplete Ionized Donor");
      if (not incmpl_param.isType<string>("Value"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor MUST specify a string Value");
      if (incmpl_param.get<string>("Value") != "Model")
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor Value MUST be Model");
      if (not incmpl_param.isSublist("Model"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor MUST define a Model");
      const ParameterList& incmpl_param_def = incmpl_param.sublist("Model");
      if (not incmpl_param_def.isType<double>("Critical Doping Value"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor Model MUST specify a "          \
          "Critical Doping Value of double type");
      if (not incmpl_param_def.isType<double>("Degeneracy Factor"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor Model MUST specify a "          \
          "Degeneracy Factor of double type");
      if ((not incmpl_param_def.isType<double>("Ionization Energy")  )  and
          (not incmpl_param_def.isType<string>("DonIncmplIoniz File")))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Donor Model MUST specify an "         \
          "Ionization Energy of double type or a two-column table in a file");
      if (not incmpl_param_def.isType<string>("Approximation"))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Error! Incomplete Ionization Acceptor Model MUST specify an "      \
          "Approximation of string type");
      withIonizDon = true;
    } // end of "Check the Incomplete Ionization Donor sublist".

    // determine the driving force for Avalanche, TrapSRH, and B2BT models
    string drForce("EffectiveField");
    if (iesParams.isParameter("Driving Force"))
      drForce = iesParams.get<string>("Driving Force");

    // driving force checking for Avalanche model
    if (iesParams.isParameter("Avalanche") and
      iesParams.get<string>("Avalanche") == "On") {
      // when the avalanche model is turned on, check if its driving force is the same
      // with the one in Physics section
      if (not myModels.isSublist(inp_names->field.avalanche_rate)) {
        // in this case the default ava model will the created with a
        // default driving force
        if (drForce != "EffectiveField")
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Error! Driving Force in Physics Block must coincide with Avalanche model  " \
            "Driving Force");
      } else {
        // retrieve avalanche model params
        const ParameterList& ava_param_def = myModels.sublist(inp_names->field.avalanche_rate);
        if (not ava_param_def.isParameter("Driving Force")) {
          // default driving force for Ava model
          if (drForce != "EffectiveField")
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "Error! Driving Force in Physics Block must coincide with Avalanche model  " \
               "Driving Force");
        } else {
          if ( not ava_param_def.isType<string>("Driving Force"))
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
               "Error! Avalanche model Driving Force must be a string");
          // check if ava model driving force coincide with the one in Physics block
          std::map<string, string> drForces = {
            {"GradQuasiFermi", "GradQuasiFermi"},
            {"EffectiveFieldParallelJ", "EffectiveField"},
            {"EffectiveFieldParallelJtot", "EffectiveField"},
            {"GradPotentialParallelJ", "GradPotential"},
            {"GradPotentialParallelJtot", "GradPotential"}
          };
          if (drForces[ava_param_def.get<string>("Driving Force")] != drForce)
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Error! Driving Force in Physics Block must coincide with Avalanche model  " \
            "Driving Force");
        }
      }
    }

    // determine the driving force for Band2Band Tunneling models
    if (iesParams.isParameter("Band2Band Tunneling") and
        iesParams.get<string>("Band2Band Tunneling") == "On") 
    {
      // when the Band2Band Tunneling model is turned on, check if its driving force is consistent
      // with the one in Physics section
      if (myModels.isSublist(names->field.bbt_rate)) 
      {
        // retrieve Band2Band Tunneling model params
        const ParameterList& bbt_param_def = myModels.sublist(names->field.bbt_rate);
        
        if (not bbt_param_def.isType<string>("Driving Force"))
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
            "Error! Band2Band Tunneling model Driving Force must be a string");

        // check if Band2Band Tunneling model driving force coincide with the one in Physics block
        std::map<string, string> drForces_bbt = {
          {"GradQuasiFermi", "GradQuasiFermi"},
          {"EffectiveField", "EffectiveField"},
          {"GradPotential", "GradPotential"},
          {"EffectiveFieldParallelJ", "EffectiveField"},
          {"EffectiveFieldParallelJtot", "EffectiveField"},
          {"GradPotentialParallelJ", "GradPotential"},
          {"GradPotentialParallelJtot", "GradPotential"}
        };

        if (drForces_bbt[bbt_param_def.get<string>("Driving Force")] != drForce)
           TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
           "Error! Driving Force in Physics Block must coincide with Band2Band Tunneling model Driving Force");
      }
    }

    // Check if Band Gap and Electron Affinity are given in the input xml.
    bool isBandGap(false), isTempDepBG(false), isAffinity(false);
    if (myModels.isSublist("Band Gap"))
    {
      isBandGap = true;
      const ParameterList& bgParamList(myModels.sublist("Band Gap"));
      if ((bgParamList.isType<string>("Value")          )  and
          (bgParamList.get<string>("Value") == "TempDep"))
        isTempDepBG = true;
    }
    if (myModels.isSublist("Electron Affinity"))
      isAffinity = true;

    // Obtain the values of bUseFD and FDFormula.
    string isFermiDirac("False");
    if (iesParams.isParameter("Fermi Dirac"))
      isFermiDirac = iesParams.get<string>("Fermi Dirac");
    bool bUseFD(false);
    if (isFermiDirac == "True")
      bUseFD = true;
    string FDFormula("Schroeder");
    if (iesParams.isParameter("FD Formula"))
      FDFormula = iesParams.get<string>("FD Formula");


    // Check if eqnSetType = Lattice and Heat Generation = Analytic.
    bool isHGAnalytic(false);
    if (eqnSetType == "Lattice")
    {
      string heatGenType = iesParams.get<string>("Heat Generation");
      if (heatGenType == "Analytic")
        isHGAnalytic = true;
    }

    // Check if eqnSetType is in eqSetList1 or eqSetList2.
    RCP<vector<string>> eqSetList1 = rcp(new vector<string>);
    RCP<vector<string>> eqSetList2 = rcp(new vector<string>);
    eqSetList1->push_back("Laplace"                  );
    eqSetList1->push_back("SGCVFEM Laplace"          );
    eqSetList1->push_back("NLPoisson"                );
    eqSetList1->push_back("SGCVFEM NLPoisson"        );
    eqSetList1->push_back("Drift Diffusion"          );
    eqSetList1->push_back("EFFPG Drift Diffusion"    );
    eqSetList1->push_back("SGCVFEM Drift Diffusion"  );
    eqSetList1->push_back("SGCharon1 Drift Diffusion");
    eqSetList1->push_back("Frequency Domain"         );
    //eqSetList2->push_back("Lattice"                  );
    eqSetList2->push_back("DDIon"                    );
    eqSetList2->push_back("DDLattice"                );
    eqSetList2->push_back("DDIonLattice"             );
    eqSetList2->push_back("EFFPG DDIonLattice"       );

    bool inEqSetList1(false);
    for (size_t eq(0); eq < eqSetList1->size(); ++eq)
    {
      if (eqnSetType == eqSetList1->at(eq))
        inEqSetList1 = true;
    }

    bool inEqSetList2(false);
    for (size_t eq(0); eq < eqSetList2->size(); ++eq)
    {
      if (eqnSetType == eqSetList2->at(eq))
        inEqSetList2 = true;
    }

    bool isCVFEM(false);
    Teuchos::ParameterList cvfem_data;
    cvfem_data.set<bool>("Is CVFEM",isCVFEM);
    if ((eqnSetType == "SGCVFEM Laplace"          )  or
        (eqnSetType == "SGCVFEM NLPoisson"        )  or
        (eqnSetType == "SGCVFEM Drift Diffusion"  )  or
        (eqnSetType == "SGCharon1 Drift Diffusion"))
    {
        isCVFEM = true;

        // define basis and IR for subcontrol vol centroid
        Teuchos::RCP<panzer::BasisIRLayout> basis =
          defaults.get<Teuchos::RCP<panzer::BasisIRLayout> >("Basis");
        Teuchos::RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();

        panzer::CellData celldata(basis->numCells(),cellTopoInfo->getCellTopology() );

        std::string cvfem_type = "volume";
        Teuchos::RCP<panzer::IntegrationRule> cvfem_vol_ir =
          Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));

        Teuchos::RCP<const panzer::PureBasis> Hgradbasis =
          Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),
                                             cellTopoInfo->getCellTopology()));
        Teuchos::RCP<panzer::BasisIRLayout> hgrad_vol_cvfem =
          Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_vol_ir));
	Teuchos::RCP<const panzer::PureBasis> Hcurlbasis =
	  Teuchos::rcp(new panzer::PureBasis("HCurl",1,basis->numCells(),
					     cellTopoInfo->getCellTopology()));
	Teuchos::RCP<panzer::BasisIRLayout> hcurl_vol_cvfem =
          Teuchos::rcp(new panzer::BasisIRLayout(Hcurlbasis,*cvfem_vol_ir));

        // Add to parameter list
        cvfem_data.set<Teuchos::RCP<panzer::BasisIRLayout>>("CVFEM Vol Basis", hgrad_vol_cvfem);
        cvfem_data.set<Teuchos::RCP<panzer::IntegrationRule>>("CVFEM Vol IR", cvfem_vol_ir);
	cvfem_data.set<Teuchos::RCP<panzer::BasisIRLayout>>("CVFEM Vol HCurlBasis", hcurl_vol_cvfem);
        cvfem_data.set<bool>("Is CVFEM",isCVFEM);
    }

    // Check if Solve Ion = True or False.
    bool bSolveIon(false);
    if (iesParams.isParameter("Solve Ion"))
    {
      string solveIon = iesParams.get<string>("Solve Ion");
      if (solveIon == "True")
        bSolveIon = true;
    }

    // Obtain Ion Charge when eqnSetType = "DDIon" or "DDIonLattice" or
    // "EFFPG DDIonLattice".
    int ionCharge(1);
    if (iesParams.isParameter("Ion Charge"))
      ionCharge = iesParams.get<int>("Ion Charge");

    //-------------------------------------------------------------------------
    // Set up the empirical model if it is called for.
    //-------------------------------------------------------------------------

    // Need to add a conditional to ensure this only happens once.               // JMG:  What needs to happen here?  12/16/2016.

    // Build the Empirical_Damage_Model evaluator when empirical model = On in
    // the Physics Blocks.
    ParameterList EmpiricalModelParameters;
    string empiricalDefectRecomb("Off");
    if (iesParams.isParameter("Empirical Defect"))
      empiricalDefectRecomb = iesParams.get<string>("Empirical Defect");
    if (empiricalDefectRecomb == "On")
      createRecombRateEmpiricalDefect(evaluators, defaults, *ir, models, userData);

    // Loop over the models for a given modelId (e.g., Silicon Parameter).
    for (ParameterList::ConstIterator modelIt = myModels.begin();
      modelIt != myModels.end(); ++modelIt)
    {
      bool found(false);
      // freqDom note: when comparing to strings from the input deck, use key_pure
      //               when naming a Field, we should use key everywhere internally
      const string key_pure(modelIt->first); 
      const string key(key_pure + m_fd_suffix);
      if ((key == names->key.material_name                  )  or
          (key == names->key.radiative_recombination        )  or
          (key == names->key.auger_recombination            )  or
          (key == names->key.trap_srh_recombination         )  or
	  (key == names->key.dynamic_traps_recombination    )  or
          (key == names->key.defect_cluster_recombination   )  or
          (key == names->key.empirical_defect_recombination )  or
          (key == names->key.particle_strike                )  or
          (key == names->key.incomplete_ionized_acceptor    )  or
          (key == names->key.incomplete_ionized_donor       )  or
	  (key == names->key.tid) ) 
        found = true;

      if (key != names->key.material_name)
      {
	//std::cout << "The frequency domain suffix is currently: '" << m_fd_suffix << "'" << std::endl;
	//std::cout << "The current key_pure value is: '" << key_pure 
        //          << "' and the modified current key value is: '" << key << "'." << std::endl; 
        //std::cout << "For reference, 'names->field.doping' has value: " << names->field.doping << std::endl;
        const Teuchos::ParameterEntry& entry(modelIt->second);
        const ParameterList& plist(Teuchos::getValue<ParameterList>(entry));

        if (plist.isType<double>("Value"))
        {
          // Constant Relative Permittivity is specified (unitless, no
          // scaling).
          if (key == names->field.rel_perm) {
            /*
            found = createConstant(evaluators, defaults, names->field.rel_perm,
              plist.get<double>("Value"));
	    */
	    if (!matProperty.hasMoleFracDependence(matName)) {
	      found = createConstant(evaluators, defaults, names->field.rel_perm,
				     plist.get<double>("Value"));
 	    } else {
	      // with mole dependence
	      found = createRelPermittivity(evaluators, defaults, matName, myModels);
	    }
	  }

          // Constant Band Gap is specified (in [eV], no scaling).
          if (key == names->field.band_gap)
          {
            found = createConstant(evaluators, defaults, names->field.band_gap,
              plist.get<double>("Value"));
        
            // set eff_band_gap to band_gap for insulators
            const string matType(matProperty.getMaterialType(matName));
            if (matType == "Insulator")
              found &= createConstant(evaluators, defaults, names->field.eff_band_gap,
                plist.get<double>("Value"));
          }

          // Constant Electron Affinity is specified (in [eV], no scaling).
          if (key == names->field.affinity)
          {
            found = createConstant(evaluators, defaults, names->field.affinity,
              plist.get<double>("Value"));

            // set eff_affinity to affinity for insulators
            const string matType(matProperty.getMaterialType(matName)); 
            if (matType == "Insulator")
              found &= createConstant(evaluators, defaults, names->field.eff_affinity,
                plist.get<double>("Value"));
          }

          // Constant Intrinsic Concentration is specified.
          if (key == names->field.intrin_conc)
            found = createIntrinsicConcentration(evaluators, defaults,
              plist.get<double>("Value"));

          // Constant Electron Mobility is specified (needs to be scaled, so
          // use a Default evaluator).
          if (key == names->field.elec_mobility)
	  {
            found = createMobility(evaluators, defaults, CARRIER_TYPE_ELECTRON,
              plist.get<double>("Value"));
	  }

          // Constant Hole Mobility is specified.
          if (key == names->field.hole_mobility)
            found = createMobility(evaluators, defaults, CARRIER_TYPE_HOLE,
              plist.get<double>("Value"));

          // Constant Ion Mobility is specified.
          if (key == names->field.ion_mobility)
            found = createMobility(evaluators, defaults, CARRIER_TYPE_ION,
              plist.get<double>("Value"), ionCharge);

          // Constant SRH Electron Lifetime.
          if (key == names->field.elec_lifetime)
            found = createSRHLifetime(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, plist.get<double>("Value"));

          // Constant SRH Hole Lifetime.
          if (key == names->field.hole_lifetime)
            found = createSRHLifetime(evaluators, defaults, CARRIER_TYPE_HOLE,
              plist.get<double>("Value"));

          // Set the property value in charon::Material_Properties.
          double propertyValue(plist.get<double>("Value"));
          matProperty.setPropertyValue(matName, key_pure, propertyValue);
          // freqDom note: we use the unsuffixed key_pure (instead of key) to grab the name from the PL
        } // end of if (plist.isType<double>("Value"))

        if (plist.isType<string>("Value"))
        {
          const string value(plist.get<string>("Value"));

          // Check to see it is user data.                                       // JMG:  What needs to happen here?  12/16/2016.
          if (value == "User Data")
          {
            //TODO
          }

          //*******************************************************************
          // If it's not user data, then try direct key/value pairs.
          //*******************************************************************

          // Gather ION_DENSITY from Exodus when bSolveIon == false.
          else if ((value == "Exodus File"        )  and
                   (key   == names->dof.iondensity)  and
                   (not bSolveIon                 ))
          {
            found = createGatherFields(evaluators, defaults, key, userData);
            createDOF(evaluators, defaults, key, ir);
          }


          // Save Initial GRAD_ELECTRIC_POTENTIAL to a Phalanx field
          else if ((key == names->dof.phi) and (value == "Initial PotentialGrad"))
            found = createInitialPotentialGrad(evaluators, defaults, userData, cvfem_data);  

          //*******************************************************************
          // Instantiate doping related evaluators
          //*******************************************************************

          // Build the Doping_StepJunction evaluator.
          else if ((key == names->field.doping) and (value == "Step Junction"))
            found = createDopingStepJunction(evaluators, defaults, plist,
              withIonizAcc, withIonizDon, myModels);

          // Build the Doping_Function evaluator.
          //else if ((key == "Doping") and (value == "Function"))  // ACH: Why is this not names->field.doping?
          else if ((key_pure == "Doping") and (value == "Function")) // ACH: this is equivalent to doing 'key == "Doping" ' with un-fd-suffixed key
	    {
            found = createDopingFunction(evaluators, defaults, withIonizAcc,
                                         withIonizDon, globalData, userData, myModels);
	    }

          //*******************************************************************
          // Instantiate Permittivity related evaluators
          //*******************************************************************
          else if ((key == names->field.rel_perm) and (value == "Nitride"))
            found = createPermittivityNitride(evaluators, defaults, matName, myModels);
          //*******************************************************************
          // Instantiate band structure related evaluators
          //*******************************************************************

          // Build the BandGap_Nitride evaluator.
          else if ((key == names->field.band_gap) and (value == "Nitride"))
            found = createBandGapNitride(evaluators, defaults, matName, myModels);

          // Build the BandGap_TempDep evaluator.
          else if ((key == names->field.band_gap) and (value == "TempDep"))
            found = createBandGapTempDep(evaluators, defaults, matName,
              isAffinity, myModels);

          // Build the IntrinsicConc_OldSlotboom evaluator.
          else if ((key == names->field.intrin_conc) and (value == "Old Slotboom"))
            found = createIntrinsicConcOldSlotboom(evaluators, defaults, matName,
              bgn, myModels);

          // Build the IntrinsicConc_Persson evaluator.
          else if ((key == names->field.intrin_conc) and (value == "Persson"))
            found = createIntrinsicConcPersson(evaluators, defaults, matName,
              bgn, myModels);

          // Build the IntrinsicConc_Slotboom evaluator.
          else if ((key == names->field.intrin_conc) and (value == "Slotboom"))
            found = createIntrinsicConcSlotboom(evaluators, defaults, matName,
              bgn, myModels);

          // Build the IntrinsicConc_Harmon evaluator.
          else if ((key == names->field.intrin_conc) and (value == "Harmon"))
            found = createIntrinsicConcHarmon(evaluators, defaults, bgn,
              myModels);

          // Build the EffectiveDOS_Simple evaluator.
          else if ((key == "Effective DOS") and (value == "Simple")) 
            found = createEffectiveDOSSimple(evaluators, defaults, matName,
              myModels);

          // Build the Effective_DOS_Nitride evaluator.
          else if ((key == "Effective DOS") and (value == "Nitride"))
            found = createEffectiveDOSNitride(evaluators, defaults, matName,
              myModels);

          //*******************************************************************
          // Instantiate SRH lifetime related evaluators.
          //*******************************************************************

          // Build the SRHLifetime_Function evaluator for Electron.
          else if ((key   == names->field.elec_lifetime)  and
                   (value == "Function"                ))
            found = createSRHLifetimeFunction(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the SRHLifetime_Function evaluator for Hole.
          else if ((key   == names->field.hole_lifetime)  and
                   (value == "Function"                ))
            found = createSRHLifetimeFunction(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          //*******************************************************************
          // Instantiate low-field mobility related evaluators.
          //*******************************************************************

          // Build the Mobility_Analytic evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "Analytic"                ))
            found = createMobilityAnalytic(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Mobility_Analytic evaluator for Hole.
          else if ((key   == names->field.hole_mobility)  and
                   (value == "Analytic"                ))
            found = createMobilityAnalytic(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          // Build the Mobility_Arora evaluator for Electron.
          else if ((key == names->field.elec_mobility) and (value == "Arora"))
            found = createMobilityArora(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels, false);

          // Build the Mobility_Arora evaluator for Hole.
          else if ((key == names->field.hole_mobility) and (value == "Arora"))
            found = createMobilityArora(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels, false);

          // Build the Mobility_Masetti evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "Masetti"                 ))
            found = createMobilityMasetti(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Mobility_Masetti evaluator for Hole.
          else if ((key   == names->field.hole_mobility)  and
                   (value == "Masetti"                 ))
            found = createMobilityMasetti(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          // Build the Mobility_UniBo evaluator for Electron.
          else if ((key == names->field.elec_mobility) and (value == "UniBo"))
            found = createMobilityUniBo(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Mobility_UniBo evaluator for Hole.
          else if ((key == names->field.hole_mobility) and (value == "UniBo"))
            found = createMobilityUniBo(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          // Build the Mobility_DopantTempDep evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "DopantTempDep"           ))
            found = createMobilityDopantTempDep(evaluators, defaults,
              bSolveIon, myModels);

          // Build the Mobility_MOSFET evaluator for Electron.
          else if ((key == names->field.elec_mobility) and (value == "MOSFET"))
	    found = createMobilityMOSFET(evaluators, defaults,
	      CARRIER_TYPE_ELECTRON, matName, myModels);


          // Build the Mobility_MOSFET evaluator for Hole.
          else if ((key == names->field.hole_mobility) and (value == "MOSFET"))
            found = createMobilityMOSFET(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

           //*******************************************************************
          // Instantiate high-field mobility related evaluators.
          //*******************************************************************

          // Build the Mobility_PhilipsThomas evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "Philips-Thomas"          ))
            found = createMobilityPhilipsThomas(evaluators, defaults,
	      CARRIER_TYPE_ELECTRON, matName, myModels,false);

          // Build the Mobility_PhilipsThomas evaluator for Hole.
          else if ((key   == names->field.hole_mobility)  and
                   (value == "Philips-Thomas"          ))
            found = createMobilityPhilipsThomas(evaluators, defaults,
	      CARRIER_TYPE_HOLE, matName, myModels,false);

          // Build the Mobility_Lucent evaluator for Electron.
          else if ((key == names->field.elec_mobility) and (value == "Lucent"))
            found = createMobilityLucent(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Mobility_Lucent evaluator for Hole.
          else if ((key == names->field.hole_mobility) and (value == "Lucent"))
            found = createMobilityLucent(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          // Build the Mobility_GaAs evaluator for Electron.
          else if ((key == names->field.elec_mobility) and (value == "GaAs"))
            found = createMobilityGaAs(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Mobility_GaAs evaluator for Hole.
          else if ((key == names->field.hole_mobility) and (value == "GaAs"))
            found = createMobilityGaAs(evaluators, defaults, CARRIER_TYPE_HOLE,
              matName, myModels);

          // Build the Albrecht mobility evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "Albrecht"                ))
            found = createMobilityAlbrecht(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Albrecht mobility evaluator for Hole.
          else if ((key   == names->field.hole_mobility)  and
                   (value == "Albrecht"                ))
            found = createMobilityAlbrecht(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          // Build the Farahmand mobility evaluator for Electron.
          else if ((key   == names->field.elec_mobility)  and
                   (value == "Farahmand"               ))
            found = createMobilityFarahmand(evaluators, defaults,
              CARRIER_TYPE_ELECTRON, matName, myModels);

          // Build the Farahmand mobility evaluator for Hole.
          else if ((key   == names->field.hole_mobility)  and
                   (value == "Farahmand"               ))
            found = createMobilityFarahmand(evaluators, defaults,
              CARRIER_TYPE_HOLE, matName, myModels);

          //*******************************************************************
          // Instantiate Avalanche generation related evaluators.
          //*******************************************************************

          // Build the Avalanche_vanOverstraeten evaluator.
          else if ((key   == names->field.avalanche_rate)  and
                   (value == "vanOverstraeten"          ))
            found = createAvalancheVanOverstraeten(evaluators, defaults,
              matName, cvfem_data, myModels);

          // Build the Avalanche_Okuto evaluator.
          else if ((key == names->field.avalanche_rate) and (value == "Okuto"))
            found = createAvalancheOkuto(evaluators, defaults, matName,
              myModels, cvfem_data);

          // Build the Avalanche_Lackner evaluator.
          else if ((key   == names->field.avalanche_rate)  and
                   (value == "Lackner"                  ))
            found = createAvalancheLackner(evaluators, defaults, matName,
              myModels, cvfem_data);

          // Build the Avalanche_UniBo evaluator.
          else if ((key == names->field.avalanche_rate) and (value == "UniBo"))
            found = createAvalancheUniBo(evaluators, defaults, matName,
              myModels, cvfem_data);

          // Build the Avalanche_UniBoNew evaluator.
          else if ((key   == names->field.avalanche_rate)  and
                   (value == "UniBoNew"               ))
            found = createAvalancheUniBoNew(evaluators, defaults, matName,
              myModels, cvfem_data);

          // Build the Avalanche_Selberherr evaluator.
          else if ((key   == names->field.avalanche_rate)  and
                   (value == "Selberherr"                   ))
            found = createAvalancheSelberherr(evaluators, defaults, matName,
              eqnSetType, myModels, cvfem_data);

	  // Build the Avalanche_CrowellSze evaluator.
	  else if ((key   == names->field.avalanche_rate)  and
                   (value == "Crowell-Sze"              ))
            found = createAvalancheCrowellSze(evaluators, defaults, matName,
              eqnSetType, myModels, cvfem_data);

          //*******************************************************************
          // Instantiate Band2Band Tunneling generation related evaluators.
          //*******************************************************************

           // Build the Band2Band_Tunneling_Local evaluator.
           else if ((key == names->field.bbt_rate) and
                    (value == "Local"))
             found = createBand2BandTunnelingLocal(evaluators, defaults, matName,
               eqnSetType, myModels, cvfem_data);

          //*******************************************************************
          // Instantiate thermal related evaluators.
          //*******************************************************************

          // Build the HeatCapacity_TempDep evaluator.
          else if ((key == names->field.heat_cap) and (value == "TempDep"))
            found = createHeatCapacityTempDep(evaluators, defaults, matName,
              myModels);

          // Build the HeatCapacity_PowerLawTempDep evaluator.
          else if ((key == names->field.heat_cap) and (value == "PowerLawTempDep"))
            found = createHeatCapacityPowerLawTempDep(evaluators, defaults, matName,
              myModels);

          // Build the ThermalConduct_TempDep evaluator.
          else if ((key == names->field.kappa) and (value == "TempDep"))
            found = createThermalConductTempDep(evaluators, defaults, matName,
              myModels);

          // Build the ThermalConduct_PowerLawTempDep evaluator.
          else if ((key == names->field.kappa) and (value == "PowerLawTempDep"))
            found = createThermalConductPowerLawTempDep(evaluators, defaults, matName,
              myModels);

          // Build the ThermalConduct_LinearTempDep evaluator.
          else if ((key == names->field.kappa) and (value == "LinearTempDep"))
            found = createThermalConductLinearTempDep(evaluators, defaults,
              myModels);

          // Build the ThermalConduct_LinearIonDep evaluator.
          else if ((key == names->field.kappa) and (value == "LinearIonDep"))
            found = createThermalConductLinearIonDep(evaluators, defaults,
              myModels);

          // Build the Analytic_HeatGeneration evaluator.
          else if ((key   == names->field.heat_gen)  and
                   (value == "Analytic"           )  and
                   (isHGAnalytic                  ))
            found = createAnalyticHeatGeneration(evaluators, defaults,
              myModels);

          // Build the SoretCoeff_TempDep evaluator.
          else if ((key   == names->field.ion_soret_coeff)  and
                   (value == "TempDep"                   ))
            found = createSoretCoeffTempDep(evaluators, defaults, matName,
              myModels);

          // Build the ThermodiffCoeff_Custom evaluator.
          else if ((key   == names->field.ion_thermodiff_coeff)  and
                   (value == "Custom"                         ))
            found = createThermodiffCoeffCustom(evaluators, defaults,
              myModels);

          //*******************************************************************
          // Instantiate the ion mobility evaluator (also compute ion velocity)
          // and the ion diff. coeff. evaluator.
          //*******************************************************************

          // Build the Mobility_RigidPointIon evaluator for Ion.
          else if ((key   == names->field.ion_mobility)  and
                   (value == "RigidPointIon"          ))
            found = createMobilityRigidPointIon(evaluators, defaults, matName,
              ionCharge, myModels);

          // Build the DiffCoeff_IonDep evaluator for Ion.
          else if ((key   == names->field.ion_diff_coeff)  and
                   (value == "IonDep"                   ))
            found = createDiffCoeffIonDep(evaluators, defaults, myModels);

          //*******************************************************************
          // Instantiate the optical generation evaluator.
          //*******************************************************************

          // Build the OptGen_Function evaluator.
          else if ((key == names->field.opt_gen) and (value == "Function"))
            found = createOptGenFunction(evaluators, defaults, myModels);

          //*******************************************************************
          // Instantiate the mole fraction evaluator.
          //*******************************************************************

          // Build the MoleFraction_Function evaluator.
          else if ((key == names->field.mole_frac) and (value == "Function"))
            found = createMoleFractionFunction(evaluators, defaults, matName, myModels);

          //*******************************************************************
          // Instantiate the bulk fixed charge evaluator.
          //*******************************************************************

          // Build the BulkFixCharge_Function evaluator.
          else if ((key == names->field.fixed_charge) and (value == "Function"))
            found = createBulkFixChargeFunction(evaluators, defaults, myModels, globalData, cvfem_data);

          //*******************************************************************
          // End of "If it's not user data, then try direct key/value pairs".
          //*******************************************************************

          // Make sure a model was found.
          if (not found)
          {
            stringstream msg;
            msg << "ClosureModelFactory doesn't know how to build a model "
                << "with:"                  << std::endl
                << "  model = " << modelId  << std::endl
                << "  key   = " << key      << std::endl
                << "  value = " << value    << std::endl;
            TEUCHOS_TEST_FOR_EXCEPTION(not found, std::logic_error, msg.str());
          }
        } // end if (plist.isType<string>("Value"))
      } // end if (key != "Material Name")

      // Throw an exception if we weren't able to find the model.
      if ((not found) and (m_throw_if_model_not_found))
      {
        stringstream msg;
        msg << "ClosureModelFactory failed to build evaluator for key \""
            << key << "\"" << std::endl << " in model \"" << modelId
            << "\".  Please correct the type or add support to the factory."
            << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION(not found, std::logic_error, msg.str());
      }
    } // end loop over the models for a given modelId (e.g., Silicon Parameter)

    //*************************************************************************
    // Register fields using values from charon::Material_Properties,
    // when they are not specified in input .xml file.
    //*************************************************************************

    // Register "Relative Permittivity" when not specified in the input .xml
    // file.
    //if (not myModels.isSublist(names->field.rel_perm))
    if (not myModels.isSublist("Relative Permittivity"))
    {
      const string matType(matProperty.getMaterialType(matName));
      if ((matType == "Semiconductor") or (matType == "Insulator"))
        createConstant(evaluators, defaults, names->field.rel_perm,
		       matProperty.getPropertyValue(matName, "Relative Permittivity"));
      //createConstant(evaluators, defaults, "Relative Permittivity"+m_fd_suffix,
      //	       //matProperty.getPropertyValue(matName, names->field.rel_perm));
      //	       matProperty.getPropertyValue(matName, "Relative Permittivity"));
    }

    //*************************************************************************
    // Instantiate default band structure related evaluators.
    //*************************************************************************

    // Register "Band Gap" when not specified in the input .xml.
    if (not myModels.isSublist("Band Gap"))
    {
      const string matType(matProperty.getMaterialType(matName));
      if ((matType == "Semiconductor") or (matType == "Insulator"))
      {
        createConstant(evaluators, defaults, names->field.band_gap,
          matProperty.getPropertyValue(matName, "Band Gap"));
        createConstant(evaluators, defaults, names->field.eff_band_gap,
          matProperty.getPropertyValue(matName, "Band Gap"));
      }
    }

    // When neither Electron Affinity nor TempDep Band Gap is given, use the
    // value from charon::Material_Properties.
    if ((not isAffinity) and (not isTempDepBG))
    {
      const string matType(matProperty.getMaterialType(matName));
      if ((matType == "Semiconductor") or (matType == "Insulator"))
      {
        createConstant(evaluators, defaults, names->field.affinity,
          matProperty.getPropertyValue(matName, "Electron Affinity"));
        createConstant(evaluators, defaults, names->field.eff_affinity,
          matProperty.getPropertyValue(matName, "Electron Affinity"));
      }
    }

    // When neither Band Gap nor Intrinsic Concentration is specified, use the
    // ni value from charon::Material_Properties.
    //if ((not myModels.isSublist(names->field.intrin_conc))  and
    if ((not myModels.isSublist("Intrinsic Concentration"))  and
        (not isBandGap                                   ))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createIntrinsicConcentration(evaluators, defaults,
				     //matProperty.getPropertyValue(matName, names->field.intrin_conc));
          matProperty.getPropertyValue(matName, "Intrinsic Concentration"));
    }

    // When Band Gap is specified but Intrinsic Concentration is NOT, calculate
    // ni = sqrt(Nc*Nv)*exp(-Eg/2kbT).
    if ((not myModels.isSublist(names->field.intrin_conc)) and (isBandGap))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createIntrinsicConcOldSlotboom(evaluators, defaults, matName, bgn);
    }

    // When Effective DOS is not specified, build the EffectiveDOS_Simple
    // evaluator using parameters from Material_Properties.
    if (not myModels.isSublist("Effective DOS"))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createEffectiveDOSSimple(evaluators, defaults, matName);
    }

    //*************************************************************************
    // Instantiate default mobility related evaluators.
    //*************************************************************************

    // Build the Mobility_Default evaluator for Electron.
    //if (not myModels.isSublist(names->field.elec_mobility))
    if (not myModels.isSublist("Electron Mobility"))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createMobility(evaluators, defaults, CARRIER_TYPE_ELECTRON,
	  //matProperty.getPropertyValue(matName, names->field.elec_mobility));
          matProperty.getPropertyValue(matName, "Electron Mobility"));
    }

    // Build the Mobility_Default evaluator for Hole.
    //if (not myModels.isSublist(names->field.hole_mobility))
    if (not myModels.isSublist("Hole Mobility"))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createMobility(evaluators, defaults, CARRIER_TYPE_HOLE,
	  //matProperty.getPropertyValue(matName, names->field.hole_mobility));
	  matProperty.getPropertyValue(matName,"Hole Mobility"));
    }

    // Build the DiffCoeff_Default evaluator for Electron.
    if (not myModels.isSublist(names->field.elec_diff_coeff))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createDiffCoeffDefault(evaluators, defaults, CARRIER_TYPE_ELECTRON,
          bUseFD, FDFormula);
    }

    // Build the DiffCoeff_Default evaluator for Hole.
    if (not myModels.isSublist(names->field.hole_diff_coeff))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createDiffCoeffDefault(evaluators, defaults, CARRIER_TYPE_HOLE, bUseFD,
          FDFormula);
    }

    // Build the Mobility_Default evaluator for Ion when bSolveIon = true.
    if ((not myModels.isSublist(names->field.ion_mobility)) and (bSolveIon))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createMobility(evaluators, defaults, CARRIER_TYPE_ION,
          matProperty.getPropertyValue(matName, names->field.ion_mobility),
          ionCharge);
    }

    // Build the DiffCoeff_Default evaluator for Ion when bSolveIon = true.
    if ((not myModels.isSublist(names->field.ion_diff_coeff)) and (bSolveIon))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
      {
        bool useFermiDirac(false);
        createDiffCoeffDefault(evaluators, defaults, CARRIER_TYPE_ION,
          useFermiDirac, FDFormula);
      }
    }

    //*************************************************************************
    // Instantiate default SRH lifetime related evaluators.
    //*************************************************************************

    string srhRecomb("Off");
    if (iesParams.isParameter("SRH"))
      srhRecomb = iesParams.get<string>("SRH");

    // Build the SRHLifetime_Constant evaluator for Electron.
    //if ((not myModels.isSublist(names->field.elec_lifetime))  and
    if ((not myModels.isSublist("Electron Lifetime"))  and
        (srhRecomb == "On"                                 ))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createSRHLifetimeConstant(evaluators, defaults, CARRIER_TYPE_ELECTRON,
	  //matProperty.getPropertyValue(matName, names->field.elec_lifetime));
          matProperty.getPropertyValue(matName, "Electron Lifetime"));
    }

    // Build the SRHLifetime_Constant evaluator for Hole.
    //if ((not myModels.isSublist(names->field.hole_lifetime))  and
    if ((not myModels.isSublist("Hole Lifetime"))  and
        (srhRecomb == "On"                                 ))
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Semiconductor")
        createSRHLifetimeConstant(evaluators, defaults, CARRIER_TYPE_HOLE,
	  //matProperty.getPropertyValue(matName, names->field.hole_lifetime));
          matProperty.getPropertyValue(matName, "Hole Lifetime"));
    }

    //*************************************************************************
    // Instantiate recombination and generation related evaluators.
    //*************************************************************************

    // Build the RecombRate_SRH evaluator when SRH = On in the Physics Blocks.
    if (srhRecomb == "On")
      createRecombRateSRH(evaluators, defaults, bUseFD);

    // Build the RecombRate_TrapSRH evaluator when Trap SRH = On in the Physics Blocks.
    string trapSrhRecomb("Off");
    if (iesParams.isParameter("Trap SRH"))
      trapSrhRecomb = iesParams.get<string>("Trap SRH");
    if (trapSrhRecomb == "On")
      createRecombRateTrapSRH(evaluators, defaults, matName, myModels, eqnSetType, drForce, cvfem_data);

    string dynamicTrapsRecomb("Off");
    if (iesParams.isParameter("Dynamic Traps"))
      dynamicTrapsRecomb = iesParams.get<string>("Dynamic Traps");
    if (dynamicTrapsRecomb == "On")
      createRecombRateDynamicTraps(evaluators, defaults, matName, myModels, eqnSetType, drForce, cvfem_data);

    // Build the RecombRate_Defect_Cluster evaluator when Defect Cluster = On
    // in the Physics Blocks.
    string defectClusterRecomb("Off");
    if (iesParams.isParameter("Defect Cluster"))
      defectClusterRecomb = iesParams.get<string>("Defect Cluster");
    if (defectClusterRecomb == "On")
      createRecombRateDefectCluster(evaluators, defaults, myModels);

    // Build the ionization_particleStrike evaluator when Particle Strike = On
    // in the Physics Blocks.
    string ionizationParticleStrike("Off");
    if (iesParams.isParameter("Particle Strike"))
      ionizationParticleStrike = iesParams.get<string>("Particle Strike");
    if (ionizationParticleStrike == "On")
      createIonizationParticleStrike(evaluators, defaults, myModels);

    string TID("Off");
    if (iesParams.isParameter("TID"))
      TID = iesParams.get<string>("TID");
    if (TID == "On") 
    {
      const string matType(matProperty.getMaterialType(matName));
      if (matType == "Insulator") {
	createTID(evaluators, defaults, matName, myModels, globalData, userData, cvfem_data);
      }
    }

    // Build the RecombRate_Radiative evaluator when Radiative = On.
    string radRecomb("Off");
    if (iesParams.isParameter("Radiative"))
      radRecomb = iesParams.get<string>("Radiative");
    if (radRecomb == "On")
      createRecombRateRadiative(evaluators, defaults, matName, bUseFD,
        myModels);

    // Build the RecombRate_Auger evaluator when Auger = On.
    string augerRecomb("Off");
    if (iesParams.isParameter("Auger"))
      augerRecomb = iesParams.get<string>("Auger");
    if (augerRecomb == "On")
      createRecombRateAuger(evaluators, defaults, matName, bUseFD, myModels);

    // When Avalanche Generation is NOT specified, use the vanOverstraeten
    // model as default.
    string avaGen("Off");
    if (iesParams.isParameter("Avalanche"))
      avaGen = iesParams.get<string>("Avalanche");
    if ((avaGen == "On"                                     )  and
        (not myModels.isSublist(names->field.avalanche_rate)))
      createAvalancheVanOverstraeten(evaluators, defaults, matName, cvfem_data);

    // When Band2Band Tunneling is NOT specified, throw  error message and stop simulation
    string bbtGen("Off");
    if (iesParams.isParameter("Band2Band Tunneling"))
      bbtGen = iesParams.get<string>("Band2Band Tunneling");
    if ((bbtGen == "On") and (not myModels.isSublist(names->field.bbt_rate)))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        std::endl << "Error: Band2Band Tunneling model is not found!");

    // Is Optical Generation on?
    string optGen("Off");
    if (iesParams.isParameter("Optical Generation"))
      optGen = iesParams.get<string>("Optical Generation");
    if ((optGen == "On") and (not myModels.isSublist(names->field.opt_gen)))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        std::endl << "Error: Optical Generation model is not found!");

    // Build the RecombRate_Total evaluator to sum up all recomb. rates.
    const string matType(matProperty.getMaterialType(matName));
    if (matType == "Semiconductor")
      createRecombRateTotal(evaluators, defaults, srhRecomb, trapSrhRecomb,
        defectClusterRecomb, empiricalDefectRecomb,
        ionizationParticleStrike, radRecomb, augerRecomb, optGen,
        avaGen, bbtGen, eqnSetType, cvfem_data);

    //*************************************************************************
    // Build the Intrinsic_FermiEnergy evaluator to make Ei available.
    //*************************************************************************
    if (matType == "Semiconductor")
      createIntrinsicFermiEnergy(evaluators, defaults);

    //*************************************************************************
    // Build the CondVale_Band evaluator to make Ec and Ev available.
    //*************************************************************************
    if ((matType == "Semiconductor") or (matType == "Insulator"))
      createCondValeBand(evaluators, defaults);

    //*************************************************************************
    // Build the Degeneracy_Factor evaluator to make the field available.
    //*************************************************************************
    if (matProperty.getMaterialType(matName) == "Semiconductor")
      createDegeneracyFactor(evaluators, defaults, bUseFD, FDFormula);

    //*************************************************************************
    // Build the Space_Charge evaluator to make the field available.
    //*************************************************************************
    if (matProperty.getMaterialType(matName) == "Semiconductor")
      createSpaceCharge(evaluators, defaults);

    //*************************************************************************
    // Instantiate material-independent evaluators.
    //*************************************************************************

    // Register "Lattice Temperature" when inEqSetList1 = true or eqnSetType =
    // "DDIon".
    if ((inEqSetList1) or (eqnSetType == "DDIon"))
      createLatticeTempConstant(evaluators, defaults, models);

    // Otherwise, throw out an error if Lattice Temperature is specified.
    else
    {
      const string key(names->field.latt_temp);
      if (models.isParameter(key))
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          "Lattice Temperature should not be given as it is position "
          << "dependent for the eqnset of " << eqnSetType << "!" << std::endl);
    }

    //*************************************************************************
    // Obtain the Reference Energy needed for heterogeneous devices, applicable
    // to drift-diffusion-only formulations.
    //*************************************************************************
    string refMaterial("");
    if (inEqSetList1) refMaterial = "Silicon";
    if (models.isParameter("Reference Material"))
      refMaterial = models.get<string>("Reference Material");
    bool foundRefMat(false);
    
    // Loop over modelIds to compute Reference Energy for all the modelIds.
    for (ParameterList::ConstIterator modelIt = models.begin();
      modelIt != models.end(); ++modelIt)
    {
      const string key(modelIt->first);

      // If the modelId specifies a sublist...
      if (models.isSublist(key))
      {
        const ParameterList& submodels = models.sublist(key);
        
        if ((submodels.isParameter("Material Name")               )  and
            (submodels.get<string>("Material Name") == refMaterial))
        {
          foundRefMat = true;
          createReferenceEnergy(evaluators, defaults, refMaterial, submodels);
        }

        else if ((submodels.isParameter("Material Name")         )  and
            (submodels.get<string>("Material Name") != refMaterial))
        {
          /*
          if (eqnSetType.find("Laplace") != string::npos) 
          {
            foundRefMat = true;
            createReferenceEnergy(evaluators, defaults, refMaterial, submodels);
          } 
          */
	  if (eqnSetType.find("Laplace") != string::npos) 
          {
            foundRefMat = true;
            createReferenceEnergy(evaluators, defaults, refMaterial, submodels);
          } else {
            if (matProperty.hasMoleFracDependence(matName)) {
	      matProperty.validateMaterialName(refMaterial);
	      foundRefMat = true;
	      createReferenceEnergy(evaluators, defaults, refMaterial, submodels);
	    }
	  }
        }

        // Parse the global MMS parameters (applied to all materials).  These
        // quantities are only relevant for scalar evaluations; no derivatives
        // are required.
        else if ((boost::iequals(key, "Global MMS Parameters"))  and
                 (typeid(EvalT) == typeid(Traits::Residual)   ))
        {

          for (ParameterList::ConstIterator smodIt=submodels.begin();
            smodIt != submodels.end(); ++smodIt)
          {
            ParameterList in =
              Teuchos::getValue<ParameterList>(smodIt->second);
            const string value(in.get<string>("Value"));

            if (boost::iequals(smodIt->first, "Global Statistics"))
              createGlobalStatistics(evaluators, in, value, ir, globalData,
                userData, fm);
            else if (boost::iequals(smodIt->first, "Analytic Solution"))
              createMMSAnalyticSolution(evaluators, defaults, value, ir, fl,
                modelId);
            else if (boost::iequals(smodIt->first, "Analytic Comparison"))
              createAnalyticComparison(evaluators, defaults, value, in, fm);
            else if (boost::iequals(smodIt->first,
              "Analytic Comparison: L2 Error"))
              createAnalyticComparisonL2Error(evaluators, defaults, value, in,
                ir, userData, fm);
            else if (boost::iequals(smodIt->first,
              "Analytic Comparison: Relative Error"))
              createAnalyticComparisonRelError(evaluators, defaults, value, in,
                fm);
            else
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                "Unknown Global MMS Parameter (List) \"" << smodIt->first <<
                "\" specified in Closure Models section of the input file");
          }
        } // end of else if for Global MMS Parameters section

        
        // Parse the global MMS parameters (applied to all materials).  These
        // quantities are only relevant for scalar evaluations; no derivatives
        // are required.
        else if ((boost::iequals(key, "Global Norm"))  and
                 (typeid(EvalT) == typeid(Traits::Residual)   ))
        {
          for (ParameterList::ConstIterator smodIt=submodels.begin();
            smodIt != submodels.end(); ++smodIt)
          {
            ParameterList in =
              Teuchos::getValue<ParameterList>(smodIt->second);
            const string value(in.get<string>("Value"));
            if (boost::iequals(smodIt->first, "Norm Calculation: L2"))
              createNormCalculationL2(evaluators, defaults, value, in,
                ir, userData, globalData, fm);
            else if (boost::iequals(smodIt->first, "Norm Calculation: H1"))
              createNormCalculationH1(evaluators, defaults, value, in,
                ir, userData, globalData, fm);
            else if (boost::iequals(smodIt->first, "Manufactured Solution"))
              createManufacturedSolution(evaluators, defaults, value, in,
                ir, userData, globalData, fm);
            else if (boost::iequals(smodIt->first, "Norm Calculation: L2 Error"))
              createNormCalculationL2Error(evaluators, defaults, value, in,
                ir, userData, globalData, fm);
            else if (boost::iequals(smodIt->first, "Norm Calculation: H1 Error"))
              createNormCalculationH1Error(evaluators, defaults, value, in,
                ir, userData, globalData, fm);
            else
              TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                "Unknown Norm Parameter (List) \"" << smodIt->first <<
                "\" specified in Closure Models section of the input file");
          }
        } // end of Global Norms section
      }  // end of if (models.isSublist(key))
    }  // end of loop over modeIds

    // Check to make sure the Reference Material was found.
    if ((not foundRefMat) and (inEqSetList1))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error: Reference Material is NOT found!");

    // Check to make sure the Reference Material was NOT specified when inEqSetList2 = true
    //if ((foundRefMat) and (inEqSetList2))
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    //    "Error: Reference Material should NOT be specified for " << eqnSetType << " !");

    // Check to make sure the Reference Material was found.
    if ((not foundRefMat) and (inEqSetList2))
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error: Reference Material is NOT found!");

    //*************************************************************************
    // Instantiate the default heat capacity and thermal conductivity
    // evaluators.
    //*************************************************************************
    string str("Lattice");
    size_t foundtemp(eqnSetType.find(str));

    // Build the HeatCapacity_TempDep evaluator when not given in the input
    // .xml file.
    if ((not myModels.isSublist(names->field.heat_cap))  and
        (foundtemp != string::npos                    ))
      createHeatCapacityTempDep(evaluators, defaults, matName);

    // Build the ThermalConduct_TempDep evaluator when not given in the input
    // .xml file.
    if ((not myModels.isSublist(names->field.kappa))  and
        (foundtemp != string::npos                 ))
      createThermalConductTempDep(evaluators, defaults, matName);

    // Build the SoretCoeff_TempDep evaluator when not given in the input .xml
    // file, and bSolveIon = true and contains Lattice.
    if ((not myModels.isSublist(names->field.ion_soret_coeff))  and
        (bSolveIon                                           )  and
        (foundtemp != string::npos                           ))
      createSoretCoeffTempDep(evaluators, defaults, matName);

    // Build the ThermodiffCoeff_Default evaluator when not given in the input
    // .xml file, and bSolveIon = true and contains Lattice.
    if ((not myModels.isSublist(names->field.ion_thermodiff_coeff))  and
        (bSolveIon                                                )  and
        (foundtemp != string::npos                                ))
      createThermodiffCoeffDefault(evaluators, defaults);

    //*************************************************************************
    // Build the FEM_GradNegPotential evaluator to make the field available.
    //*************************************************************************
    createFEMGradNegPotential(evaluators, defaults);
  }  // end of if (modelName == "Charon Parameter->Closure Models")

  // **************************************************************************
  // Build and Register Charon Parameters->Initial Conditions.
  // **************************************************************************
  else if (modelName == m_model_name_prefix + "Initial Conditions")
  {
    // Loop over the Initial Conditions.
    for (ParameterList::ConstIterator modelIt = myModels.begin();
      modelIt != myModels.end(); ++modelIt)
    {
      bool found(false);
      const string key(modelIt->first);
     
      const Teuchos::ParameterEntry& entry(modelIt->second);
      const ParameterList& plist(Teuchos::getValue<ParameterList>(entry));

      if (plist.isType<double>("Value"))
        found = createConstant(evaluators, defaults, key, plist.get<double>("Value"));

      TEUCHOS_TEST_FOR_EXCEPTION(not plist.isParameter("Value"),
        std::logic_error, "Error!  The model \"" + modelId +
        "\" with the key \"" + key + "\" does not have a valid \"Value\" " +
        "associated in the corresponding ParameterList!");

      if (plist.isType<string>("Value"))
      {
        const string value(plist.get<string>("Value"));

        // Check to see it is user data.                                  
        // JMG:  What needs to happen here?  12/16/2016.
        if (value == "User Data")
        {
          //TODO
        }

        //*********************************************************************
        // If it's not user data, then try direct key/value pairs.
        //*********************************************************************

        // Read from Exodus File.
        else if (value == "Exodus File")
          found = createGatherScaledFields(evaluators, defaults, key, userData);

        // Read "Input DOF Name" from Exodus File and map it to "DOF Name".
        else if (value == "Remap")
        {
          const string inputDOFName(plist.get<string>("Input DOF Name"));
          found = createGatherScaledFields(evaluators, defaults, inputDOFName,userData);
          found &= createICRemap(evaluators, defaults, key, inputDOFName);
        }

        // Support HB; create default values for all un-specified DOF frequencies 
        // TODO: support un-remapped key specification
        else if (value == "Frequency Domain Defaults" && m_fd_suffix == "_TP0_")
	{
          // gather DOF names *not* provided in the input deck via inclusion-...
          int maxHarmonic = plist.get<int>("Maximum harmonic");
          std::vector<string> keys(3*2*(maxHarmonic+1));
          for(int h = 0.0 ; h < maxHarmonic + 1 ; h++){
            keys[6*h+0] = ("ELECTRIC_POTENTIAL_CosH"+std::to_string(float(h))+"_");
            keys[6*h+1] = ("ELECTRIC_POTENTIAL_SinH"+std::to_string(float(h))+"_");
            keys[6*h+2] = ("ELECTRON_DENSITY_CosH"+std::to_string(float(h))+"_");
            keys[6*h+3] = ("ELECTRON_DENSITY_SinH"+std::to_string(float(h))+"_");
            keys[6*h+4] = ("HOLE_DENSITY_CosH"+std::to_string(float(h))+"_");
            keys[6*h+5] = ("HOLE_DENSITY_SinH"+std::to_string(float(h))+"_");
          }
          //for(auto key : keys)
          //  std::cout << "We might default IC for " << key << std::endl;

          // ... -exclusion
          std::vector<string> providedKeys;
          for (ParameterList::ConstIterator modelIt = myModels.begin();
               modelIt != myModels.end(); ++modelIt)
            if((modelIt->first).find("ELECTRIC_POTENTIAL_") != std::string::npos ||
               (modelIt->first).find("ELECTRON_DENSITY_")   != std::string::npos ||
               (modelIt->first).find("HOLE_DENSITY_")       != std::string::npos ){
              std::cout << "Found provided DOF for: " << modelIt->first << std::endl;
              providedKeys.push_back(modelIt->first);
            }
          for(auto providedKey : providedKeys)
            keys.erase(std::remove(keys.begin(), keys.end(), providedKey), keys.end()); 

          //for(auto key : keys)
          //  std::cout << "We default IC for " << key << std::endl;

          // take the default values to be 0.0, for now
          std::vector<double> values(keys.size(),0.0); // set them all to zero
          found = createFreqDomConstants(evaluators, defaults, keys, values);
	}
        else if (value == "Frequency Domain Defaults" && m_fd_suffix != "_TP0_")
          found = true;

        // Calculate initial carrier density from the equilibrium potential,
        // i.e., n0 = nie*exp(-Ei/kbT) and p0 = nie*exp(Ei/kbT).
        else if (value == "Equilibrium Density")
	{
	  found = createICEquilibriumDensity(evaluators, defaults, key, plist);
	}

        else if (value.substr(0,21) == "Equilibrium Potential")
	{
	  string SemBlockIds;
	  Teuchos::RCP<std::vector<string>> semBlocks = rcp(new vector<string>); 
	  if(value.substr(0,21).length() < value.length()) {
	    Teuchos::RCP<std::map<string,string>> block_mat = 
	      userData.get<Teuchos::RCP<std::map<string,string>>>("block2mat");
	    Material_Properties& matProperty = Material_Properties::getInstance();
	    const string matType(matProperty.getMaterialType(
			(*block_mat)[defaults.get<string>("Block ID")]));      
	    if (matType != "Insulator") {
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
					 "Error: SemBlockIds can only be specified in Insulators!");
	    }
	    size_t pos = value.find("SemBlockIds");
	    if (pos != string::npos) {
	      string tmp =  value.substr(pos,value.length());
	      size_t pos1 = tmp.find("=");
	      if (pos1 != string::npos) {
		SemBlockIds = tmp.substr(pos1+1,tmp.length());
	        // trim whitespaces
		size_t p = SemBlockIds.find_first_not_of(" \t");
		SemBlockIds.erase(0, p);
		p = SemBlockIds.find_last_not_of(" \t");
		if (string::npos != p)
		  SemBlockIds.erase(p+1);
		// strip off paranthesis
		SemBlockIds.erase(0,1);
		SemBlockIds.erase(SemBlockIds.length()-1,1); 
		stringstream ss(SemBlockIds); 
		while (ss.good()) { 
		  string substr; 
		  getline(ss, substr, ','); 
		  semBlocks->push_back(substr); 
		} 
		
		// validate blocks
		for(size_t kk=0; kk < semBlocks->size(); kk++) {
		  const string matType1(matProperty.getMaterialType(
				(*block_mat)[(*semBlocks)[kk]])); 
		  if (matType1 != "Semiconductor") {
		    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
			  "Error: " << (*semBlocks)[kk] << " must be a semiconductor block!");
		  }
		}
	      } else {
		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
					 "Error: SemBlockIds= must be specified!");
	      }
	    } else {
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
					 "Error: SemBlockIds must be specified!");
	    }
	    found = createICEquilibriumPotential(evaluators, defaults, key, userData, semBlocks);
	  } else {
            found = createICEquilibriumPotential(evaluators, defaults, key, userData, semBlocks);
	  }
	}
      
        // Build the IC_Gauss evaluator.
        else if (value == "Gauss")
          found = createICGauss(evaluators, defaults, key, plist);

        // Build the IC_Function evaluator.
        else if (value == "Function")
          found = createICFunction(evaluators, defaults, key, plist);
      }  // end of if (plist.isType<string>("Value"))
      if (not found)
      {
        stringstream msg;
        msg << "charon::ClosureModelFactory failed to build evaluator for "
            << "key \"" << key << "\"\n in model \"" << modelId << "\".  "
            << "Please correct the type or add support to the \nfactory."
            << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION(not found, std::logic_error, msg.str());
      }
    } // end of for loop
  }  // end of if (modelName == "Charon Parameter->Initial Conditions")


  // **************************************************************************
  // Throw exception for other modelName.
  // **************************************************************************
  else
  {
    if (m_throw_if_model_not_found)
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Invalid closure model name: " << modelName << "!" << std::endl);
  }
  return evaluators;
} // end of buildClosureModels()

///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateEmpiricalDefect()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateEmpiricalDefect(
  EvaluatorVector                evaluators,
  const Teuchos::ParameterList&  defaults,
  const panzer::IntegrationRule& ir,
  const Teuchos::ParameterList&  models,
  const Teuchos::ParameterList&  /* userData */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::PulseDamage_Spec;
  using charon::EmpiricalDamage_Data;
  using charon::RecombRate_Empirical_Defect;
  using panzer::Point;
  using panzer::BASIS;
  const string key("Empirical Defect Recombination");
  string param("eb x low");
  double ebxlo(-1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebxlo = models.sublist(key).get<double>(param);
  param = "eb x high";
  double ebxhi(1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebxhi = models.sublist(key).get<double>(param);
  param = "eb y low";
  double ebylo(-1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebylo = models.sublist(key).get<double>(param);
  param = "eb y high";
  double ebyhi(1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebyhi = models.sublist(key).get<double>(param);
  param = "eb z low";
  double ebzlo(-1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebzlo = models.sublist(key).get<double>(param);
  param = "eb z high";
  double ebzhi(1e3);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    ebzhi = models.sublist(key).get<double>(param);
  param = "thermal velocity";
  double thermalVelocity(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    thermalVelocity = models.sublist(key).get<double>(param);
  param = "cross section";
  double crossSection(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    crossSection = models.sublist(key).get<double>(param);
  param = "pulse data file";
  string pulseDataFile;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulseDataFile = models.sublist(key).get<string>(param);
  param = "mu data file";
  string muDataFile;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    muDataFile = models.sublist(key).get<string>(param);
  param = "pulse type";
  string pulseType;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulseType = models.sublist(key).get<string>(param);
  param = "pulse start";
  double pulseStart(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulseStart = models.sublist(key).get<double>(param);
  param = "pulse end";
  double pulseEnd(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulseEnd = models.sublist(key).get<double>(param);
  param = "pulse magnitude";
  double pulseMagnitude(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulseMagnitude = models.sublist(key).get<double>(param);
  param = "eb voltage override";
  double ebVoltageOverride(0);
  bool   ebOverrideBool(false);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
  {
    ebVoltageOverride = models.sublist(key).get<double>(param);
    ebOverrideBool    = true;
  }
  param = "cb voltage override";
  double cbVoltageOverride(0);
  bool cbOverrideBool(false);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
  {
    cbVoltageOverride = models.sublist(key).get<double>(param);
    cbOverrideBool    = true;
  }
  param = "pulse resolution";
  int pulsePoints(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pulsePoints = models.sublist(key).get<int>(param);
  param = "file pulse sampling scheme";
  string fDP;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    fDP = models.sublist(key).get<string>(param);
  param = "pulse is rate";
  bool pIR=false;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    pIR = models.sublist(key).get<bool>(param);
  RCP<BIRL> layout = defaults.get<RCP<BIRL>>("Basis");
  ParameterList in(key);
  in.set("eb x low",                   ebxlo                      );
  in.set("eb x high",                  ebxhi                      );
  in.set("eb y low",                   ebylo                      );
  in.set("eb y high",                  ebyhi                      );
  in.set("eb z low",                   ebzlo                      );
  in.set("eb z high",                  ebzhi                      );
  in.set("thermal velocity",           thermalVelocity          );
  in.set("cross section",              crossSection             );
  in.set("pulse data file",            pulseDataFile              );
  in.set("mu data file",               muDataFile               );
  in.set("pulse type",                 pulseType                );
  in.set("pulse start",                pulseStart               );
  in.set("pulse end",                  pulseEnd                 );
  in.set("pulse magnitude",            pulseMagnitude           );
  in.set("pulse resolution",           pulsePoints              );
  in.set("eb voltage override",        ebVoltageOverride          );
  in.set("eb voltage override bool",   ebOverrideBool             );
  in.set("cb voltage override",        cbVoltageOverride          );
  in.set("cb voltage override bool",   cbOverrideBool             );
  in.set("Names",                      names                      );
  in.set("IR",                         defaults.get<RCP<IR>>("IR"));
  in.set("Basis",                      layout                     );
  in.set("Scaling Parameters",         m_scale_params             );
  in.set("pulse is rate",              pIR                        );
  in.set("file pulse sampling scheme", fDP);
  in.set("damage spec object",
    models.sublist(key).get<RCP<PulseDamage_Spec> >("damage spec object"));
  in.set("empirical damage data",
         models.sublist(key).get<RCP<EmpiricalDamage_Data> >("empirical damage data"));

  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_Empirical_Defect<EvalT,
      Traits, Point>(*(layout->getBasis()), ir, in));
    evaluators->push_back(e);
  }

  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_Empirical_Defect<EvalT,
      Traits, BASIS>(*(layout->getBasis()), ir, in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRecombRateEmpiricalDefect()

///////////////////////////////////////////////////////////////////////////////
//
//  createConstant()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createConstant(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const double&                 value) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR);
  using panzer::Constant;
  ParameterList in;
  in.set("Name",  key  );
  in.set("Value", value);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createConstant()

///////////////////////////////////////////////////////////////////////////////
//
//  createFreqDomConstants()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createFreqDomConstants(
  EvaluatorVector                 evaluators,
  const Teuchos::ParameterList&   defaults,
  const std::vector<std::string>& keys,
  const std::vector<double>&      values) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR);
  using panzer::Constant;
  ParameterList in;

  for(unsigned int i = 0 ; i < keys.size() ; i++){
    in.set("Name",  keys[i]  );
    in.set("Value", values[i]);
    { // at IP
      in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
      RCP<Evaluator<Traits>> e = rcp(new Constant<EvalT, Traits>(in));
      evaluators->push_back(e);
    }
    { // at BASIS
      in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
      RCP<Evaluator<Traits>> e = rcp(new Constant<EvalT, Traits>(in));
      evaluators->push_back(e);
    }
  }
  return true;
} // end of createFreqDomConstants()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicConcentration()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicConcentration(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const double&                 value) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::IntrinsicConc_Default;
  ParameterList in;
  in.set("Value", value);
  in.set("Names", names);
  in.set("Scaling Parameters",m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicConcentration()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobility()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobility(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const double&                 value,
  const int&                    ionCharge /* = 0 */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::FEM_Velocity;
  using charon::Mobility_Default;
  ParameterList in;
  in.set("Value", value);
  in.set("Names", names);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      break;
    case CARRIER_TYPE_ION:
      in.set("Carrier Type", "Ion");
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge (center of an edge)
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    RCP<BIRL> basis = defaults.get<RCP<BIRL>>("Basis");
    in.set("Edge Data Layout", basis->getCellTopologyInfo()->edge_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  if (type == CARRIER_TYPE_ION)
  {
    // Instantiate the Ion Velocity evaluator together with the Ion
    // Mobility.
    ParameterList pl("Ion Velocity");
    pl.set<RCP<const Names>>("Names", names);
    pl.set("IR", defaults.get<RCP<IR>>("IR"));
    pl.set("Carrier Type", "Ion");
    pl.set<int>("Ion Charge", ionCharge);
    { // at IP
      RCP<Evaluator<Traits>> e = rcp(new FEM_Velocity<EvalT, Traits>(pl));
      evaluators->push_back(e);
    }
    { // at center of a primary edge
      pl.set<bool>("Is Edge Data Layout", true);
      pl.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
      RCP<Evaluator<Traits>> e = rcp(new FEM_Velocity<EvalT, Traits>(pl));
      evaluators->push_back(e);
    }
  }
  return true;
} // end of createMobility()

///////////////////////////////////////////////////////////////////////////////
//
//  createSRHLifetime()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createSRHLifetime(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const double&                 value) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::SRHLifetime_Constant;
  ParameterList in;
  in.set("Value", value);
  in.set("Names", names);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createSRHLifetime()

///////////////////////////////////////////////////////////////////////////////
//
//  createGatherFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createGatherFields(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const Teuchos::ParameterList& userData) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using panzer_stk::GatherFields;
  using panzer_stk::STK_Interface;
  ParameterList in;
  RCP<vector<string>> dofName = rcp(new vector<string>);
  dofName->push_back(key);
  in.set("Field Names", dofName                         );
  in.set("Basis",       defaults.get<RCP<BIRL>>("Basis"));
  RCP<const STK_Interface> mesh =
    userData.sublist("Panzer Data").get<RCP<STK_Interface>>("STK Mesh");
  RCP<Evaluator<Traits>> e = rcp(new GatherFields<EvalT, Traits>(mesh, in));
  evaluators->push_back(e);
  return true;
} // end of createGatherFields()

///////////////////////////////////////////////////////////////////////////////
//
//  createDOF()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDOF(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          key,
  const Teuchos::RCP<panzer::IntegrationRule> ir) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using panzer::DOF;
  ParameterList in;
  in.set("Name",  key                             );
  in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  in.set("IR",    ir                              );
  RCP<Evaluator<Traits>> e = rcp(new DOF<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createDOF()


////////////////////////////////////////////////////////////////////////////
//
//  createDOFGradient()
//
///////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDOFGradient(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using panzer::DOFGradient; 
  ParameterList in;

  if (cvfem_data.get<bool>("Is CVFEM"))
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis")); // HGrad basis, not HCurl basis
  }
  else
  {
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  in.set("Name",          key                             );
  in.set("Gradient Name", names->grad_dof.phi             );
  // in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  // in.set("IR",            ir                              );
  RCP<Evaluator<Traits>> e = rcp(new DOFGradient<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createDOFGradient()


//////////////////////////////////////////////////////////////////////////
//
//  createInitialPotentialGrad()
//
//////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createInitialPotentialGrad(
  EvaluatorVector                evaluators,
  const Teuchos::ParameterList&  defaults,
  const Teuchos::ParameterList&  userData,
  const Teuchos::ParameterList&  cvfem_data)  const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Initial_PotentialGrad;
  ParameterList in;

  if (cvfem_data.get<bool>("Is CVFEM"))
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis")); // HGrad basis, not HCurl basis
  }
  else
  {
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  in.set("Names", names                           );
  // in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  // in.set("IR",    ir                              );
  in.set("Max Worksets",userData.get<std::size_t>("Max Worksets"));
  RCP<Evaluator<Traits>> e = rcp(new Initial_PotentialGrad<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createInitialPotentialGrad()


////////////////////////////////////////////////////////////////////////////
//
//  createDopingStepJunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDopingStepJunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& plist,
  const bool&                   withIonizAcc,
  const bool&                   withIonizDon,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Doping_Function;
  using charon::Doping_StepJunction;
  ParameterList in;
  in.set("Names",             names                                 );
  in.set("IR",                defaults.get<RCP<IR>>("IR")           );
  in.set("Basis",             defaults.get<RCP<BIRL>>("Basis")      );
  in.set("Acceptor Value",    plist.get<double>("Acceptor Value")   );
  in.set("Donor Value",       plist.get<double>("Donor Value")      );
  in.set("Configuration",     plist.get<string>("Configuration")    );
  in.set("Direction",         plist.get<string>("Direction")        );
  in.set("Junction Location", plist.get<double>("Junction Location"));
  in.set("Scaling Parameters",m_scale_params);
  RCP<Evaluator<Traits>> e = rcp(new Doping_StepJunction<EvalT, Traits>(in));
  evaluators->push_back(e);
  if (withIonizAcc)
    in.sublist("IncmplIonizAcc Doping ParameterList") =
      models.sublist("Incomplete Ionized Acceptor").sublist("Model");
  if (withIonizDon)
    in.sublist("IncmplIonizDon Doping ParameterList") =
      models.sublist("Incomplete Ionized Donor").sublist("Model");
  e = rcp(new Doping_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createDopingStepJunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createDopingFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDopingFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const bool&                   withIonizAcc,
  const bool&                   withIonizDon,
  const Teuchos::RCP<panzer::GlobalData>& globalData,
  const Teuchos::ParameterList& userData,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::DopingRaw_Function;
  using charon::Doping_Function;
  ParameterList in;
  in.set("Names", names                           );
  in.set("IR",    defaults.get<RCP<IR>>("IR")     );
  in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  in.set("Scaling Parameters",m_scale_params      );
  in.set("Max Worksets",userData.get<std::size_t>("Max Worksets"));
  in.sublist("Doping ParameterList") = models.sublist("Doping");
  // doping parameterlist homotopy stuff

  if ((in.sublist("Doping ParameterList").isType<std::string>("Doping Homotopy")            )  and
      (in.sublist("Doping ParameterList").get<std::string>("Doping Homotopy") == "Parameter"))
    {
      in.sublist("Doping ParameterList").set<Teuchos::RCP<panzer::ParamLib>>(
        "ParamLib", globalData->pl);
    }

  if ((in.sublist("Doping ParameterList").isParameter("SweepingIsOn")))
    if (in.sublist("Doping ParameterList").get<bool>("SweepingIsOn") == true)
      {
	in.sublist("Doping ParameterList").set<Teuchos::RCP<panzer::ParamLib>>("ParamLib", globalData->pl);
      }



  if (withIonizAcc)
    in.sublist("IncmplIonizAcc Doping ParameterList") =
      models.sublist("Incomplete Ionized Acceptor").sublist("Model");
  if (withIonizDon)
    in.sublist("IncmplIonizDon Doping ParameterList") =
      models.sublist("Incomplete Ionized Donor").sublist("Model");

  RCP<Evaluator<Traits>> e = rcp(new DopingRaw_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  e = rcp(new Doping_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createDopingFunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createSpaceCharge()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createSpaceCharge(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Space_Charge;
  ParameterList in;
  in.set("Names",         names);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new Space_Charge<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new Space_Charge<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createSpaceCharge()



///////////////////////////////////////////////////////////////////////////////
//
//  createRelPermittivity()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRelPermittivity(
    EvaluatorVector               evaluators,
    const Teuchos::ParameterList& defaults,
    const std::string&            matName,
    const Teuchos::ParameterList& models) const 
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  
  ParameterList in;
  in.set("Names",            names         );
  in.set("Material Name",    matName       );
  in.sublist("Relative Permittivity ParameterList") = 
    models.sublist("Relative Permittivity");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Relative_Permittivity<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new Relative_Permittivity<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRelPermittivity()



///////////////////////////////////////////////////////////////////////////////
//
//  createPermittivityNitride()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createPermittivityNitride(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Permittivity_Nitride;
  ParameterList in;
  in.set("Names",            names         );
  in.set("Material Name",    matName       );
  // in.set("Scaling Parameters",m_scale_params);
  in.sublist("Relative Permittivity ParameterList") = 
    models.sublist("Relative Permittivity");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Permittivity_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new Permittivity_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createPermittivityNitride()

///////////////////////////////////////////////////////////////////////////////
//
//  createBandGapNitride()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createBandGapNitride(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::BandGap_Nitride;
  ParameterList in;
  in.set("Names",            names         );
  in.set("Material Name",    matName       );
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Bandgap ParameterList") = models.sublist("Band Gap");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new BandGap_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new BandGap_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createBandGapNitride()

///////////////////////////////////////////////////////////////////////////////
//
//  createBandGapTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createBandGapTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const bool&                   isAffinity,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::BandGap_TempDep;
  ParameterList in;
  in.set("Names",            names         );
  in.set("Material Name",    matName       );
  in.set("Compute Affinity", not isAffinity);
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Bandgap ParameterList") = models.sublist("Band Gap");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new BandGap_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new BandGap_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createBandGapTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicConcOldSlotboom()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicConcOldSlotboom(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            bgn,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::IntrinsicConc_OldSlotboom;
  ParameterList in;
  in.set("Names",              names  );
  in.set("Material Name",      matName);
  in.set("Band Gap Narrowing", bgn    );
  in.set("Scaling Parameters",m_scale_params);
  if (models.isSublist(names->field.intrin_conc))
    in.sublist("Intrinsic Conc ParameterList") =
      models.sublist(inp_names->field.intrin_conc);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_OldSlotboom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_OldSlotboom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicConcOldSlotboom()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicConcPersson()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicConcPersson(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            bgn,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::IntrinsicConc_Persson;
  ParameterList in;
  in.set("Names",              names  );
  in.set("Material Name",      matName);
  in.set("Band Gap Narrowing", bgn    );
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Intrinsic Conc ParameterList") =
    models.sublist(inp_names->field.intrin_conc);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Persson<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Persson<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicConcPersson()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicConcSlotboom()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicConcSlotboom(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            bgn,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::IntrinsicConc_Slotboom;
  ParameterList in;
  in.set("Names",              names  );
  in.set("Material Name",      matName);
  in.set("Band Gap Narrowing", bgn    );
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Intrinsic Conc ParameterList") =
    models.sublist(inp_names->field.intrin_conc);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Slotboom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Slotboom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicConcSlotboom()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicConcHarmon()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicConcHarmon(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            bgn,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::IntrinsicConc_Harmon;
  ParameterList in;
  in.set("Names",              names);
  in.set("Band Gap Narrowing", bgn  );
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Intrinsic Conc ParameterList") =
    models.sublist(inp_names->field.intrin_conc);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Harmon<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new IntrinsicConc_Harmon<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicConcHarmon()

///////////////////////////////////////////////////////////////////////////////
//
//  createEffectiveDOSSimple()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createEffectiveDOSSimple(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::EffectiveDOS_Simple;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  if (models.isSublist("Effective DOS"))
    in.sublist("Effective DOS ParameterList") =
      models.sublist("Effective DOS");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new EffectiveDOS_Simple<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new EffectiveDOS_Simple<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createEffectiveDOSSimple()

///////////////////////////////////////////////////////////////////////////////
//
//  createEffectiveDOSNitride()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createEffectiveDOSNitride(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::EffectiveDOS_Nitride;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  if (models.isSublist("Effective DOS"))
    in.sublist("Effective DOS ParameterList") =
      models.sublist("Effective DOS");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new EffectiveDOS_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new EffectiveDOS_Nitride<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createEffectiveDOSNitride()

///////////////////////////////////////////////////////////////////////////////
//
//  createSRHLifetimeFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createSRHLifetimeFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::SRHLifetime_Function;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      in.sublist("Lifetime ParameterList") =
        models.sublist(inp_names->field.elec_lifetime);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      in.sublist("Lifetime ParameterList") =
        models.sublist(inp_names->field.hole_lifetime);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Function<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Function<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createSRHLifetimeFunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityAnalytic()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityAnalytic(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Analytic;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Analytic<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Analytic<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge
    // information
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Analytic<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityAnalytic()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityArora()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityArora(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const bool useSuppliedParameterList) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Arora;
  ParameterList in;
 
  if(useSuppliedParameterList)
  {
    // the Arora mobility model does not take "IR" as an input parameter, so remove it.
    ParameterList bulkMobility;
    bulkMobility.setParameters(defaults);
    bulkMobility.remove("IR");
    in.setParameters(bulkMobility);
    in.setName("Arora");
  }
  else
  {
    in.set("Names",         names  );
    in.set("Material Name", matName);
    in.set("Scaling Parameters",m_scale_params);
  } 

  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      if(!useSuppliedParameterList)
        in.sublist("Mobility ParameterList") =
          models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      if(!useSuppliedParameterList)
        in.sublist("Mobility ParameterList") =
          models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Arora<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Arora<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge information
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    if (!useSuppliedParameterList) 
      in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Arora<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityArora()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityMasetti()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityMasetti(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Masetti;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Masetti<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Masetti<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge
    // information
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Masetti<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityMasetti()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityUniBo()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityUniBo(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_UniBo;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_UniBo<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_UniBo<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge
    // information
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_UniBo<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityUniBo()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityDopantTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityDopantTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const bool&                   bSolveIon,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_DopantTempDep;
  ParameterList in;
  if (bSolveIon)
    in.set("Dopant Name", names->dof.iondensity);
  else
    in.set("Dopant Name", names->field.donor);
  in.set("Names", names);
  in.set("Scaling Parameters",m_scale_params);
  in.sublist("Mobility ParameterList") =
    models.sublist(inp_names->field.elec_mobility);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new Mobility_DopantTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new Mobility_DopantTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityDopantTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityMOSFET()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityMOSFET(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_MOSFET;
  ParameterList in;
  ParameterList efin;
  in.set("Names",         names  );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName);
  in.set("Scaling Parameters",m_scale_params);

  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      efin.set("Carrier Type", "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      efin.set("Carrier Type", "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }

  //LCM : This is another thing that needs to change
  // create a copy of the parameter list as it may contain things that validate in the bulk or perpendicular models, but not in the MOSFET model
  ParameterList inCopy = in;
  in.sublist("Mobility ParameterList").remove("Bulk Mobility");
  in.sublist("Mobility ParameterList").remove("Perpendicular Field Model");
  in.sublist("Mobility ParameterList").set<std::string>("Bulk Mobility",inCopy.sublist("Mobility ParameterList").sublist("Bulk Mobility").get<std::string>("Value"));
  in.sublist("Mobility ParameterList").set<std::string>("Perpendicular Field Model",inCopy.sublist("Mobility ParameterList").sublist("Perpendicular Field Model").get<std::string>("Value"));
  { // at IP
    //in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));

    //Need to make sure we get the electric field for this evaluator when the discretization is SUPG
    efin.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    efin.set("IR", defaults.get<RCP<IR>>("IR"));
    efin.set("Scaling Parameters",m_scale_params);
    efin.set("Names",         names  );
    RCP<Evaluator<Traits> > ef = rcp(new FEM_ElectricField<EvalT, Traits>(efin));
    evaluators->push_back(ef);

    RCP<Evaluator<Traits>> e = rcp(new Mobility_MOSFET<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    //in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_MOSFET<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge information
    //in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_MOSFET<EvalT, Traits>(in));
    evaluators->push_back(e);
  }

  //The following blocks of code look bizarre.  But it's the only way to get what we want in the fewest lines of code.
  //Need to copy parameters out of the mosfet parameter list and create a parameter list that would be consistent
  //with bulk and Shirahata mobility models that are not created through the MOSFET mobility model.
  //Create the bulk mobility parameter list
  ParameterList bulkMobility;
  bulkMobility.setParameters(inCopy);
  bulkMobility.setName("MOSFET Bulk Mobility");
  bulkMobility.remove("Mobility ParameterList");
  if(inCopy.sublist("Mobility ParameterList").isSublist("Bulk Mobility"))
    bulkMobility.sublist("Mobility ParameterList").setParameters(inCopy.sublist("Mobility ParameterList").sublist("Bulk Mobility"));
  bulkMobility.sublist("Mobility ParameterList").set<std::string>("Compound Mobility", "On", "Rename the evaluated field to reflect compound mobility");
  if(inCopy.sublist("Mobility ParameterList").sublist("Bulk Mobility").isParameter("Enable High Field"))
  {
    bulkMobility.sublist("Mobility ParameterList").remove("Enable High Field");
    bulkMobility.sublist("Mobility ParameterList").set<std::string>("G Parameter Set", "Klaassen", "Use Klaassen or Meyer parameter set");
    bulkMobility.sublist("Mobility ParameterList").set<std::string>("Driving Force", "ElectricField", "Different high field calculation methods");
    bulkMobility.sublist("Mobility ParameterList").set<std::string>("High Field", "On", "Turn on/off high field dependence");
  }

  ParameterList perpMobility;
  perpMobility.setParameters(inCopy);
  perpMobility.setName("MOSFET Perpendicular Field Model");
  perpMobility.remove("Mobility ParameterList");
  if(inCopy.sublist("Mobility ParameterList").isSublist("Perpendicular Field Model"))
  {
    perpMobility.sublist("Mobility ParameterList").setParameters(inCopy.sublist("Mobility ParameterList").sublist("Perpendicular Field Model"));
    perpMobility.sublist("Mobility ParameterList").set("Value",inCopy.sublist("Mobility ParameterList").sublist("Perpendicular Field Model").get<std::string>("Value"));
  }
/*
  //If bulk is Klaassen and there is a perpendicular field model, disable lattice scattering in Klaassen
  if(perpMobility.isSublist("Mobility ParameterList"))
  {
    if(perpMobility.sublist("Mobility ParameterList").get<std::string>("Value") != "")
      bulkMobility.sublist("Mobility ParameterList").set<bool>("Disable Lattice Scattering",true);
  }
*/
  //Create the bulk mobility model
  if(bulkMobility.sublist("Mobility ParameterList").get<std::string>("Value") == "Klaassen")
  {
    // If bulk is Klaassen and there is a perpendicular field model, disable lattice scattering in Klaassen
    if(perpMobility.isSublist("Mobility ParameterList"))
    {
      if(perpMobility.sublist("Mobility ParameterList").get<std::string>("Value") != "")
        bulkMobility.sublist("Mobility ParameterList").set<bool>("Disable Lattice Scattering",true);
    }
    bulkMobility.sublist("Mobility ParameterList").set<std::string>("G Parameter Set", "Klaassen");
    createMobilityPhilipsThomas(evaluators, bulkMobility, type, matName, models, true);
  }
  
  if (bulkMobility.sublist("Mobility ParameterList").get<std::string>("Value") == "Arora")
    createMobilityArora(evaluators, bulkMobility, type, matName, models, true); 

  //Create the perp mobility model
  if(perpMobility.isSublist("Mobility ParameterList"))
  {
    if(perpMobility.sublist("Mobility ParameterList").get<std::string>("Value") == "Shirahata")
      createMobilityShirahata(evaluators, perpMobility, type, matName, models, true);
  }

  return true;
} // end of createMobilityMOSFET()


///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityShirahata()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityShirahata(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const bool useSuppliedParameterList) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Shirahata;
  ParameterList in;
  ParameterList efin;
  if(useSuppliedParameterList)
    {
      in.setParameters(defaults);
      in.setName("Shirahata Mobility");
    }
  else
    {
      in.set("Names",         names  );
      in.set("IR",            defaults.get<RCP<IR>>("IR")     );
      in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
      in.set("Material Name", matName);
      in.set("Scaling Parameters",m_scale_params);
    }
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      efin.set("Carrier Type", "Electron");
      if(!useSuppliedParameterList)
	in.sublist("Mobility ParameterList") =
	  models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      efin.set("Carrier Type", "Hole");
      if(!useSuppliedParameterList)
	in.sublist("Mobility ParameterList") =
	  models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }


  // create the Shirahata model here.
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Is Edge Data Layout", false);
    in.set("Nodal Evaluator", false);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));

    //Need to make sure we get the electric field for this evaluator when the discretization is NOT SUPG
    efin.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    efin.set("IR", defaults.get<RCP<IR>>("IR"));
    efin.set("Scaling Parameters",m_scale_params);
    efin.set("Names",         names  );
    RCP<Evaluator<Traits> > ef = rcp(new FEM_ElectricField<EvalT, Traits>(efin));
    evaluators->push_back(ef);

    RCP<Evaluator<Traits>> e = rcp(new Mobility_Shirahata<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    in.set("Nodal Evaluator", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Shirahata<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary element Edge, pass in Basis to retrieve edge
    // information
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Nodal Evaluator", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Shirahata<EvalT, Traits>(in));
    evaluators->push_back(e);
  }

  return true;
} // end of createMobilityShirahata()


///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityPhilipsThomas()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityPhilipsThomas(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const bool useSuppliedParameterList) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_PhilipsThomas;
  ParameterList in;
  ParameterList efin;
  if(useSuppliedParameterList)
    {
      in.setParameters(defaults);
      in.setName("Philips-Thomas");
    }
  else
    {
      in.set("Names",         names                           );
      in.set("IR",            defaults.get<RCP<IR>>("IR")     );
      in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
      in.set("Material Name", matName                         );
      in.set("Scaling Parameters",m_scale_params              );
    }
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type",  "Electron");
      efin.set("Carrier Type", "Electron");
      if(!useSuppliedParameterList)
	in.sublist("Mobility ParameterList") =
	  models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type",  "Hole");
      efin.set("Carrier Type", "Hole");
      if(!useSuppliedParameterList)
	in.sublist("Mobility ParameterList") =
	  models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }

  { // at IP, when isEdgedl = false
    in.set("Is Edge Data Layout", false);

    //Need to make sure we get the electric field for this evaluator when the discretization is NOT SUPG
    efin.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    efin.set("IR", defaults.get<RCP<IR>>("IR"));
    efin.set("Scaling Parameters",m_scale_params);
    efin.set("Names",         names  );
    RCP<Evaluator<Traits> > ef = rcp(new FEM_ElectricField<EvalT, Traits>(efin));
    evaluators->push_back(ef);

      RCP<Evaluator<Traits>> e =
      rcp(new Mobility_PhilipsThomas<EvalT, Traits>(in));
    evaluators->push_back(e);
  }

  { // at primary Edge, when isEdgedl = true
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e =
      rcp(new Mobility_PhilipsThomas<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityPhilipsThomas()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityLucent()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityLucent(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Lucent;
  ParameterList in;
  in.set("Names",         names                           );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName                         );
  in.set("Scaling Parameters",m_scale_params              );
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type",  "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type",  "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP, when isEdgedl = false
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Lucent<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary Edge, when isEdgedl = true
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Lucent<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityLucent()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityGaAs()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityGaAs(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_GaAs;
  ParameterList in;
  in.set("Names",         names                           );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName                         );
  in.set("Scaling Parameters",m_scale_params              );
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type",  "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type",  "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP, when isEdgedl = false
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_GaAs<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary Edge, when isEdgedl = true
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_GaAs<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityGaAs()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityAlbrecht()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityAlbrecht(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Albrecht;
  ParameterList in;
  in.set("Names",         names                           );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName                         );
  in.set("Scaling Parameters",m_scale_params              );
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type",  "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type",  "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP, when isEdgedl = false
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Albrecht<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary Edge, when isEdgedl = true
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Albrecht<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityAlbrecht()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityFarahmand()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityFarahmand(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_Farahmand;
  ParameterList in;
  in.set("Names",         names                           );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName                         );
  in.set("Scaling Parameters",m_scale_params              );
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type",  "Electron");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.elec_mobility);
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type",  "Hole");
      in.sublist("Mobility ParameterList") =
        models.sublist(inp_names->field.hole_mobility);
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  { // at IP, when isEdgedl = false
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Farahmand<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at primary Edge, when isEdgedl = true
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e = rcp(new Mobility_Farahmand<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityFarahmand()

///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheVanOverstraeten()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheVanOverstraeten(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& cvfem_data,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */ ) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::Avalanche_vanOverstraeten;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
     in.set("Scalar Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_scalar);
     in.set("Vector Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_vector);
  } else {
    in.set("Scalar Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Vector Data Layout", defaults.get<RCP<IR>>("IR")->dl_vector);
  }
  if (models.isSublist(names->field.avalanche_rate))
    in.sublist("Avalanche ParameterList") =
      models.sublist(names->field.avalanche_rate);
  else
    in.setName("Default Avalanche Generation");
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_vanOverstraeten<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheVanOverstraeten()

///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheOkuto()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheOkuto(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::Avalanche_Okuto;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
     in.set("Scalar Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_scalar);
     in.set("Vector Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_vector);
  } else {
    in.set("Scalar Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Vector Data Layout", defaults.get<RCP<IR>>("IR")->dl_vector);
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_Okuto<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheOkuto()

///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheLackner()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheLackner(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::Avalanche_Lackner;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
     in.set("Scalar Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_scalar);
     in.set("Vector Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_vector);
  } else {
    in.set("Scalar Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Vector Data Layout", defaults.get<RCP<IR>>("IR")->dl_vector);
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_Lackner<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheLackner()

///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheUniBo()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheUniBo(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::Avalanche_UniBo;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
     in.set("Scalar Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_scalar);
     in.set("Vector Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_vector);
  } else {
    in.set("Scalar Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Vector Data Layout", defaults.get<RCP<IR>>("IR")->dl_vector);
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_UniBo<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheUniBo()

///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheUniBoNew()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheUniBoNew(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::Avalanche_UniBoNew;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
     in.set("Scalar Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_scalar);
     in.set("Vector Data Layout", cvfem_data.get<RCP<IR>>("CVFEM Vol IR")->dl_vector);
  } else {
    in.set("Scalar Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Vector Data Layout", defaults.get<RCP<IR>>("IR")->dl_vector);
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_UniBoNew<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheUniBoNew()


///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheSelberherr()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheSelberherr(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            eqnSetType,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Avalanche_Selberherr;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Equation Set Type",  eqnSetType                            );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis")); // HGrad basis, not HCurl basis
  } else {
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis")); 
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_Selberherr<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheSelberherr()


///////////////////////////////////////////////////////////////////////////////
//
//  createAvalancheCrowellSze()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAvalancheCrowellSze(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            eqnSetType,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Avalanche_CrowellSze;
  ParameterList in;
  in.set("Names",              names                                 );
  in.set("Material Name",      matName                               );
  in.set("Equation Set Type",  eqnSetType                            );
  in.set("Scaling Parameters", m_scale_params                        );
  if (cvfem_data.get<bool>("Is CVFEM")) {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis")); // HGrad basis, not HCurl basis
  } else {
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis")); 
  }
  in.sublist("Avalanche ParameterList") =
    models.sublist(names->field.avalanche_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e =
      rcp(new Avalanche_CrowellSze<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAvalancheCrowellSze()


///////////////////////////////////////////////////////////////////////////////
//
//  createBand2BandTunnelingLocal()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createBand2BandTunnelingLocal(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const std::string&            eqnSetType,
  const Teuchos::ParameterList& models,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Band2Band_Tunneling_Local;
  ParameterList in;
  in.set("Names", names);
  in.set("Material Name", matName);
  in.set("Equation Set Type", eqnSetType);
  in.set("Scaling Parameters", m_scale_params);
  if (cvfem_data.get<bool>("Is CVFEM")) 
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis")); // HGrad basis, not HCurl basis
  }
  else 
  {
    in.set("IR", defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  in.sublist("Band2Band Tunneling ParameterList") =
    models.sublist(names->field.bbt_rate);
  { // at IP for now
    RCP<Evaluator<Traits>> e = rcp(new Band2Band_Tunneling_Local<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createBand2BandTunnelingLocal()


///////////////////////////////////////////////////////////////////////////////
//
//  createHeatCapacityTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createHeatCapacityTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::HeatCapacity_TempDep;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters", m_scale_params);
  ParameterList pl;
  if (models.isSublist(names->field.heat_cap))
    pl = models.sublist(names->field.heat_cap);
  else
  {
    pl.setName("Heat Capacity ParameterList");
    pl.set("Value", "TempDep");
  }
  in.sublist("Heat Capacity ParameterList") = pl;
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new HeatCapacity_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new HeatCapacity_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createHeatCapacityTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createHeatCapacityPowerLawTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createHeatCapacityPowerLawTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::HeatCapacity_PowerLawTempDep;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters", m_scale_params);
  ParameterList pl;
  if (models.isSublist(names->field.heat_cap))
    pl = models.sublist(names->field.heat_cap);
  else
  {
    pl.setName("Heat Capacity ParameterList");
    pl.set("Value", "PowerLawTempDep");
  }
  in.sublist("Heat Capacity ParameterList") = pl;
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new HeatCapacity_PowerLawTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new HeatCapacity_PowerLawTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createHeatCapacityPowerLawTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createThermalConductTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermalConductTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermalConduct_TempDep;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters", m_scale_params);
  ParameterList pl;
  if (models.isSublist(names->field.kappa))
    pl = models.sublist(names->field.kappa);
  else
  {
    pl.setName("Thermal Conductivity ParameterList");
    pl.set("Value", "TempDep");
  }
  in.sublist("Thermal Conductivity ParameterList") = pl;
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermalConductTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createThermalConductPowerLawTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermalConductPowerLawTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermalConduct_PowerLawTempDep;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters", m_scale_params);
  ParameterList pl;
  if (models.isSublist(names->field.kappa))
    pl = models.sublist(names->field.kappa);
  else
  {
    pl.setName("Thermal Conductivity ParameterList");
    pl.set("Value", "PowerLawTempDep");
  }
  in.sublist("Thermal Conductivity ParameterList") = pl;
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_PowerLawTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_PowerLawTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermalConductPowerLawTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createThermalConductLinearTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermalConductLinearTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermalConduct_LinearTempDep;
  ParameterList in;
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Thermal Conductivity ParameterList") =
    models.sublist(names->field.kappa);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_LinearTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_LinearTempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermalConductLinearTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createThermalConductLinearIonDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermalConductLinearIonDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermalConduct_LinearIonDep;
  ParameterList in;
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Thermal Conductivity ParameterList") =
    models.sublist(names->field.kappa);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_LinearIonDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermalConduct_LinearIonDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermalConductLinearIonDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createAnalyticHeatGeneration()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAnalyticHeatGeneration(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Analytic_HeatGeneration;
  ParameterList in;
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Heat Generation ParameterList") =
    models.sublist(names->field.heat_gen);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new Analytic_HeatGeneration<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new Analytic_HeatGeneration<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createAnalyticHeatGeneration()

///////////////////////////////////////////////////////////////////////////////
//
//  createSoretCoeffTempDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createSoretCoeffTempDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models /* = Teuchos::ParameterList() */) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::SoretCoeff_TempDep;
  ParameterList in;
  in.set("Names",         names  );
  in.set("Material Name", matName);
  in.set("Scaling Parameters", m_scale_params);
  ParameterList pl;
  if (models.isSublist(names->field.ion_soret_coeff))
    pl = models.sublist(names->field.ion_soret_coeff);
  else
  {
    pl.setName("Soret Coefficient ParameterList");
    pl.set("Value", "TempDep");
  }
  in.sublist("Soret Coefficient ParameterList") = pl;
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new SoretCoeff_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new SoretCoeff_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  {// at Edge (center of an edge)
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new SoretCoeff_TempDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createSoretCoeffTempDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createThermodiffCoeffCustom()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermodiffCoeffCustom(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermodiffCoeff_Custom;
  ParameterList in;
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Thermodiffusion Coefficient ParameterList") =
    models.sublist(names->field.ion_thermodiff_coeff);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Custom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Custom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  {// at Edge (center of an edge)
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Custom<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermodiffCoeffCustom()

///////////////////////////////////////////////////////////////////////////////
//
//  createMobilityRigidPointIon()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMobilityRigidPointIon(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const int&                    ionCharge,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Mobility_RigidPointIon;
  ParameterList in;
  in.set("Names",         names                           );
  in.set("IR",            defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",         defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName                         );
  in.set("Ion Charge",    ionCharge                       );
  in.set("Scaling Parameters", m_scale_params             );
  in.sublist("Mobility ParameterList") =
    models.sublist(names->field.ion_mobility);
  { // at IP
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e =
      rcp(new Mobility_RigidPointIon<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at center of a primary edge
    in.set("Is Edge Data Layout", true);
    RCP<Evaluator<Traits>> e =
      rcp(new Mobility_RigidPointIon<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createMobilityRigidPointIon()

///////////////////////////////////////////////////////////////////////////////
//
//  createDiffCoeffIonDep()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDiffCoeffIonDep(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::DiffCoeff_IonDep;
  ParameterList in;
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Diffusion ParameterList") =
    models.sublist(names->field.ion_diff_coeff);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new DiffCoeff_IonDep<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createDiffCoeffIonDep()

///////////////////////////////////////////////////////////////////////////////
//
//  createOptGenFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createOptGenFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::OptGen_Function;
  ParameterList in;
  in.set("Names", names                           );
  in.set("IR",    defaults.get<RCP<IR>>("IR")     );
  in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  in.set("Scaling Parameters", m_scale_params     );
  in.sublist("Optical Generation ParameterList") =
    models.sublist(names->field.opt_gen);
  RCP<Evaluator<Traits>> e = rcp(new OptGen_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createOptGenFunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createMoleFractionFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMoleFractionFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::MoleFraction_Function;
  ParameterList in;
  in.set("Names", names                           );
  in.set("IR",    defaults.get<RCP<IR>>("IR")     );
  in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  in.set("Material Name", matName);
  in.sublist("Mole Fraction ParameterList") =
    models.sublist(names->field.mole_frac);
  RCP<Evaluator<Traits>> e = rcp(new MoleFraction_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createMoleFractionFunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createBulkFixChargeFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createBulkFixChargeFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models,
  const Teuchos::RCP<panzer::GlobalData>& globalData,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::BulkFixCharge_Function;
  const string key(names->field.fixed_charge);
  ParameterList in(key);
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  in.sublist("Bulk FixCharge ParameterList") = models.sublist(key);
  ParameterList subsubList = in.sublist("Bulk FixCharge ParameterList").sublist("Function 1");
  if ((subsubList.isType<std::string>("Varying Charge Density")            )  and
      (subsubList.get<std::string>("Varying Charge Density") == "Parameter"))
    {
      in.sublist("Bulk FixCharge ParameterList").set<Teuchos::RCP<panzer::ParamLib>>(
        "ParamLib", globalData->pl);
    }
  if (cvfem_data.get<bool>("Is CVFEM"))  // SGCVFEM
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis"));  // HGrad basis, not HCurl basis (edge basis vectors)
  }
  else  // SUPG-FEM and EFFPG-FEM
  {
    in.set("IR",    defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  RCP<Evaluator<Traits>> e = rcp(new BulkFixCharge_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createBulkFixChargeFunction()

///////////////////////////////////////////////////////////////////////////////
//
//  createDiffCoeffDefault()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDiffCoeffDefault(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const bool&                   bUseFD,
  const std::string&            FDFormula) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::DiffCoeff_Default;
  ParameterList in;
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      break;
    case CARRIER_TYPE_ION:
      in.set("Carrier Type", "Ion");
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  in.set("Names",       names    );
  in.set("Fermi Dirac", bUseFD   );
  in.set("FD Formula",  FDFormula);
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new DiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e = rcp(new DiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at Edge (center of an edge)
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e = rcp(new DiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createDiffCoeffDefault()

///////////////////////////////////////////////////////////////////////////////
//
//  createSRHLifetimeConstant()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createSRHLifetimeConstant(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const CarrierType&            type,
  const double&                 value) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::SRHLifetime_Constant;
  ParameterList in;
  switch (type)
  {
    case CARRIER_TYPE_ELECTRON:
      in.set("Carrier Type", "Electron");
      break;
    case CARRIER_TYPE_HOLE:
      in.set("Carrier Type", "Hole");
      break;
    default:
      stringstream msg;
      msg << __PRETTY_FUNCTION__ << "was called with an invalid CarrierType: "
          << type;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, msg.str());
      return false;
  }
  in.set("Value", value);
  in.set("Names", names);
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new SRHLifetime_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createSRHLifetimeConstant()

///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateSRH()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateSRH(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const bool&                   bUseFD) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::RecombRate_SRH;
  ParameterList in(names->field.srh_recomb);
  in.set("Names",       names );
  in.set("Fermi Dirac", bUseFD);
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_SRH<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_SRH<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRecombRateSRH()


///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateDynamicTraps()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateDynamicTraps(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const std::string&            eqnSetType,
  const std::string&            drForce,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::RecombRate_DynamicTraps;
  const string key("Dynamic Traps Recombination");
  ParameterList in(key);
  in.set("Names",              names                 );
  in.set("Material Name",      matName               );
  in.set("Equation Set Type",  eqnSetType            );
  in.set("Driving Force",      drForce               );
  in.set("Scaling Parameters", m_scale_params        );
  if (cvfem_data.get<bool>("Is CVFEM"))  // SGCVFEM
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis"));  // HGrad basis, not HCurl basis (edge basis vectors)
  }
  else  // SUPG-FEM and EFFPG-FEM
  {
    in.set("IR",    defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  if (not models.isSublist(key))
  {
    stringstream msg;
    msg << "Error!  " << key
        << " ParameterList must be specified when Dynamic Traps = On!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }
  in.sublist("Dynamic Traps ParameterList") = models.sublist(key);
   
  RCP<Evaluator<Traits>> e = rcp(new RecombRate_DynamicTraps<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createRecombRateDynamicTraps()


///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateTrapSRH()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateTrapSRH(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const std::string&            eqnSetType,
  const std::string&            drForce,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::RecombRate_TrapSRH;
  const string key(names->field.trap_srh_recomb);
  ParameterList in(key);
  in.set("Names",              names                 );
  in.set("Material Name",      matName               );
  in.set("Equation Set Type",  eqnSetType            );
  in.set("Driving Force",      drForce               );
  in.set("Scaling Parameters", m_scale_params        );
  if (cvfem_data.get<bool>("Is CVFEM"))  // SGCVFEM
  {
    in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  // IR for subcv centroid
    in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis"));  // HGrad basis, not HCurl basis (edge basis vectors)
  }
  else  // SUPG-FEM and EFFPG-FEM
  {
    in.set("IR",    defaults.get<RCP<IR>>("IR"));
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
  }
  if (not models.isSublist(key))
  {
    stringstream msg;
    msg << "Error!  " << key
        << " ParameterList must be specified when Trap SRH = On!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }
  in.sublist("Trap SRH ParameterList") = models.sublist(key);
  RCP<Evaluator<Traits>> e = rcp(new RecombRate_TrapSRH<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createRecombRateTrapSRH()

///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateDefectCluster()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateDefectCluster(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::clusterInterpolator;
  using charon::RecombRate_Defect_Cluster;
  const string key("Defect Cluster Recombination");
  string param("Shepard Power");
  double shepardPneg(2.0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    shepardPneg = models.sublist(key).get<double>(param);
  param = "Cascade Density";
  double cascadeDensity(0.0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    cascadeDensity = models.sublist(key).get<double>(param);
  param = "Influence Radius";
  double influenceRadius(-1.0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    influenceRadius = models.sublist(key).get<double>(param);
  param = "Number of Input Files";
  int numberOfFiles(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    numberOfFiles = models.sublist(key).get<int>(param);
  param = "Interpolant Method";
  string methodName;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    methodName = models.sublist(key).get<string>(param);
  param = "Input File Type";
  string inputType;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    inputType = models.sublist(key).get<string>(param);
  param = "cluster interpolator";
  RCP<clusterInterpolator> cInterp;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    cInterp = models.sublist(key).get<RCP<clusterInterpolator>>(param);
  param = "Cluster Netlist Template";
  string templateName;
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    templateName = models.sublist(key).get<string>(param);
  ParameterList in(key);
  in.set("Cascade Density",       cascadeDensity                  );
  in.set("Influence Radius",      influenceRadius                 );
  in.set("Number of Input Files", numberOfFiles                   );
  in.set("Input File Type",       inputType                       );
  in.set("Interpolant Method",    methodName                      );
  in.set("Shepard Power",         shepardPneg                     );
  in.set("Names",                 names                           );
  in.set("IR",                    defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",                 defaults.get<RCP<BIRL>>("Basis"));
  in.set("Scaling Parameters", m_scale_params                     );
  in.set<RCP<clusterInterpolator>>("cluster interpolator", cInterp);
  in.set<string>("Cluster Netlist Template", templateName);
  { // at IP
    in.set("Is IP Set", true);
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new RecombRate_Defect_Cluster<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Is IP Set", false);
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new RecombRate_Defect_Cluster<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRecombRateDefectCluster()

///////////////////////////////////////////////////////////////////////////////
//
//  createIonizationParticleStrike()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIonizationParticleStrike(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Ionization_ParticleStrike;
  const string key("Particle Strike");
  string param("Start Point X");
  double startX(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    startX = models.sublist(key).get<double>(param);
  param = "Start Point Y";
  double startY(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    startY = models.sublist(key).get<double>(param);
  param = "Start Point Z";
  double startZ(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    startZ = models.sublist(key).get<double>(param);
  param = "End Point X";
  double endX(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    endX = models.sublist(key).get<double>(param);
  param = "End Point Y";
  double endY(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    endY = models.sublist(key).get<double>(param);
  param = "End Point Z";
  double endZ(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    endZ = models.sublist(key).get<double>(param);
  param = "Strike Radius";
  double strikeRadius(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    strikeRadius = models.sublist(key).get<double>(param);
  param = "Generation Rate";
  double generationRate(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    generationRate = models.sublist(key).get<double>(param);
  param = "Total Charge";
  double totalCharge(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    totalCharge = models.sublist(key).get<double>(param);
  param = "Start Time";
  double startTime(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    startTime = models.sublist(key).get<double>(param);
  param = "End Time";
  double endTime(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    endTime = models.sublist(key).get<double>(param);
  param = "Temporal Waveform";
  string TW("Square");
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    TW = models.sublist(key).get<string>(param);
  ParameterList in(key);
  in.set("Names",             names                           );
  in.set("IR",                defaults.get<RCP<IR>>("IR")     );
  in.set("Basis",             defaults.get<RCP<BIRL>>("Basis"));
  in.set("Start Point X",     startX                          );
  in.set("Start Point Y",     startY                          );
  in.set("Start Point Z",     startZ                          );
  in.set("End Point X",       endX                            );
  in.set("End Point Y",       endY                            );
  in.set("End Point Z",       endZ                            );
  in.set("Strike Radius",     strikeRadius                    );
  in.set("Generation Rate",   generationRate                  );
  in.set("Total Charge",      totalCharge                     );
  in.set("Start Time",        startTime                       );
  in.set("End Time",          endTime                         );
  in.set("Temporal Waveform", TW                              );
  in.set("Scaling Parameters", m_scale_params                 );
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new Ionization_ParticleStrike<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new Ionization_ParticleStrike<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIonizationParticleStrike()


///////////////////////////////////////////////////////////////////////////////
//
//  createTID()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createTID(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const Teuchos::ParameterList& models,
  const Teuchos::RCP<panzer::GlobalData>& globalData,
  const Teuchos::ParameterList& userData,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::KimptonTID;
  using Teuchos::Comm;

  const string key("TID");
 
  bool withKimptonModel = false;
  if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")) {
    withKimptonModel = true;
  }
  if (withKimptonModel) { // Kimpton model
    const string key1("Kimpton_TID");
    ParameterList in1(key1);

    string param("Dose");
    double dose(0);
    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isParameter(param) )
      dose = models.sublist(key).sublist("Kimpton Model").get<double>(param);
    in1.set("Dose",dose);

    param = "Effective Dose Enhancement Factor";
    double DEF(1);
    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isParameter(param) )
      DEF = models.sublist(key).sublist("Kimpton Model").get<double>(param);
    in1.set("Effective Dose Enhancement Factor",DEF);
    
    param = "Electron-Hole Pair Formation Energy";
    double Eform(0);
    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isParameter(param) ) 
      Eform = models.sublist(key).sublist("Kimpton Model").get<double>(param);
    in1.set("Electron-Hole Pair Formation Energy",Eform);

    param = "Electric Field Power Dependency";
    double pow_dep(0);
    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isParameter(param) )
      pow_dep = models.sublist(key).sublist("Kimpton Model").get<double>(param);
    in1.set("Electric Field Power Dependency",pow_dep);

    in1.set("Names", names);
    in1.set("Scaling Parameters", m_scale_params);
    in1.set("Material Name", matName);

    param = "Gate Contact";
    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isParameter(param) ) {
      // this parameter is checked in checkTID_param function in Charon_Main.cpp
      // 'Gate Contact' is a valid gate contact with varying voltage (LOCA)
      in1.set("ParamLib", globalData->pl);
      in1.set("Voltage Sweep",true);
      // freeze voltage
      param = "Freeze Voltage";
      if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	  and models.sublist(key).sublist("Kimpton Model").isParameter(param) ) {
	double V_freeze = models.sublist(key).sublist("Kimpton Model").get<double>(param);
	in1.set("Freeze Voltage", V_freeze);
      } else {
	stringstream msg;
	msg << "'Freeze Voltage' parameter in TID section is missing!" << std::endl;
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
      }


      }

    if (cvfem_data.get<bool>("Is CVFEM")) { // SGCVFEM
      // IR for subcv centroid
      in1.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));  
      // HGrad basis, not HCurl basis (edge basis vectors)
      in1.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis"));  
      in1.set("HCurlBasis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol HCurlBasis"));
      in1.set("Is CVFEM", true);
    } else { // SUPG-FEM and EFFPG-FEM
      in1.set("IR",    defaults.get<RCP<IR>>("IR"));
      in1.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
      in1.set("Is CVFEM", false);
    }

    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isSublist("Interface Traps") ) {
      in1.set("WithInterfaceTraps", true);
      // sideset name
      param = "Sideset ID";
      string sidesetID; 
      if (models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").isParameter(param)) {
	sidesetID = (models.sublist(key).sublist("Kimpton Model").
		     sublist("Interface Traps")).get<string>(param);
	in1.set("Sideset ID", sidesetID);
      } 
      // initial interface trap filling factor 
      double fill_facti(0);
      param = "Initial Filling Factor"; 
      if (models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").isParameter(param)) {
	fill_facti = (models.sublist(key).sublist("Kimpton Model").
		      sublist("Interface Traps")).get<double>(param);
	in1.set("Interface Initial Filling Factor", fill_facti);
      } 
      // interface trap capture cross section


      param = "Capture Cross Section Sweep is On";
      bool crossSectionSweep = false;
      if ( models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").isParameter(param))
	crossSectionSweep = models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").get<bool>(param);

      if (crossSectionSweep)
	{
	  bool ITTrapsSweep = (models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps")).get<bool>(param);
	  std::string ITCCSweep = "Interface Trap Capture Cross Section Sweep is On";
	  in1.set(ITCCSweep, ITTrapsSweep);
	  std::string initialCC = "Interface Trap Initial Capture Cross Section";
	  std::string finalCC = "Interface Trap Final Capture Cross Section";
	  double initialCrossSection = (models.sublist(key).sublist("Kimpton Model").
					sublist("Interface Traps")).get<double>("Initial Capture Cross Section");
	  double finalCrossSection = (models.sublist(key).sublist("Kimpton Model").
				      sublist("Interface Traps")).get<double>("Final Capture Cross Section");
	  in1.set(initialCC, initialCrossSection);
	  in1.set(finalCC, finalCrossSection);
	  in1.set("ParamLib", globalData->pl);
		  in1.set("Interface Trap Capture Cross Section Sweep","Parameter");
		  
	}
      else
	{
	  double sigma_i(0);
	  std::string param1 = "Capture Cross Section";
	  if (models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").isParameter(param1)) {
	    sigma_i = (models.sublist(key).sublist("Kimpton Model").
		       sublist("Interface Traps")).get<double>(param1);
	    in1.set("Interface Trap Capture Cross Section", sigma_i);
	  } 
	}
      
      // interface trap density 
      param = "Total Density Sweep is On";
      bool totalDensitySweep = false;
      if ( models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").isParameter(param))
	totalDensitySweep = models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps").get<bool>(param);

      if ( totalDensitySweep)
	{
	  bool ITDensitySweep = (models.sublist(key).sublist("Kimpton Model").sublist("Interface Traps")).get<bool>(param);
	  std::string ITTDSweep = "Interface Trap Total Density Sweep is On";
	  in1.set(ITTDSweep, ITDensitySweep);
	  std::string initialTD = "Interface Trap Initial Total Density";
	  std::string finalTD = "Interface Trap Final Total Density";
	  double initialTotalDensity = (models.sublist(key).sublist("Kimpton Model").
					sublist("Interface Traps")).get<double>("Initial Total Density");
	  double finalTotalDensity = (models.sublist(key).sublist("Kimpton Model").
				      sublist("Interface Traps")).get<double>("Final Total Density");
	  in1.set(initialTD, initialTotalDensity);
	  in1.set(finalTD, finalTotalDensity);
	  in1.set("ParamLib", globalData->pl);
	  in1.set("Interface Trap Total Density Sweep","Parameter");
	}
      else
	{
	  double Nti(0);
	  std::string param1 = "Total Density";
	      Nti = (models.sublist(key).sublist("Kimpton Model").
		     sublist("Interface Traps")).get<double>(param1);
	      in1.set("Interface Trap Density", Nti);
	}
      

    } // interface traps

    if (models.isSublist(key) and models.sublist(key).isSublist("Kimpton Model")
	and models.sublist(key).sublist("Kimpton Model").isSublist("Volume Traps") ) {
      in1.set("WithVolumeTraps", true);
      // volume trap density 
      double Ntv(0);
      param = "Total Density";
      if (models.sublist(key).sublist("Kimpton Model").sublist("Volume Traps").isParameter(param)) {
	Ntv = (models.sublist(key).sublist("Kimpton Model").
	       sublist("Volume Traps")).get<double>(param);
	in1.set("Volume Trap Density", Ntv);
      } 
      // initial volume trap filling factor 
      double fill_factv(0);
      param = "Initial Filling Factor"; 
      if (models.sublist(key).sublist("Kimpton Model").sublist("Volume Traps").isParameter(param)) {
	fill_factv = (models.sublist(key).sublist("Kimpton Model").
		      sublist("Volume Traps")).get<double>(param);
	in1.set("Volume Initial Filling Factor", fill_factv);
      } 
      // volume trap capture cross section
      double sigma_v(0);
      param = "Capture Cross Section";
      if (models.sublist(key).sublist("Kimpton Model").sublist("Volume Traps").isParameter(param)) {
	sigma_v = (models.sublist(key).sublist("Kimpton Model").
		   sublist("Volume Traps")).get<double>(param);
	in1.set("Volume Trap Capture Cross Section", sigma_v);
      }
      // volume trap critical capture cross section
      double sigma_v_crit(0);
      param = "Critical Capture Cross Section";
      if (models.sublist(key).sublist("Kimpton Model").sublist("Volume Traps").isParameter(param)) {
	sigma_v_crit = (models.sublist(key).sublist("Kimpton Model").
			sublist("Volume Traps")).get<double>(param);
	in1.set("Volume Trap Critical Capture Cross Section", sigma_v_crit);
      }
    } // volume traps

    std::string blockID = defaults.get<string>("Block ID");
    in1.set("Block ID",blockID);
    
    Teuchos::RCP<const panzer_stk::STK_Interface> mesh = 
      userData.sublist("Panzer Data").
      get<Teuchos::RCP<panzer_stk::STK_Interface>>("STK Mesh");
    in1.set("Mesh", mesh);
    
    in1.set("Comm",userData.get<RCP<const Comm<int> > >("Comm"));

    RCP<Evaluator<Traits>> e = rcp(new KimptonTID<EvalT, Traits>(in1));
    evaluators->push_back(e);
  }
  
  return true;
} // end of createTID()


///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateRadiative()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateRadiative(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const bool&                   bUseFD,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Material_Properties;
  using charon::RecombRate_Radiative;
  const string key(names->field.rad_recomb);
  string param("Coefficient");
  double recombCoeff(0.0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    recombCoeff = models.sublist(key).get<double>(param);
  else
  {
    Material_Properties& matProperty = Material_Properties::getInstance();
    recombCoeff = matProperty.getPropertyValue(matName,
      "Radiative Recombination Coefficient");
  }
  ParameterList in(key);
  in.set("Coefficient", recombCoeff);
  in.set("Names",       names      );
  in.set<bool>("Fermi Dirac", bUseFD);
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new RecombRate_Radiative<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new RecombRate_Radiative<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRecombRateRadiative()

///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateAuger()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateAuger(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            matName,
  const bool&                   bUseFD,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Material_Properties;
  using charon::RecombRate_Auger;
  Material_Properties& matProperty = Material_Properties::getInstance();
  const string key(names->field.auger_recomb);
  string param("Electron Auger Coefficient");
  double eAugerCoeff(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    eAugerCoeff = models.sublist(key).get<double>(param);
  else
    eAugerCoeff = matProperty.getPropertyValue(matName, param);
  param = "Hole Auger Coefficient";
  double hAugerCoeff(0);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    hAugerCoeff = models.sublist(key).get<double>(param);
  else
    hAugerCoeff = matProperty.getPropertyValue(matName, param);
  param = "With Generation";
  bool includeGen(false);
  if ((models.isSublist(key)) and (models.sublist(key).isParameter(param)))
    includeGen = models.sublist(key).get<bool>("With Generation");
  ParameterList in(key);
  in.set("Electron Auger Coefficient", eAugerCoeff);
  in.set("Hole Auger Coefficient",     hAugerCoeff);
  in.set("With Generation",            includeGen );
  in.set("Names",                      names      );
  in.set("Scaling Parameters", m_scale_params     );
  in.set<bool>("Fermi Dirac", bUseFD);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_Auger<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_Auger<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createRecombRateAuger()

///////////////////////////////////////////////////////////////////////////////
//
//  createRecombRateTotal()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createRecombRateTotal(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            srh,
  const std::string&            trapSrh,
  const std::string&            defectCluster,
  const std::string&            empiricalDefect,
  const std::string&            ionizationParticleStrike,
  const std::string&            rad,
  const std::string&            auger,
  const std::string&            optGen,
  const std::string&            avaGen,
  const std::string&            bbtGen,
  const std::string&            eqnSetType,
  const Teuchos::ParameterList& cvfem_data) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::RecombRate_Total;
  ParameterList in("Total Recombination");
  in.set("SRH",                 srh                     );
  in.set("Trap SRH",            trapSrh                 );
  in.set("Defect Cluster",      defectCluster           );
  in.set("Empirical Defect",    empiricalDefect         );
  in.set("Particle Strike",     ionizationParticleStrike);
  in.set("Radiative",           rad                     );
  in.set("Auger",               auger                   );
  in.set("Optical Generation",  optGen                  );
  in.set("Avalanche",           avaGen                  );
  in.set("Band2Band Tunneling", bbtGen                  );
  in.set("Equation Set Type",   eqnSetType              );
  in.set("Names",               names                   );
  { 
    // at finite element IP for SUPG-FEM and EFFPG-FEM, but at SubCV centroid for SGCVFEM formulations
    if (cvfem_data.get<bool>("Is CVFEM"))  // for SGCVFEM
    {
       in.set("IR", cvfem_data.get<RCP<IR>>("CVFEM Vol IR"));
       in.set("Basis", cvfem_data.get<RCP<BIRL>>("CVFEM Vol Basis"));
    } 
    else  // for SUPG-FEM and EFFPG-FEM
    {
      in.set("IR", defaults.get<RCP<IR>>("IR"));
      in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    }
    RCP<Evaluator<Traits>> e = rcp(new RecombRate_Total<EvalT, Traits>(in));
    evaluators->push_back(e);
  }

  return true;
} // end of createRecombRateTotal()

///////////////////////////////////////////////////////////////////////////////
//
//  createIntrinsicFermiEnergy()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createIntrinsicFermiEnergy(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Intrinsic_FermiEnergy;
  ParameterList in("Intrinsic Fermi Energy");
  in.set("Names",             names     );
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new Intrinsic_FermiEnergy<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new Intrinsic_FermiEnergy<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createIntrinsicFermiEnergy()

///////////////////////////////////////////////////////////////////////////////
//
//  createCondValeBand()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createCondValeBand(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::CondVale_Band;
  ParameterList in("Conduction and Valence Band");
  in.set("Names",             names     );
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new CondVale_Band<EvalT, Traits>(in));
    evaluators->push_back(e);
    e = rcp(new QuasiFermiLevels<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new CondVale_Band<EvalT, Traits>(in));
    evaluators->push_back(e);
    e = rcp(new QuasiFermiLevels<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createCondValeBand()

///////////////////////////////////////////////////////////////////////////////
//
//  createDegeneracyFactor()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createDegeneracyFactor(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const bool&                   bUseFD,
  const std::string&            FDFormula) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Degeneracy_Factor;
  ParameterList in("Degeneracy Factor");
  in.set("Names",       names    );
  in.set("Fermi Dirac", bUseFD   );
  in.set("FD Formula",  FDFormula);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Degeneracy_Factor<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new Degeneracy_Factor<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createDegeneracyFactor()

///////////////////////////////////////////////////////////////////////////////
//
//  createLatticeTempConstant()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createLatticeTempConstant(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const Teuchos::ParameterList& models) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::LatticeTemp_Constant;
  using charon::Material_Properties;
  Material_Properties& matProperty = Material_Properties::getInstance();
  //const std::string key(names->field.latt_temp);
  const std::string key("Lattice Temperature");
  double propertyValue(0);
  if (models.isParameter(key))
  {
    propertyValue = models.get<double>(key);
    matProperty.setPropertyValue(key, propertyValue);
  }
  else
    propertyValue = matProperty.getPropertyValue(key);
  ParameterList in;
  in.set("Names", names        );
  in.set("Value", propertyValue);
  in.set("Scaling Parameters", m_scale_params);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e =
      rcp(new LatticeTemp_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e =
      rcp(new LatticeTemp_Constant<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createLatticeTempConstant()

///////////////////////////////////////////////////////////////////////////////
//
//  createReferenceEnergy()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createReferenceEnergy(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            refMaterial,
  const Teuchos::ParameterList& submodels) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::Reference_Energy;
  ParameterList in("Reference Energy");
  in.set("Names",              names      );
  in.set("Reference Material", refMaterial);
  in.set("Scaling Parameters", m_scale_params);
  if (submodels.isSublist("Electron Affinity"))
  {
    const ParameterList& chiParamList = submodels.sublist("Electron Affinity");
    if (chiParamList.isType<double>("Value"))
      in.set("Constant Electron Affinity", chiParamList.get<double>("Value"));
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "User-defined Electron Affinity must be a constant!");
  }
  if (submodels.isSublist("Band Gap"))
  {
    const ParameterList& bgParamList = submodels.sublist("Band Gap");
    if (bgParamList.isType<double>("Value"))
      in.set("Constant Band Gap", bgParamList.get<double>("Value"));
    else if (bgParamList.isType<string>("Value"))
      in.sublist("Bandgap ParameterList") = bgParamList;
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Wrong type of Band Gap->Value!");
  }
  if (submodels.isSublist("Effective DOS"))
    in.sublist("Effective DOS ParameterList") =
      submodels.sublist("Effective DOS");
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    RCP<Evaluator<Traits>> e = rcp(new Reference_Energy<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    RCP<Evaluator<Traits>> e = rcp(new Reference_Energy<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createReferenceEnergy()

///////////////////////////////////////////////////////////////////////////////
//
//  createGlobalStatistics()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createGlobalStatistics(
  EvaluatorVector                             evaluators,
  Teuchos::ParameterList&                     in,
  const std::string&                          value,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::RCP<panzer::GlobalData>      globalData,
  const Teuchos::ParameterList&               userData,
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_NONE);
  using panzer::GlobalStatistics;
  using Teuchos::Comm;
  in.set("Names",       value                                     );
  in.set("IR",          ir                                        );
  in.set("Global Data", globalData                                );
  in.set("Comm",        userData.get<RCP<const Comm<int>>>("Comm"));
  RCP<GlobalStatistics<EvalT, Traits>> e =
    rcp(new GlobalStatistics<EvalT, Traits>(in));
  evaluators->push_back(e);
  fm.template requireField<EvalT>(e->getRequiredFieldTag());
  return true;
} // end of createGlobalStatistics()

///////////////////////////////////////////////////////////////////////////////
//
//  createMMSAnalyticSolution()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createMMSAnalyticSolution(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const panzer::FieldLayoutLibrary&           fl,
  const std::string&                          modelId) const
{
  CHARON_USINGS_ETC(DEFINE_NAMES);
  using charon::MMS_NLP_GLH_1_AnalyticSolution;
  using charon::MMS_DD_RDH_1_AnalyticSolution;
  using charon::MMS_DD_RDH_2_AnalyticSolution;
  using panzer::FieldLayoutLibrary;
  RCP<const FieldLayoutLibrary> fieldLayoutLibrary = Teuchos::rcpFromRef(fl);
  const string prefix("Analytic_");
  bool found(false);
  ParameterList pl;
  if (boost::iequals(value, "mms_nlp_glh_1"))
  {
    pl.set("Scaling Parameters", m_scale_params);
    RCP<Evaluator<Traits>> e = rcp(new MMS_NLP_GLH_1_AnalyticSolution<EvalT,
      Traits>(prefix, *names, fieldLayoutLibrary, ir, pl));
    evaluators->push_back(e);
    found = true;
  }
  else if (boost::iequals(value, "mms_dd_rdh_1"))
  {
    pl.set("Scaling Parameters", m_scale_params);
    RCP<Evaluator<Traits>> e = rcp(new MMS_DD_RDH_1_AnalyticSolution<EvalT,
      Traits>(prefix, *names, fieldLayoutLibrary, ir, pl));
    evaluators->push_back(e);
    found = true;
  }
  else if (boost::iequals(value, "mms_dd_rdh_2"))
  {
    pl.set("Scaling Parameters", m_scale_params);
    RCP<Evaluator<Traits>> e = rcp(new MMS_DD_RDH_2_AnalyticSolution<EvalT,
      Traits>(prefix, *names, fieldLayoutLibrary, ir, pl));
    evaluators->push_back(e);
    found = true;
  }
  if (not found)
  {
    std::stringstream err_msg;
    err_msg << "ClosureModelFactory failed to build evaluator for analytic "
            << "solution \"" << value << "\" in model \"" << modelId
            << "\"." << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, err_msg.str());
  }
  return found;
} // end of createMMSAnalyticSolution()

///////////////////////////////////////////////////////////////////////////////
//
//  createAnalyticComparison()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAnalyticComparison(
  EvaluatorVector                    evaluators,
  const Teuchos::ParameterList&      defaults,
  const std::string&                 value,
  Teuchos::ParameterList&            in,
  PHX::FieldManager<panzer::Traits>& fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using charon::AnalyticComparison;
  using PHX::FieldTag;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "Error_"                                    );
    RCP<Evaluator<Traits>> e = rcp(new AnalyticComparison<EvalT, Traits>(in));
    evaluators->push_back(e);
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));
  }
  return true;
} // end of createAnalyticComparison()

///////////////////////////////////////////////////////////////////////////////
//
//  createAnalyticComparisonL2Error()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAnalyticComparisonL2Error(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using charon::AnalyticComparison_L2Error;
  using PHX::FieldTag;
  using Teuchos::Comm;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    // in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Basis",           defaults.get<RCP<BIRL>>("Basis")            );
    in.set("IR",              ir                                          );
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "L2Error_"                                  );
    in.set("Comm",            userData.get<RCP<const Comm<int>>>("Comm")  );
    RCP<Evaluator<Traits>> e =
      rcp(new AnalyticComparison_L2Error<EvalT, Traits>(in));
    evaluators->push_back(e);
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));
    
    
  }
  return true;
} // end of createAnalyticComparisonL2Error()

///////////////////////////////////////////////////////////////////////////////
//
//  createAnalyticComparisonRelError()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createAnalyticComparisonRelError(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using charon::AnalyticComparison_RelError;
  using PHX::FieldTag;
  using Teuchos::Comm;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  if (not in.isParameter("Use Absolute"))
    in.set<bool>("Use Absolute", true);
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    in.set("DataLayout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Name", *fieldName);
    in.set("Analytic Prefix", "Analytic_");
    in.set("Error Prefix", "RelError_");
    RCP<Evaluator<Traits>> e =
      rcp(new AnalyticComparison_RelError<EvalT, Traits>(in));
    evaluators->push_back(e);
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));
  }
  return true;
} // end of createAnalyticComparisonRelError()


///////////////////////////////////////////////////////////////////////////////
//
//  createManufacturedSolution()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createManufacturedSolution(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData, //NEW ADDITION
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_NONE);
  using boost::char_separator;
  using boost::tokenizer;
  using PHX::FieldTag;
  using Teuchos::Comm;


  using charon::DD_RDH_1_AnalyticSolution;
  bool found(false);
  if (boost::iequals(value, "mms_dd_rdh_1"))
  {
    RCP<Evaluator<Traits>> e = rcp(new DD_RDH_1_AnalyticSolution<EvalT,
      Traits>(ir));
    evaluators->push_back(e);
    found = true;
  }
 
  if (not found)
  {
    std::stringstream err_msg;
    err_msg << "ClosureModelFactory failed to build evaluator for analytic "
            << "solution \"" << value << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, err_msg.str());
  }
  return found;

} // end of createManufacturedSolution()


///////////////////////////////////////////////////////////////////////////////
//
//  createNormCalculationL2Error()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createNormCalculationL2Error(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData, //NEW ADDITION
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using PHX::FieldTag;
  using Teuchos::Comm;
  using charon::Norm_L2Error;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    // in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Basis",           defaults.get<RCP<BIRL>>("Basis")            );
    in.set("IR",              ir                                          );
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "L2Error_"                                  );
    in.set("Comm",            userData.get<RCP<const Comm<int>>>("Comm")  );
    
    // For norm scalar exo output
    Teuchos::ParameterList NORM;
    in.sublist("Norm ParameterList") = NORM;
    in.sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib>>("ParamLib", globalData->pl);

    
    // NEW ADDITIONS FOR Norm_L2Error
    RCP<Evaluator<Traits>> e = rcp(new Norm_L2Error<EvalT, Traits>(in));
    evaluators->push_back(e);
    
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));

        
    // Just throwing in an analytic solution to test error calculation
    //RCP<Evaluator<Traits>> e2 = rcp(new DD_RDH_1_AnalyticSolution<EvalT, Traits>("analytic_solution", ir));
    //evaluators->push_back(e2);
    
  }
  return true;
} // end of createNormCalculationL2Error()


///////////////////////////////////////////////////////////////////////////////
//
//  createNormCalculationH1Error()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createNormCalculationH1Error(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData, //NEW ADDITION
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using PHX::FieldTag;
  using Teuchos::Comm;
  using charon::Norm_H1Error;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    // in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Basis",           defaults.get<RCP<BIRL>>("Basis")            );
    in.set("IR",              ir                                          );
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "H1Error_"                                  );
    in.set("Comm",            userData.get<RCP<const Comm<int>>>("Comm")  );
    
    // For norm scalar exo output
    Teuchos::ParameterList NORM;
    in.sublist("Norm ParameterList") = NORM;
    in.sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib>>("ParamLib", globalData->pl);

    
    // NEW ADDITIONS FOR Norm_L2Error
    RCP<Evaluator<Traits>> e = rcp(new Norm_H1Error<EvalT, Traits>(in));
    evaluators->push_back(e);
    
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));

        
    // Just throwing in an analytic solution to test error calculation
    //RCP<Evaluator<Traits>> e2 = rcp(new AnalyticSolution<EvalT, Traits>("analytic_solution", ir));
    //RCP<Evaluator<Traits>> e2 = rcp(new DD_RDH_1_AnalyticSolution<EvalT, Traits>("analytic_solution", ir));
    //evaluators->push_back(e2);
    
  }
  return true;
} // end of createNormCalculationH1Error()


///////////////////////////////////////////////////////////////////////////////
//
//  createNormCalculationL2()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createNormCalculationL2(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData, //NEW ADDITION
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using PHX::FieldTag;
  using Teuchos::Comm;
  using charon::Norm_L2;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    // in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Basis",           defaults.get<RCP<BIRL>>("Basis")            );
    in.set("IR",              ir                                          );
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "L2Norm_"                                  );
    in.set("Comm",            userData.get<RCP<const Comm<int>>>("Comm")  );
    
    // For norm scalar exo output
    Teuchos::ParameterList NORM;
    in.sublist("Norm ParameterList") = NORM;
    in.sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib>>("ParamLib", globalData->pl);

    
    // NEW ADDITIONS FOR Norm_L2
    RCP<Evaluator<Traits>> e = rcp(new Norm_L2<EvalT, Traits>(in));
    evaluators->push_back(e);
    
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));

    
  }
  return true;
} // end of createNormCalculationL2()


///////////////////////////////////////////////////////////////////////////////
//
//  createNormCalculationH1()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createNormCalculationH1(
  EvaluatorVector                             evaluators,
  const Teuchos::ParameterList&               defaults,
  const std::string&                          value,
  Teuchos::ParameterList&                     in,
  const Teuchos::RCP<panzer::IntegrationRule> ir,
  const Teuchos::ParameterList&               userData,
  const Teuchos::RCP<panzer::GlobalData>&     globalData, //NEW ADDITION
  PHX::FieldManager<panzer::Traits>&          fm) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using boost::char_separator;
  using boost::tokenizer;
  using PHX::FieldTag;
  using Teuchos::Comm;
  using charon::Norm_H1;
  typedef tokenizer<char_separator<char>>::const_iterator FieldNameIterator;
  char_separator<char> sep(": , ");
  tokenizer<char_separator<char>> fieldNames(value, sep);
  for (FieldNameIterator fieldName = fieldNames.begin();
    fieldName != fieldNames.end(); ++fieldName)
  {
    // in.set("DataLayout",      defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Basis",           defaults.get<RCP<BIRL>>("Basis")            );
    in.set("IR",              ir                                          );
    in.set("Name",            *fieldName                                  );
    in.set("Analytic Prefix", "Analytic_"                                 );
    in.set("Error Prefix",    "H1Norm_"                                  );
    in.set("Comm",            userData.get<RCP<const Comm<int>>>("Comm")  );
    
    // For norm scalar exo output
    Teuchos::ParameterList NORM;
    in.sublist("Norm ParameterList") = NORM;
    in.sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib>>("ParamLib", globalData->pl);

    
    // NEW ADDITIONS FOR Norm_H1
    RCP<Evaluator<Traits>> e = rcp(new Norm_H1<EvalT, Traits>(in));
    evaluators->push_back(e);
    
    vector<RCP<FieldTag>> fieldTags = e->evaluatedFields();
    for (size_t fieldTag(0); fieldTag < fieldTags.size(); ++fieldTag)
      fm.template requireField<EvalT>(*(fieldTags[fieldTag]));

    
  }
  return true;
} // end of createNormCalculationH1()


///////////////////////////////////////////////////////////////////////////////
//
//  createThermodiffCoeffDefault()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createThermodiffCoeffDefault(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_IR_NAMES);
  using charon::ThermodiffCoeff_Default;
  ParameterList in;
  in.set("Names", names);
  { // at IP
    in.set("Data Layout", defaults.get<RCP<IR>>("IR")->dl_scalar);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at BASIS
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", false);
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  { // at Edge (center of an edge)
    in.set("Data Layout", defaults.get<RCP<BIRL>>("Basis")->functional);
    in.set("Is Edge Data Layout", true);
    in.set("Basis", defaults.get<RCP<BIRL>>("Basis"));
    RCP<Evaluator<Traits>> e =
      rcp(new ThermodiffCoeff_Default<EvalT, Traits>(in));
    evaluators->push_back(e);
  }
  return true;
} // end of createThermodiffCoeffDefault()

///////////////////////////////////////////////////////////////////////////////
//
//  createFEMGradNegPotential()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createFEMGradNegPotential(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults) const
{
  CHARON_USINGS_ETC(DEFINE_IR_NAMES);
  using charon::FEM_GradNegPotential;
  ParameterList in("Negative Potential Gradient");
  in.set("Names", names                      );
  in.set("IR",    defaults.get<RCP<IR>>("IR"));
  RCP<Evaluator<Traits>> e = rcp(new FEM_GradNegPotential<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createFEMGradNegPotential()

///////////////////////////////////////////////////////////////////////////////
//
//  createGatherScaledFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createGatherScaledFields(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            inputDOFName,
  const Teuchos::ParameterList& userData) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using charon::GatherScaledFields;
  using panzer_stk::STK_Interface;
  RCP<vector<string>> dofName = rcp(new vector<string>);
  dofName->push_back(inputDOFName);
  ParameterList in;
  in.set("Field Names",        dofName                         );
  in.set("Basis",              defaults.get<RCP<BIRL>>("Basis"));
  in.set("Scaling Parameters", m_scale_params                  );
  RCP<const STK_Interface> mesh =
    userData.sublist("Panzer Data").get<RCP<STK_Interface>>("STK Mesh");
  RCP<Evaluator<Traits>> e =
    rcp(new GatherScaledFields<EvalT, Traits>(mesh, in));
  evaluators->push_back(e);
  return true;
} // end of createGatherScaledFields()

///////////////////////////////////////////////////////////////////////////////
//
//  createICRemap()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createICRemap(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const std::string&            inputDOFName) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_NAMES);
  using charon::IC_Remap;
  ParameterList in;
  in.set("DOF Name",       key );
  in.set("Input DOF Name", inputDOFName);
  in.set("Basis",          defaults.get<RCP<BIRL>>("Basis"));
  in.set("Names",          names);
  RCP<Evaluator<Traits>> e = rcp(new IC_Remap<EvalT, Traits>(in));
  evaluators->push_back(e);

  return true;
} // end of createICRemap()

///////////////////////////////////////////////////////////////////////////////
//
//  createICEquilibriumDensity()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createICEquilibriumDensity(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const Teuchos::ParameterList& plist) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_NAMES);
  using charon::IC_Equilibrium_Density;
  ParameterList in;
  in.set("DOF Name", key                             );
  in.set("Basis",    defaults.get<RCP<BIRL>>("Basis"));
  in.set("Names",    names                           );
  in.set("Scaling Parameters",m_scale_params         );
  in.sublist("Equilibrium ParameterList") = plist;
  RCP<Evaluator<Traits>> e =
    rcp(new IC_Equilibrium_Density<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createICEquilibriumDensity()

///////////////////////////////////////////////////////////////////////////////
//
//  createICEquilibriumPotential()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createICEquilibriumPotential(
  EvaluatorVector                 evaluators,
  const Teuchos::ParameterList&   defaults,
  const std::string&              key,
  const Teuchos::ParameterList&   userData,
  const Teuchos::RCP<std::vector<std::string>> semBlocks) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL_NAMES);
  using charon::IC_Equilibrium_Potential;

  ParameterList in;
  in.set("DOF Name",          key                              );
  in.set("Basis",             defaults.get<RCP<BIRL>>("Basis") );
  in.set("Names",             names                            );
  in.set("Scaling Parameters",m_scale_params                   );
  in.set("Equation Set Type", defaults.get<std::string>("Type"));
  const ParameterList& iesParams = defaults.sublist("Options");
  string isFermiDirac("False");
  if (iesParams.isParameter("Fermi Dirac"))
    isFermiDirac = iesParams.get<string>("Fermi Dirac");
  bool bUseFD(false);
  if (isFermiDirac == "True")
    bUseFD = true;
  in.set("Fermi Dirac", bUseFD);
  Teuchos::RCP<const panzer_stk::STK_Interface> mesh = 
    userData.sublist("Panzer Data").
    get<Teuchos::RCP<panzer_stk::STK_Interface>>("STK Mesh");
  in.set("Mesh", mesh);
  Teuchos::RCP<std::vector<string>> sc_cnts = 
    userData.get<Teuchos::RCP<std::vector<string>>>("SchottkyContacts");
  in.set("SchottkyContacts",sc_cnts);
  Teuchos::RCP<std::vector<string>> sc_blks = 
    userData.get<Teuchos::RCP<std::vector<string>>>("SchottkyBlocks");
  in.set("SchottkyBlocks",sc_blks);
  Teuchos::RCP<std::vector<double>> sc_wf = 
    userData.get<Teuchos::RCP<std::vector<double>>>("SchottkyWF");
  in.set("SchottkyWF",sc_wf);
  Teuchos::RCP<std::vector<double>> sc_vapp = 
    userData.get<Teuchos::RCP<std::vector<double>>>("SchottkyVapp");
  in.set("SchottkyVapp",sc_vapp);
  Teuchos::RCP<std::vector<string>> g_cnts = 
    userData.get<Teuchos::RCP<std::vector<string>>>("GateContacts");
  in.set("GateContacts",g_cnts);
  Teuchos::RCP<std::vector<string>> g_blks = 
    userData.get<Teuchos::RCP<std::vector<string>>>("GateBlocks");
  in.set("GateBlocks",g_blks);
  Teuchos::RCP<std::vector<double>> g_wf = 
    userData.get<Teuchos::RCP<std::vector<double>>>("GateWF");
  in.set("GateWF",g_wf);
  Teuchos::RCP<std::vector<double>> g_vapp = 
    userData.get<Teuchos::RCP<std::vector<double>>>("GateVapp");
  in.set("GateVapp",g_vapp);

  Teuchos::RCP<std::map<string,string>> block_mat = 
    userData.get<Teuchos::RCP<std::map<string,string>>>("block2mat");
  in.set("Material",(*block_mat)[defaults.get<string>("Block ID")]);

  in.set("Block ID", defaults.get<string>("Block ID"));

  in.set("Semiconductor Blocks",semBlocks);
  Teuchos::RCP<std::vector<string>> semMaterials = 
    rcp(new std::vector<string>); 
  if(semBlocks->size() > 0) {
    for (size_t kk=0; kk<semBlocks->size(); kk++) {
      semMaterials->push_back((*block_mat)[(*semBlocks)[kk]]);
    }
  } 
  in.set("Semiconductor Materials",semMaterials);
 
  RCP<Evaluator<Traits>> e =
    rcp(new IC_Equilibrium_Potential<EvalT, Traits>(in));
  evaluators->push_back(e);

  return true;
} // end of createICEquilibriumPotential()

///////////////////////////////////////////////////////////////////////////////
//
//  createICGauss()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createICGauss(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const Teuchos::ParameterList& plist) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using charon::IC_Gauss;
  ParameterList in;
  in.set("DOF Name", key                             );
  in.set("Basis",    defaults.get<RCP<BIRL>>("Basis"));
  in.sublist("Gauss ParameterList") = plist;
  RCP<Evaluator<Traits>> e = rcp(new IC_Gauss<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createICGauss()

///////////////////////////////////////////////////////////////////////////////
//
//  createICFunction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
bool charon::ClosureModelFactory<EvalT>::
createICFunction(
  EvaluatorVector               evaluators,
  const Teuchos::ParameterList& defaults,
  const std::string&            key,
  const Teuchos::ParameterList& plist) const
{
  CHARON_USINGS_ETC(DEFINE_BIRL);
  using charon::IC_Function;
  ParameterList in;
  in.set("DOF Name", key                             );
  in.set("Basis",    defaults.get<RCP<BIRL>>("Basis"));
  in.sublist("Function ParameterList") = plist;
  RCP<Evaluator<Traits>> e = rcp(new IC_Function<EvalT, Traits>(in));
  evaluators->push_back(e);
  return true;
} // end of createICFunction()


///////////////////////////////////////////////////////////////////////////////
//
//  setupMoleFraction()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT>
void charon::ClosureModelFactory<EvalT>::
setupMoleFraction(const Teuchos::ParameterList& myModels) const {
  using std::stringstream;
  using Teuchos::ParameterList;

  Material_Properties& matProperty = Material_Properties::getInstance();
  const string& matName(myModels.get<string>("Material Name"));
  if (matProperty.hasMoleFracDependence(matName)) {
    if (!myModels.isSublist("Mole Fraction")) {
      stringstream msg;
      msg << "'Mole Fraction' must be specified for mole-fraction"
	  << "\ndependent material '" << matName << "'" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }

    // check if mole fractions are specified correctly according to arity
    const ParameterList& moleFracParamList = myModels.sublist("Mole Fraction");
    for (ParameterList::ConstIterator model_it = moleFracParamList.begin();
	 model_it != moleFracParamList.end(); ++model_it) {
      const string key = model_it->first;
      if (key.compare(0, 8, "Function") == 0) {
	const Teuchos::ParameterEntry& entry = model_it->second;
	const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
	const string funcType = funcParamList.get<string>("Function Type");
	if (matProperty.getArityType(matName) == "Binary" or
	    matProperty.getArityType(matName) == "Ternary") {
	  if (funcType == "Uniform") {
	    if (!funcParamList.isParameter("xMoleFrac Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	  } else if (funcType == "Linear") {
	    if (!funcParamList.isParameter("xMoleFrac Start Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac Start Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	    if (!funcParamList.isParameter("xMoleFrac End Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac End Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	  }
	} else if (matProperty.getArityType(matName) == "Quaternary") {
	  if (funcType == "Uniform") {
	    if (!funcParamList.isParameter("xMoleFrac Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	    if (!funcParamList.isParameter("yMoleFrac Value")) {
	      stringstream msg;
	      msg << "'yMoleFrac Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	  } else if (funcType == "Linear") {
	    if (!funcParamList.isParameter("xMoleFrac Start Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac Start Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	    if (!funcParamList.isParameter("xMoleFrac End Value")) {
	      stringstream msg;
	      msg << "'xMoleFrac End Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	    if (!funcParamList.isParameter("yMoleFrac Start Value")) {
	      stringstream msg;
	      msg << "'yMoleFrac Start Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	    if (!funcParamList.isParameter("yMoleFrac End Value")) {
	      stringstream msg;
	      msg << "'yMoleFrac End Value' must be specified for mole-fraction"
		  << "\ndependent material '" << matName << "'" << std::endl;
	      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
	    }
	  }
	}
      }
    } // arity error check

    // on demand mole-fraction material setup
    matProperty.registerMoleFracMaterial(matName);

    // bandgap mole fraction parameters setup
    if (myModels.isSublist("Band Gap")) {
      bool withTempDepBG = false;
      const ParameterList& bgParamList(myModels.sublist("Band Gap"));
      if (bgParamList.isType<string>("Value") and
	  bgParamList.get<string>("Value") == "TempDep")
	withTempDepBG = true;
      if (withTempDepBG) {
	double Eg300_b = 0.0, Eg300_c = 0.0;
	double Eg300_0 = 0.0, Eg300_1 = 0.0;
	double Chi300_b = 0.0, Chi300_c = 0.0;
	double Chi300_0 = 0.0, Chi300_1 = 0.0;
	double alpha_b = 0.0, alpha_c = 0.0;
	double alpha_0 = 0.0, alpha_1 = 0.0;
	double beta_b = 0.0, beta_c = 0.0;
	double beta_0 = 0.0, beta_1 = 0.0;
	// retrieve user input if specified
	if (bgParamList.isSublist("Mole Fraction Parameters")) {
	  const ParameterList& mf_params = bgParamList.sublist("Mole Fraction Parameters");
	  if (mf_params.isSublist("Eg300")) {
	    const ParameterList& Eg300_params = mf_params.sublist("Eg300");
	    if (Eg300_params.isParameter("b"))
	      Eg300_b = Eg300_params.get<double>("b");
	    if (Eg300_params.isParameter("c"))
	      Eg300_c = Eg300_params.get<double>("c");  
	  }
	  if (mf_params.isSublist("Chi300")) {
	    const ParameterList& Chi300_params = mf_params.sublist("Chi300");
	    if (Chi300_params.isParameter("b"))
	      Chi300_b = Chi300_params.get<double>("b");
	    if (Chi300_params.isParameter("c"))
	      Chi300_c = Chi300_params.get<double>("c");
	  }
	  if (mf_params.isSublist("alpha")) {
	    const ParameterList& alpha_params = mf_params.sublist("alpha");
	    if (alpha_params.isParameter("b"))
	      alpha_b = alpha_params.get<double>("b");
	    if (alpha_params.isParameter("c"))
	      alpha_c = alpha_params.get<double>("c");
	  }
	  if (mf_params.isSublist("beta")) {
	    const ParameterList& beta_params = mf_params.sublist("beta");
	    if (beta_params.isParameter("b"))
	      beta_b = beta_params.get<double>("b");
	    if (beta_params.isParameter("c"))
	      beta_c = beta_params.get<double>("c");
	  }
	  if (mf_params.isParameter("Eg300(x=0)")) {
	    Eg300_0 = mf_params.get<double>("Eg300(x=0)");
	  }
	  if (mf_params.isParameter("Eg300(x=1)")) {
	    Eg300_1 = mf_params.get<double>("Eg300(x=1)");
	  }
	  if (mf_params.isParameter("Chi300(x=0)")) {
	    Chi300_0 = mf_params.get<double>("Chi300(x=0)");
	  }
	  if (mf_params.isParameter("Chi300(x=1)")) {
	    Chi300_1 = mf_params.get<double>("Chi300(x=1)");
	  }
	  if (mf_params.isParameter("alpha(x=0)")) {
	    alpha_0 = mf_params.get<double>("alpha(x=0)");
	  }
	  if (mf_params.isParameter("alpha(x=1)")) {
	    alpha_1 = mf_params.get<double>("alpha(x=1)");
	  }
	  if (mf_params.isParameter("beta(x=0)")) {
	    beta_0 = mf_params.get<double>("beta(x=0)");
	  }
	  if (mf_params.isParameter("beta(x=1)")) {
	    beta_1 = mf_params.get<double>("beta(x=1)");
	  }
	}
	matProperty.setupMoleFracBandgapParams(matName,Eg300_b, Eg300_c,Chi300_b,Chi300_c,
					       alpha_b,alpha_c,beta_b,beta_c, Eg300_0, Eg300_1,
					       Chi300_0, Chi300_1, alpha_0, alpha_1, beta_0, beta_1);
      } // withTempDepBG
    } // bandgap mole fraction parameters setup

    // relative permittivity mole fraction parameters setup
    if (myModels.isSublist("Relative Permittivity")) {
      const ParameterList& permParamList(myModels.sublist("Relative Permittivity"));
      double Value_b = 0.0, Value_c = 0.0;
      double Value_0 = 0.0, Value_1 = 0.0;
      // retrieve user input if specified
      if (permParamList.isSublist("Mole Fraction Parameters")) {
	const ParameterList& mf_params = permParamList.sublist("Mole Fraction Parameters");
	if (mf_params.isSublist("Value")) {
	  const ParameterList& Value_params = mf_params.sublist("Value");
	  if (Value_params.isParameter("b"))
	    Value_b = Value_params.get<double>("b");
	  if (Value_params.isParameter("c"))
	    Value_c = Value_params.get<double>("c");  
	  if (mf_params.isParameter("Value(x=0)")) {
	    Value_0 = mf_params.get<double>("Value(x=0)");
	  }
	  if (mf_params.isParameter("Value(x=1)")) {
	    Value_1 = mf_params.get<double>("Value(x=1)");
	  }
	}
      }
      matProperty.setupMoleFracRelPermittivityParams(matName,Value_b, Value_c, 
						     Value_0, Value_1);
    } // relative permittivity mole fraction parameters setup

    // Mihai
    // effective DOS mole fraction parameters setup
    if (myModels.isSublist("Effective DOS")) {
      const ParameterList& effDOSParamList(myModels.sublist("Effective DOS"));
      if (effDOSParamList.isParameter("Value") and 
	  effDOSParamList.isType<string>("Value") and 
	  effDOSParamList.get<string>("Value") == "Simple") {
	double Nc300_b = 0.0, Nc300_c = 0.0;
	double Nc300_0 = 0.0, Nc300_1 = 0.0;
        double Nv300_b = 0.0, Nv300_c = 0.0;
	double Nv300_0 = 0.0, Nv300_1 = 0.0;
        double Nc_F_b = 0.0, Nc_F_c = 0.0;
        double Nc_F_0 = 0.0, Nc_F_1 = 0.0;
        double Nv_F_b = 0.0, Nv_F_c = 0.0;
        double Nv_F_0 = 0.0, Nv_F_1 = 0.0;
	// retrieve user input if specified
	if (effDOSParamList.isSublist("Mole Fraction Parameters")) {
	  const ParameterList& mf_params = effDOSParamList.sublist("Mole Fraction Parameters");
	  if (mf_params.isSublist("Nc300")) {
	    const ParameterList& Nc300_params = mf_params.sublist("Nc300");
	    if (Nc300_params.isParameter("b"))
	      Nc300_b = Nc300_params.get<double>("b");
	    if (Nc300_params.isParameter("c"))
	      Nc300_c = Nc300_params.get<double>("c");  
	  }
	  if (mf_params.isSublist("Nc_F")) {
	    const ParameterList& Nc_F_params = mf_params.sublist("Nc_F");
	    if (Nc_F_params.isParameter("b"))
	      Nc_F_b = Nc_F_params.get<double>("b");
	    if (Nc_F_params.isParameter("c"))
	      Nc_F_c = Nc_F_params.get<double>("c");  
	  }
          if (mf_params.isSublist("Nv300")) {
	    const ParameterList& Nv300_params = mf_params.sublist("Nv300");
	    if (Nv300_params.isParameter("b"))
	      Nv300_b = Nv300_params.get<double>("b");
	    if (Nv300_params.isParameter("c"))
	      Nv300_c = Nv300_params.get<double>("c");  
	  }
          if (mf_params.isSublist("Nv_F")) {
	    const ParameterList& Nv_F_params = mf_params.sublist("Nv_F");
	    if (Nv_F_params.isParameter("b"))
	      Nv_F_b = Nv_F_params.get<double>("b");
	    if (Nv_F_params.isParameter("c"))
	      Nv_F_c = Nv_F_params.get<double>("c");  
	  }
	  if (mf_params.isParameter("Nc300(x=0)")) {
	    Nc300_0 = mf_params.get<double>("Nc300(x=0)");
	  }
	  if (mf_params.isParameter("Nc300(x=1)")) {
	    Nc300_1 = mf_params.get<double>("Nc300(x=1)");
	  }
	  if (mf_params.isParameter("Nc_F(x=0)")) {
	    Nc_F_0 = mf_params.get<double>("Nc_F(x=0)");
	  }
	  if (mf_params.isParameter("Nc_F(x=1)")) {
	    Nc_F_1 = mf_params.get<double>("Nc_F(x=1)");
	  }
          if (mf_params.isParameter("Nv300(x=0)")) {
	    Nv300_0 = mf_params.get<double>("Nv300(x=0)");
	  }
	  if (mf_params.isParameter("Nv300(x=1)")) {
	    Nv300_1 = mf_params.get<double>("Nv300(x=1)");
	  }
	  if (mf_params.isParameter("Nv_F(x=0)")) {
	    Nv_F_0 = mf_params.get<double>("Nv_F(x=0)");
	  }
	  if (mf_params.isParameter("Nv_F(x=1)")) {
	    Nv_F_1 = mf_params.get<double>("Nv_F(x=1)");
	  }
	}
	matProperty.setupMoleFracDOSParams(matName,Nc300_b,Nc300_c,Nv300_b,Nv300_c,
					   Nc_F_b,Nc_F_c,Nv_F_b,Nv_F_c,Nc300_0,Nc300_1,
					   Nv300_0,Nv300_1,Nc_F_0,Nc_F_1,Nv_F_0,Nv_F_1);
      }
    } // effective DOS mole fraction parameters setup
    // Mihai


  } // mole fraction material
} // end of setupMoleFraction()



#endif // CHARON_CLOSURE_MODEL_FACTORY_IMPL_HPP
