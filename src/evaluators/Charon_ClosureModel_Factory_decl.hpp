
#ifndef   CHARON_CLOSURE_MODEL_FACTORY_DECL_HPP
#define   CHARON_CLOSURE_MODEL_FACTORY_DECL_HPP

//#include "Charon_Mobility_Analytic.hpp"
#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

#include "Panzer_ClosureModel_Factory.hpp"

namespace panzer
{
  class InputEquationSet;
}

namespace charon
{
  template<typename EvalT>
  class ClosureModelFactory
    :
    public panzer::ClosureModelFactory<EvalT>
  {
    public:
      enum Define
      {
        DEFINE_NONE,
        DEFINE_BIRL,
        DEFINE_IR,
        DEFINE_BIRL_IR,
        DEFINE_NAMES,
        DEFINE_BIRL_NAMES,
        DEFINE_IR_NAMES,
        DEFINE_BIRL_IR_NAMES
      };
      enum CarrierType
      {
        CARRIER_TYPE_UNDEFINED,
        CARRIER_TYPE_ELECTRON,
        CARRIER_TYPE_HOLE,
        CARRIER_TYPE_ION
      };
      typedef
        Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
        EvaluatorVector;

      ClosureModelFactory(
        const Teuchos::RCP<charon::Scaling_Parameters> & scaleParams,
        bool throwIfModelNotFound = true,
        const std::string& modelNamePrefix = "Charon Parameters->",
        const std::string fd_suffix = "")
        :
        m_scale_params(scaleParams),
        m_throw_if_model_not_found(throwIfModelNotFound),
        m_model_name_prefix(modelNamePrefix),
        m_fd_suffix(fd_suffix)
      {
      } // end of Constructor

      EvaluatorVector
      buildClosureModels(
        const std::string&                           model_id,
        const Teuchos::ParameterList&                models,
        const panzer::FieldLayoutLibrary&            fl,
        const Teuchos::RCP<panzer::IntegrationRule>& ir,
        const Teuchos::ParameterList&                default_params,
        const Teuchos::ParameterList&                user_data,
        const Teuchos::RCP<panzer::GlobalData>&      global_data,
        PHX::FieldManager<panzer::Traits>&           fm) const;

    private:
      Teuchos::RCP<charon::Scaling_Parameters> m_scale_params;
      bool m_throw_if_model_not_found;
      std::string m_model_name_prefix;
      std::string m_fd_suffix;

      bool createRecombRateEmpiricalDefect(
        EvaluatorVector                evaluators,
        const Teuchos::ParameterList&  defaults,
        const panzer::IntegrationRule& ir,
        const Teuchos::ParameterList&  models,
        const Teuchos::ParameterList&  userData) const;
      bool createConstant(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const double&                 value) const;
      bool createFreqDomConstants(
        EvaluatorVector                 evaluators,
        const Teuchos::ParameterList&   defaults,
        const std::vector<std::string>& keys,
        const std::vector<double>&      values) const;
      bool createIntrinsicConcentration(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const double&                 value) const;
      bool createMobility(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const double&                 value,
        const int&                    ionCharge = 0) const;
      bool createSRHLifetime(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const double&                 value) const;
      bool createGatherFields(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const Teuchos::ParameterList& userData) const;
      bool createDOF(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          key,
        const Teuchos::RCP<panzer::IntegrationRule> ir) const;
      bool createDOFGradient(
        EvaluatorVector                evaluators,
        const Teuchos::ParameterList&  defaults,
        const std::string&             key,
        const Teuchos::ParameterList&  cvfem_data) const;
      bool createDopingStepJunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& plist,
        const bool&                   withIonizAcc,
        const bool&                   withIonizDon,
        const Teuchos::ParameterList& models) const;
      bool createDopingFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const bool&                   withIonizAcc,
        const bool&                   withIonizDon,
        const Teuchos::RCP<panzer::GlobalData>& globalData,
        const Teuchos::ParameterList& user_data,
        const Teuchos::ParameterList& models) const;
      bool createRelPermittivity(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createPermittivityNitride(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createBandGapNitride(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createBandGapTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const bool&                   isAffinity,
        const Teuchos::ParameterList& models) const;
      bool createIntrinsicConcOldSlotboom(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const std::string&            bgn,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createIntrinsicConcPersson(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const std::string&            bgn,
        const Teuchos::ParameterList& models) const;
      bool createIntrinsicConcSlotboom(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const std::string&            bgn,
        const Teuchos::ParameterList& models) const;
      bool createIntrinsicConcHarmon(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            bgn,
        const Teuchos::ParameterList& models) const;
      bool createEffectiveDOSSimple(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createEffectiveDOSNitride(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createSRHLifetimeFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityAnalytic(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityArora(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const bool useSuppliedParameterList) const;
      bool createMobilityMasetti(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityUniBo(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityMOSFET(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityKlaassen(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
	const bool usesuppliedParamterList) const;
      bool createMobilityShirahata(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
	const bool usesuppliedParamterList) const;
      bool createMobilityDopantTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const bool&                   bSolveIon,
        const Teuchos::ParameterList& models) const;
      bool createMobilityPhilipsThomas(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
	const bool useSuppliedParameterList) const;
      bool createMobilityLucent(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityGaAs(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityAlbrecht(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createMobilityFarahmand(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createAvalancheVanOverstraeten(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& cvfem_data,
        const Teuchos::ParameterList& models = Teuchos::ParameterList() ) const;
      bool createAvalancheOkuto(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createAvalancheLackner(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createAvalancheUniBo(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createAvalancheUniBoNew(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createAvalancheSelberherr(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const std::string&            eqnSetType,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createAvalancheCrowellSze(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const std::string&            eqnSetType,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createBand2BandTunnelingLocal(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string& matName,
        const std::string& eqnSetType,
        const Teuchos::ParameterList& models,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createHeatCapacityTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createHeatCapacityPowerLawTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createThermalConductTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createThermalConductPowerLawTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createThermalConductLinearTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createThermalConductLinearIonDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createAnalyticHeatGeneration(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createSoretCoeffTempDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models = Teuchos::ParameterList()) const;
      bool createThermodiffCoeffCustom(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createMobilityRigidPointIon(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const int&                    ionCharge,
        const Teuchos::ParameterList& models) const;
      bool createDiffCoeffIonDep(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createOptGenFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createMoleFractionFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models) const;
      bool createBulkFixChargeFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models,
	const Teuchos::RCP<panzer::GlobalData>& globalData,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createDiffCoeffDefault(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const bool&                   bUseFD,
        const std::string&            FDFormula) const;
      bool createSRHLifetimeConstant(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const CarrierType&            type,
        const double&                 value) const;
      bool createRecombRateSRH(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const bool&                   bUseFD) const;
      bool createRecombRateTrapSRH(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const std::string&            eqnSetType,
        const std::string&            drForce,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createRecombRateDynamicTraps(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const Teuchos::ParameterList& models,
        const std::string&            eqnSetType,
        const std::string&            drForce,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createRecombRateDefectCluster(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createIonizationParticleStrike(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createTID(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
	const std::string&            matName,
        const Teuchos::ParameterList& models,
	const Teuchos::RCP<panzer::GlobalData>& globalData,
	const Teuchos::ParameterList& userData,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createRecombRateRadiative(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const bool&                   bUseFD,
        const Teuchos::ParameterList& models) const;
      bool createRecombRateAuger(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            matName,
        const bool&                   bUseFD,
        const Teuchos::ParameterList& models) const;
      bool createRecombRateTotal(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            srh,
        const std::string&            trap_srh,
        const std::string&            defect_cluster,
        const std::string&            empirical_defect,
        const std::string&            ionization_particle_strike,
        const std::string&            rad,
        const std::string&            auger,
        const std::string&            opt_gen,
        const std::string&            ava_gen,
        const std::string&            bbtGen,
        const std::string&            eqnSetType,
        const Teuchos::ParameterList& cvfem_data) const;
      bool createIntrinsicFermiEnergy(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults) const;
      bool createCondValeBand(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults) const;
      bool createDegeneracyFactor(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const bool&                   bUseFD,
        const std::string&            FDFormula) const;
      bool createLatticeTempConstant(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& models) const;
      bool createReferenceEnergy(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            refMaterial,
        const Teuchos::ParameterList& submodels) const;
      bool createGlobalStatistics(
        EvaluatorVector                             evaluators,
        Teuchos::ParameterList&                     in,
        const std::string&                          value,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::RCP<panzer::GlobalData>      globalData,
        const Teuchos::ParameterList&               userData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createMMSAnalyticSolution(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const panzer::FieldLayoutLibrary&           fl,
        const std::string&                          modelId) const;
      bool createAnalyticComparison(
        EvaluatorVector                    evaluators,
        const Teuchos::ParameterList&      defaults,
        const std::string&                 value,
        Teuchos::ParameterList&            in,
        PHX::FieldManager<panzer::Traits>& fm) const;
      bool createAnalyticComparisonL2Error(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createAnalyticComparisonRelError(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createManufacturedSolution(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        const Teuchos::RCP<panzer::GlobalData>&      globalData,
        PHX::FieldManager<panzer::Traits>&          fm) const;        
      bool createNormCalculationL2Error(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        const Teuchos::RCP<panzer::GlobalData>&      globalData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createNormCalculationH1Error(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        const Teuchos::RCP<panzer::GlobalData>&      globalData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createNormCalculationL2(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        const Teuchos::RCP<panzer::GlobalData>&      globalData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createNormCalculationH1(
        EvaluatorVector                             evaluators,
        const Teuchos::ParameterList&               defaults,
        const std::string&                          value,
        Teuchos::ParameterList&                     in,
        const Teuchos::RCP<panzer::IntegrationRule> ir,
        const Teuchos::ParameterList&               userData,
        const Teuchos::RCP<panzer::GlobalData>&      globalData,
        PHX::FieldManager<panzer::Traits>&          fm) const;
      bool createThermodiffCoeffDefault(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults) const;
      bool createFEMGradNegPotential(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults) const;
      bool createSpaceCharge(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults) const;
      bool createGatherScaledFields(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            inputDOFName,
        const Teuchos::ParameterList& userData) const;
      bool createICRemap(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const std::string&            inputDOFName) const;
      bool createICEquilibriumDensity(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const Teuchos::ParameterList& plist) const;
      bool createICFreqDomEquilibriumDensity(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const Teuchos::ParameterList& plist) const;
      bool createICEquilibriumPotential(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key, 
	const Teuchos::ParameterList& userData,
	const Teuchos::RCP<std::vector<std::string>> semBlocks) const;
      bool createICGauss(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const Teuchos::ParameterList& plist) const;
      bool createICFunction(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const std::string&            key,
        const Teuchos::ParameterList& plist) const;
      bool createInitialPotentialGrad(
        EvaluatorVector               evaluators,
        const Teuchos::ParameterList& defaults,
        const Teuchos::ParameterList& userData,
        const Teuchos::ParameterList& cvfem_data) const; 
      void setupMoleFraction(const Teuchos::ParameterList& myModels) const; 

  }; // end of class ClosureModelFactory

} // end of namespace charon

#endif // CHARON_CLOSURE_MODEL_FACTORY_DECL_HPP
