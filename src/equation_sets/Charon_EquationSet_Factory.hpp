
#ifndef __Charon_EquationSet_Factory_hpp__
#define __Charon_EquationSet_Factory_hpp__

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_GlobalData.hpp"

// Follow Steps 1 through 3 below to add a new equation set

// Step #1: Add equation set declaration here
#include "Charon_EquationSet_Laplace.hpp"
#include "Charon_EquationSet_NLPoisson.hpp"
#include "Charon_EquationSet_DriftDiffusion.hpp"
#include "Charon_EquationSet_EFFPG_DriftDiffusion.hpp"
#include "Charon_EquationSet_SGCVFEM_DriftDiffusion.hpp"
#include "Charon_EquationSet_SGCVFEM_Laplace.hpp"
#include "Charon_EquationSet_SGCVFEM_NLPoisson.hpp"
#include "Charon_EquationSet_SGCharon1_DriftDiffusion.hpp"
#include "Charon_EquationSet_DDLattice.hpp"
#include "Charon_EquationSet_Lattice.hpp"
#include "Charon_EquationSet_DDIon.hpp"
#include "Charon_EquationSet_DDIonLattice.hpp"
#include "Charon_EquationSet_EFFPG_DDIonLattice.hpp"
#include "Charon_EquationSet_FreqDom.hpp"

namespace charon {

  // Step #2: Add TempalteBuilder declaration here
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_Laplace,
    EquationSet_Laplace)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_SGCVFEM_Laplace,
    EquationSet_SGCVFEM_Laplace)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_NLPoisson,
    EquationSet_NLPoisson)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_SGCVFEM_NLPoisson,
    EquationSet_SGCVFEM_NLPoisson)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_DriftDiffusion,
    EquationSet_DriftDiffusion)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(
    charon::EquationSet_EFFPG_DriftDiffusion, EquationSet_EFFPG_DriftDiffusion)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(
    charon::EquationSet_SGCVFEM_DriftDiffusion,
    EquationSet_SGCVFEM_DriftDiffusion)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(
    charon::EquationSet_SGCharon1_DriftDiffusion,
    EquationSet_SGCharon1_DriftDiffusion)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_DDLattice,
    EquationSet_DDLattice)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_Lattice,
    EquationSet_Lattice)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_DDIon,
    EquationSet_DDIon)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_DDIonLattice,
    EquationSet_DDIonLattice)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_EFFPG_DDIonLattice,
    EquationSet_EFFPG_DDIonLattice)
  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(charon::EquationSet_FreqDom,
    EquationSet_FreqDom)

  class EquationSet_Factory : public panzer::EquationSetFactory
  {
    bool throwOnFailure_;
  public:

    EquationSet_Factory(bool throwOnFailure=true) : throwOnFailure_(throwOnFailure) {}

    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
                     const int& default_integration_order,
                     const panzer::CellData& cell_data,
                     const Teuchos::RCP<panzer::GlobalData>& global_data,
                     const bool build_transient_support) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set=
        Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);

      bool found = false;

      // Step #3: Add object builder here
      PANZER_BUILD_EQSET_OBJECTS("Laplace", EquationSet_Laplace)
      PANZER_BUILD_EQSET_OBJECTS("SGCVFEM Laplace", EquationSet_SGCVFEM_Laplace)
      PANZER_BUILD_EQSET_OBJECTS("NLPoisson", EquationSet_NLPoisson)
      PANZER_BUILD_EQSET_OBJECTS("SGCVFEM NLPoisson",
        EquationSet_SGCVFEM_NLPoisson)
      PANZER_BUILD_EQSET_OBJECTS("Drift Diffusion", EquationSet_DriftDiffusion)
      PANZER_BUILD_EQSET_OBJECTS("EFFPG Drift Diffusion",
        EquationSet_EFFPG_DriftDiffusion)
      PANZER_BUILD_EQSET_OBJECTS("SGCVFEM Drift Diffusion",
        EquationSet_SGCVFEM_DriftDiffusion)
      PANZER_BUILD_EQSET_OBJECTS("SGCharon1 Drift Diffusion",
        EquationSet_SGCharon1_DriftDiffusion)
      PANZER_BUILD_EQSET_OBJECTS("DDLattice", EquationSet_DDLattice)
      PANZER_BUILD_EQSET_OBJECTS("Lattice", EquationSet_Lattice)
      PANZER_BUILD_EQSET_OBJECTS("DDIon", EquationSet_DDIon)
      PANZER_BUILD_EQSET_OBJECTS("DDIonLattice", EquationSet_DDIonLattice)
      PANZER_BUILD_EQSET_OBJECTS("EFFPG DDIonLattice",
        EquationSet_EFFPG_DDIonLattice)
      PANZER_BUILD_EQSET_OBJECTS("Frequency Domain", EquationSet_FreqDom)

      if (!found && throwOnFailure_)
      {
        std::string msg = "Error - the \"Equation Set\" called \"" +
          params->get<std::string>("Type") +
          "\" is not a valid equation set identifier. Please supply the correct factory.\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }

      if(!found)
        return Teuchos::null;

      return eq_set;
    }

  };

}

#endif
