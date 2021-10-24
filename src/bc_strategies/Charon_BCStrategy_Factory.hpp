
#ifndef CHARON_BCSTRATEGY_FACTORY_HPP
#define CHARON_BCSTRATEGY_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"

// Follow the 3 steps below to add a boundary condiiton

// Step #1: Add bcstrategy definition
#include "Charon_BCStrategy_Dirichlet_Constant.hpp"
#include "Charon_BCStrategy_Dirichlet_OhmicContact.hpp"
#include "Charon_BCStrategy_Dirichlet_ContactOnInsulator.hpp"
#include "Charon_BCStrategy_Dirichlet_LinearRamp.hpp"
#include "Charon_BCStrategy_Dirichlet_Sinusoid.hpp"
#include "Charon_BCStrategy_Dirichlet_Periodic.hpp"
#include "Charon_BCStrategy_Dirichlet_BJT1DBaseContact.hpp"
#include "Charon_BCStrategy_Dirichlet_CurrentConstraint.hpp"
#include "Charon_BCStrategy_Dirichlet_ThermalContact.hpp"
#include "Charon_BCStrategy_Dirichlet_MMS.hpp"
#include "Charon_BCStrategy_Neumann_ThermalContact.hpp"
#include "Charon_BCStrategy_Neumann_SurfaceCharge.hpp"
#include "Charon_BCStrategy_Neumann_Constant.hpp"
#include "Charon_BCStrategy_Interface_Simple.hpp"
#include "Charon_BCStrategy_Interface_NeumannMatch.hpp"
#include "Charon_BCStrategy_Interface_Heterojunction.hpp"
#include "Charon_BCStrategy_FreqDom.hpp"
#include "Charon_BCStrategy_Dirichlet_SchottkyContact.hpp"
#include "Charon_BCStrategy_Neumann_SchottkyContact.hpp"

namespace charon {

  // Step #2: Add template builder
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_Constant, BCStrategy_Dirichlet_Constant)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_OhmicContact,
    BCStrategy_Dirichlet_OhmicContact)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_ContactOnInsulator,
    BCStrategy_Dirichlet_ContactOnInsulator)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_LinearRamp, BCStrategy_Dirichlet_LinearRamp)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_Sinusoid, BCStrategy_Dirichlet_Sinusoid)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_MMS, BCStrategy_Dirichlet_MMS)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_Periodic, BCStrategy_Dirichlet_Periodic)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_BJT1DBaseContact,
    BCStrategy_Dirichlet_BJT1DBaseContact)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_CurrentConstraint,
    BCStrategy_Dirichlet_CurrentConstraint)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_ThermalContact,
    BCStrategy_Dirichlet_ThermalContact)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Neumann_ThermalContact,
    BCStrategy_Neumann_ThermalContact)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Neumann_SurfaceCharge,
    BCStrategy_Neumann_SurfaceCharge)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Neumann_Constant, BCStrategy_Neumann_Constant)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Interface_Simple, BCStrategy_Interface_Simple)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Interface_NeumannMatch,
    BCStrategy_Interface_NeumannMatch)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Interface_Heterojunction,
    BCStrategy_Interface_Heterojunction)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_FreqDom,
    BCStrategy_FreqDom)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Dirichlet_SchottkyContact,
    BCStrategy_Dirichlet_SchottkyContact)
  PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(
    charon::BCStrategy_Neumann_SchottkyContact,
    BCStrategy_Neumann_SchottkyContact)

  struct BCFactory : public panzer::BCStrategyFactory {
    bool throwOnFailure_;
    BCFactory(bool throwOnFailure=true)
      : throwOnFailure_(throwOnFailure) {}

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm =
        Teuchos::rcp(new panzer::BCStrategy_TemplateManager<panzer::Traits>);

      bool found = false;

      // Step #3: Create the object
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant",
        BCStrategy_Dirichlet_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Ohmic Contact",
        BCStrategy_Dirichlet_OhmicContact)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Contact On Insulator",
        BCStrategy_Dirichlet_ContactOnInsulator)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Linear Ramp",
        BCStrategy_Dirichlet_LinearRamp)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Sinusoid",
        BCStrategy_Dirichlet_Sinusoid)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("MMS", BCStrategy_Dirichlet_MMS)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Periodic",
        BCStrategy_Dirichlet_Periodic)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("BJT1D Base Contact",
        BCStrategy_Dirichlet_BJT1DBaseContact)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Constant Current",
        BCStrategy_Dirichlet_CurrentConstraint)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Resistor Contact",
        BCStrategy_Dirichlet_CurrentConstraint)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Thermal Contact",
        BCStrategy_Dirichlet_ThermalContact)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Thermal Contact",
        BCStrategy_Neumann_ThermalContact)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Surface Charge",
        BCStrategy_Neumann_SurfaceCharge)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Constant",
        BCStrategy_Neumann_Constant)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Interface Simple",
        BCStrategy_Interface_Simple)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Interface Neumann Match",
        BCStrategy_Interface_NeumannMatch)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Interface Heterojunction",
        BCStrategy_Interface_Heterojunction)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Frequency Domain",
        BCStrategy_FreqDom)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Dirichlet Schottky Contact",
        BCStrategy_Dirichlet_SchottkyContact)
      PANZER_BUILD_BCSTRATEGY_OBJECTS("Neumann Schottky Contact",
        BCStrategy_Neumann_SchottkyContact)

      if(throwOnFailure_) {
        TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error,
                                 "Error - the BC Strategy called \"" << bc.strategy() <<
                                 "\" is not a valid identifier in the BCStrategyFactory.  Either add a valid implementation " <<
                                 "to the factory or fix the input file. " <<
                                 "The relevant boundary condition is:\n\n" << bc << std::endl);
      }

      // return null if the throw on failure has been
      // disabled (running in library mode)
      if(!found)
         return Teuchos::null;

      return bcs_tm;

    }

  };

}

#endif
