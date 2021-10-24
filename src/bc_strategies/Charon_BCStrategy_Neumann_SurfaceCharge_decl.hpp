

#ifndef CHARON_BCSTRATEGY_NEUMANN_SURFACECHARGE_DECL_HPP
#define CHARON_BCSTRATEGY_NEUMANN_SURFACECHARGE_DECL_HPP

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Panzer_BCStrategy_Neumann_DefaultImpl.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PureBasis.hpp"
#include "Phalanx_FieldManager.hpp"


/**
 * @brief The "Neumann Surface Charge" boundary condition has the following settings :
 *
 * <ParameterList>
 *   <Parameter name="Type" type="string" value="Neumann"/>
 *   <Parameter name="Sideset ID" type="string" value="sheet"/>
 *   <Parameter name="Element Block ID" type="string" value="silicon"/>
 *   <Parameter name="Equation Set Name" type="string" value="ELECTRIC_POTENTIAL"/>
 *   <Parameter name="Strategy" type="string" value="Neumann Surface Charge"/>
 *   <ParameterList name="Data">
 *     <Parameter name="Fixed Charge" type="double" value="1e12"/>
 *     <ParameterList name="Polarization">
 *       <Parameter name="Type" type="string" value="piezo" <!-- piezo or both !-->
 *       <Parameter name="Top" type="string" value="AlN"/>
 *       <Parameter name="Bottom" type="string" value="GaN"/>
 *       <Parameter name="Xcomp" type="double" value="0.2"/>
 *       <Parameter name="Scale" type="double" value="1.0"/>
 *     </ParameterList>
 *     <ParameterList name="Surface Trap">  <!-- 50 types of traps can be specified !-->
 *       <Parameter name="Electron Effective Mass" type="double" value="0.32" />  <!-- m0 !-->
 *       <Parameter name="Hole Effective Mass" type="double" value="0.16" />      <!-- m0 !-->
 *       <ParameterList name="Trap 0">
 *         <!-- eV from mid-bandgap, above (+) or below (-) mid-bandgap !-->
 *         <Parameter name="Trap Energy" type="double" value="0.3" />   <!-- eV !-->
 *         <Parameter name="Trap Density" type="double" value="1e11" />  <!-- cm^{-2} !-->
 *         <Parameter name="Trap Type" type="string" value="Acceptor" /> <!-- electron capture, -e when occupied !-->
 *         <Parameter name="Electron Cross Section" type="double" value="1e-12" />  <!-- cm^2 !-->
 *         <Parameter name="Hole Cross Section" type="double" value="1e-12" />      <!-- cm^2 !-->
 *       </ParameterList> 
 *       <ParameterList name="Trap 1">
 *         <!-- eV from mid-bandgap, above (+) or below (-) mid-bandgap !-->
 *         <Parameter name="Trap Energy" type="double" value="-0.2" />  <!-- eV from mid-bandgap !-->
 *         <Parameter name="Trap Density" type="double" value="2e11" />  <!-- cm^{-2} !-->
 *         <Parameter name="Trap Type" type="string" value="Donor" />    <!-- hole capture, +e when occupied, 0 when unoccupied !-->
 *         <Parameter name="Electron Cross Section" type="double" value="1e-12" />  <!-- cm^2 !-->
 *         <Parameter name="Hole Cross Section" type="double" value="1e-14" />      <!-- cm^2 !-->
 *       </ParameterList>          
 *     </ParameterList>
 *     <ParameterList name="Surface Recombination"> 
 *       <Parameter name="Electron Surface Velocity" type="double" value="1e5" />  <!-- cm/s !-->
 *       <Parameter name="Hole Surface Velocity" type="double" value="1e3" />      <!-- cm/s !-->
 *       <Parameter name="Energy Level" type="double" value="0.0" />      <!-- eV from mid-bandgap, above (+) or below (-) mid-bandgap !-->
 *     </ParameterList>
 *   </ParameterList>
 * </ParameterList>
 */
 

namespace charon {

  template <typename EvalT>
  class BCStrategy_Neumann_SurfaceCharge : public panzer::BCStrategy_Neumann_DefaultImpl<EvalT> {

  public:

    BCStrategy_Neumann_SurfaceCharge(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data);

    void setup(const panzer::PhysicsBlock& side_pb,
               const Teuchos::ParameterList& user_data);

    void buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                    const panzer::PhysicsBlock& side_pb,
                                    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
                                    const Teuchos::ParameterList& models,
                                    const Teuchos::ParameterList& user_data) const;

    virtual void
    buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                   const panzer::PhysicsBlock& side_pb,
                                                   const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                   const Teuchos::ParameterList& user_data) const;

    virtual void postRegistrationSetup(typename panzer::Traits::SetupData d,
                                       PHX::FieldManager<panzer::Traits>& vm);

    virtual void evaluateFields(typename panzer::Traits::EvalData d);

  private:
    std::string top, bottom;
    double xcomp;
    double scale;
    
    void initialize(Teuchos::RCP<const Teuchos::ParameterList> plist);
    
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
    
    Teuchos::RCP<Teuchos::ParameterList> surfTrapPList, surfRecombPList,
                                         polarPList; 
    
    double fixedCharge;

    //    std::string varyingCharge;
    
    bool bFixCharge, bVaryingCharge, bSurfTrap, bSurfRecomb, bPolar; 
   
    std::string fluxSurfCharge, fluxSurfRecomb; 

  };

}

#endif
