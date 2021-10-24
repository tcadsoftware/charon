
#include "Charon_EFFPG_WorksetFactory.hpp"

#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

namespace charon {

void EFFPG_WorksetFactory::setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh)
{
   mesh_ = mesh;
}

/** Build sets of boundary condition worksets
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> > EFFPG_WorksetFactory::
getSideWorksets(const panzer::BC & bc,
                const panzer::PhysicsBlock & pb) const
{
  panzer::WorksetNeeds needs = getNeedsFromPhysicsBlock(pb);

  return getSideWorksets(bc,needs);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > EFFPG_WorksetFactory::
getSideWorksets(const panzer::BC & bc,
                const panzer::WorksetNeeds & needs) const
{
  panzer::WorksetDescriptor desc(bc.elementBlockID(),bc.sidesetID());
  return getSideWorksets(desc,needs);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > EFFPG_WorksetFactory::
getSideWorksets(const panzer::BC & bc,
                const panzer::PhysicsBlock & pb_a,
                const panzer::PhysicsBlock & pb_b) const
{
  TEUCHOS_ASSERT(bc.bcType() == panzer::BCT_Interface);

  panzer::WorksetNeeds needs_a = getNeedsFromPhysicsBlock(pb_a);
  panzer::WorksetNeeds needs_b = getNeedsFromPhysicsBlock(pb_b);

  panzer::WorksetDescriptor desc(pb_a.elementBlockID(),pb_b.elementBlockID(),bc.sidesetID(),bc.sidesetID());
  return getSideWorksets(desc,needs_a,needs_b);
}

Teuchos::RCP<std::vector<panzer::Workset> > EFFPG_WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::PhysicsBlock & pb) const
{
  panzer::WorksetNeeds needs = getNeedsFromPhysicsBlock(pb);

  return getWorksets(worksetDesc,needs);
}

/////////////////////////////////////////////////////////////////////////////////////

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
EFFPG_WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
                const panzer::WorksetNeeds & needs_a,
                const panzer::WorksetNeeds & needs_b) const
{
  return panzer_stk::buildBCWorksets(*mesh_,
                                     needs_a, desc.getElementBlock(0),
                                     needs_b, desc.getElementBlock(1),
                                     desc.getSideset(0));
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
EFFPG_WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
  const panzer::WorksetNeeds & needs) const
{
  panzer::WorksetNeeds needs2;
  needs2.cellData = needs.cellData;
  needs2.int_rules = needs.int_rules;
  needs2.bases = needs.bases;

  // add pure hcurl basis to needs
  std::size_t cells_size = needs.cellData.numCells();
  Teuchos::RCP<const panzer::PureBasis> Hcurlbasis = Teuchos::rcp(new panzer::PureBasis("HCurl",1,cells_size,needs.cellData.getCellTopology()));
  needs2.bases.push_back(Hcurlbasis);

  return panzer_stk::buildBCWorksets(*mesh_,needs2,desc.getElementBlock(),desc.getSideset());
}

Teuchos::RCP<std::vector<panzer::Workset> > EFFPG_WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::WorksetNeeds & needs) const
{
  Teuchos::RCP< std::vector<panzer::Workset> > worksets;

  panzer::WorksetNeeds needs2;
  needs2.cellData = needs.cellData;
  needs2.int_rules = needs.int_rules;
  needs2.bases = needs.bases;

  // add pure hcurl basis to needs
   std::size_t cells_size = needs.cellData.numCells();
   Teuchos::RCP<const panzer::PureBasis> Hcurlbasis = Teuchos::rcp(new panzer::PureBasis("HCurl",1,cells_size,needs.cellData.getCellTopology()));
   needs2.bases.push_back(Hcurlbasis);

  if(!worksetDesc.useSideset()) {
     return panzer_stk::buildWorksets(*mesh_,worksetDesc.getElementBlock(), needs2);
  }
  else if(worksetDesc.useSideset() && worksetDesc.sideAssembly()) {
    // uses cascade by default, each subcell has its own workset
    return panzer_stk::buildWorksets(*mesh_,needs2,worksetDesc.getSideset(),worksetDesc.getElementBlock(),true);
  }
  else {
    TEUCHOS_ASSERT(false);
  }
}

panzer::WorksetNeeds EFFPG_WorksetFactory::
getNeedsFromPhysicsBlock(const panzer::PhysicsBlock & pb) const
{
  using Teuchos::RCP;

  panzer::WorksetNeeds needs;

  needs.cellData = pb.cellData();

  const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();
  for(std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
      ir_itr != int_rules.end(); ++ir_itr)
    needs.int_rules.push_back(ir_itr->second);

  const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();
  for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
      b_itr != bases.end(); ++b_itr)
    needs.bases.push_back(b_itr->second);

  return needs;
}

}
