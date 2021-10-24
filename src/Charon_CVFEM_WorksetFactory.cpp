
#include "Charon_CVFEM_WorksetFactory.hpp"

#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_SetupUtilities.hpp"

namespace charon {

void CVFEM_WorksetFactory::setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh)
{
   mesh_ = mesh;
}

/** Build sets of boundary condition worksets
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> > CVFEM_WorksetFactory::
getSideWorksets(const panzer::BC & bc,
              const panzer::PhysicsBlock & pb) const
{
  panzer::WorksetNeeds needs = getNeedsFromPhysicsBlock(pb);
  return getSideWorksets(bc,needs);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > CVFEM_WorksetFactory::
getSideWorksets(const panzer::BC & bc,
                const panzer::WorksetNeeds & needs) const
{
  panzer::WorksetDescriptor desc(bc.elementBlockID(),bc.sidesetID());
  return getSideWorksets(desc,needs);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > CVFEM_WorksetFactory::
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

Teuchos::RCP<std::vector<panzer::Workset> > CVFEM_WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::PhysicsBlock & pb) const
{
  panzer::WorksetNeeds needs = getNeedsFromPhysicsBlock(pb);

  return getWorksets(worksetDesc,needs);
}

//////////////////////////////////////////////////////////////////////////////////////

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
CVFEM_WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
                const panzer::WorksetNeeds & needs_a,
                const panzer::WorksetNeeds & needs_b) const
{
  Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets;
  worksets = panzer_stk::buildBCWorksets(*mesh_, needs_a, desc.getElementBlock(0),
                                                 needs_b, desc.getElementBlock(1),
                                                 desc.getSideset(0));
  // Must have this check, otherwise, worksets may be equal to 0, which causes 
  // seg fault, on some processors when running in parallel.
  if (worksets != Teuchos::null)
  {
    for (std::map<unsigned,panzer::Workset>::iterator wkst = (*worksets).begin(); wkst != (*worksets).end(); ++wkst) 
    {
      panzer::Workset& works = wkst->second;

      // add points to primary side
      addCVPointsAndBasisBoundary(works.num_cells,needs_a, works.details(0));
 
      // add points to secondary side
      addCVPointsAndBasisBoundary(works.num_cells, needs_b, works.details(1),Teuchos::rcpFromRef(works.details(0)));
    }
  }

  return worksets;
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
CVFEM_WorksetFactory::
getSideWorksets(const panzer::WorksetDescriptor & desc,
  const panzer::WorksetNeeds & needs) const
{
  Teuchos::RCP<std::map<unsigned,panzer::Workset> > worksets;
  worksets = panzer_stk::buildBCWorksets(*mesh_,needs,desc.getElementBlock(),desc.getSideset());

  // Must have this check, otherwise, worksets may be equal to 0, which causes 
  // seg fault, on some processors when running in parallel.
  if (worksets != Teuchos::null)   
  {
    for (std::map<unsigned,panzer::Workset>::iterator wkst = (*worksets).begin(); wkst != (*worksets).end(); ++wkst) 
    {
      panzer::Workset& works = wkst->second;
      
      // add boundary points 
      addCVPointsAndBasisBoundary(works.num_cells,needs,works.details(0));
    }
  }

  return worksets; 
}

Teuchos::RCP<std::vector<panzer::Workset> > CVFEM_WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::WorksetNeeds & needs) const
{
  Teuchos::RCP< std::vector<panzer::Workset> > worksets;

  if(!worksetDesc.useSideset()) {
    worksets = panzer_stk::buildWorksets(*mesh_,worksetDesc.getElementBlock(), needs);

     for (std::vector<panzer::Workset>::iterator wkst = (*worksets).begin(); wkst != (*worksets).end(); ++wkst) {
         addCVPointsAndBasis(wkst->num_cells,needs,*wkst);
     }

  }
  else if(worksetDesc.useSideset() && worksetDesc.sideAssembly()) {
    // uses cascade by default, each subcell has its own workset
     worksets = panzer_stk::buildWorksets(*mesh_,needs,worksetDesc.getSideset(),worksetDesc.getElementBlock(),true);

     for (std::vector<panzer::Workset>::iterator wkst = (*worksets).begin(); wkst != (*worksets).end(); ++wkst) {
         addCVPointsAndBasis(wkst->num_cells,needs,*wkst);
     }

  }
  else {
    TEUCHOS_ASSERT(false);
  }
    return worksets;
}

void CVFEM_WorksetFactory::addCVPointsAndBasis(std::size_t num_cells,
                                               const panzer::WorksetNeeds &needs,
                                               panzer::WorksetDetails &details) const
{

    // vector of integration degrees and basis names
    Teuchos::RCP<std::vector<int> > ir_degrees = Teuchos::rcp(new std::vector<int>(0));
    Teuchos::RCP<std::vector<std::string> > basis_names = Teuchos::rcp(new std::vector<std::string>(0));

    // integration rules currently in workset
    std::vector<Teuchos::RCP<const panzer::IntegrationRule> > int_rules = needs.int_rules;

    // pure bases currently in workset
    std::vector<Teuchos::RCP<const panzer::PureBasis> > bases = needs.bases;

    panzer::CellData cell_data(num_cells,needs.cellData.getCellTopology());

    // create pure HCurl basis - used for stabilization
    std::size_t cells_size = needs.cellData.numCells();
    Teuchos::RCP<panzer::PureBasis> Hcurlbasis = Teuchos::rcp(new panzer::PureBasis("HCurl",1,cells_size,needs.cellData.getCellTopology()));

    // add to bases
    bases.push_back(Hcurlbasis);

    // create volume integration rule
    std::string cvfem_type = "volume";
    Teuchos::RCP<const panzer::IntegrationRule> cvfem_vol_rule = Teuchos::rcp(new panzer::IntegrationRule(needs.cellData,cvfem_type));
    details.ir_degrees->push_back(cvfem_vol_rule->order());

    // create side integration rule
    cvfem_type = "side";
    Teuchos::RCP<const panzer::IntegrationRule> cvfem_side_rule = Teuchos::rcp(new panzer::IntegrationRule(needs.cellData,cvfem_type));
    details.ir_degrees->push_back(cvfem_side_rule->order());

    // Check for an empty workset. This happens, for example, when a
    // processor doesn't have any cells from a physics block.
    bool empty_wkst = (details.numOwnedCells() <= 0);

    // loop over bases to set vector of basis_names
    for(std::size_t b=0;b<bases.size();b++) {

        // get combined basis integration rule layouts
        Teuchos::RCP<panzer::BasisIRLayout> b_vol_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_vol_rule));
        Teuchos::RCP<panzer::BasisIRLayout> b_side_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_side_rule));

        details.basis_names->push_back(b_vol_layout->name());
        details.basis_names->push_back(b_side_layout->name());

    }

    // get integration values
     Teuchos::RCP<panzer::IntegrationValues2<double> > cvfem_iv_vol =
               Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));;
     cvfem_iv_vol->setupArrays(cvfem_vol_rule);
     cvfem_iv_vol->evaluateValues(details.cell_vertex_coordinates);

     Teuchos::RCP<panzer::IntegrationValues2<double> > cvfem_iv_side =
               Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));;
     cvfem_iv_side->setupArrays(cvfem_side_rule);
     cvfem_iv_side->evaluateValues(details.cell_vertex_coordinates);

     details.int_rules.push_back(cvfem_iv_vol);
     details.int_rules.push_back(cvfem_iv_side);

     // loop over pure bases
     for(std::size_t b=0;b<bases.size();b++) {

        // get combined basis integration rule layouts
        Teuchos::RCP<panzer::BasisIRLayout> b_vol_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_vol_rule));
        Teuchos::RCP<panzer::BasisIRLayout> b_side_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_side_rule));

        // add basis values for cvfem volume points
        std::size_t vol_index =
             std::distance(details.ir_degrees->begin(),
                       std::find(details.ir_degrees->begin(),
                                 details.ir_degrees->end(),
                                 cvfem_vol_rule->order()));

        Teuchos::RCP<panzer::BasisValues2<double> > bv2_vol =
                   Teuchos::rcp(new panzer::BasisValues2<double>("",true,false));

        bv2_vol->setupArrays(b_vol_layout);

        bv2_vol->evaluateValuesCV(details.int_rules[vol_index]->ref_ip_coordinates,
                                  details.int_rules[vol_index]->jac,
                                  details.int_rules[vol_index]->jac_det,
                                  details.int_rules[vol_index]->jac_inv,
                                  details.cell_vertex_coordinates,
                                  !empty_wkst,
                                  num_cells);

        details.bases.push_back(bv2_vol);

       // add basis values for cvfem side points
         std::size_t side_index =
             std::distance(details.ir_degrees->begin(),
                      std::find(details.ir_degrees->begin(),
                                details.ir_degrees->end(),
                                cvfem_side_rule->order()));

         Teuchos::RCP<panzer::BasisValues2<double> > bv2_side =
                     Teuchos::rcp(new panzer::BasisValues2<double>("",true,false));

         bv2_side->setupArrays(b_side_layout);
         bv2_side->evaluateValuesCV(details.int_rules[side_index]->ref_ip_coordinates,
                                    details.int_rules[side_index]->jac,
                                    details.int_rules[side_index]->jac_det,
                                    details.int_rules[side_index]->jac_inv,
                                    details.cell_vertex_coordinates,
                                    !empty_wkst,
                                    num_cells);

         details.bases.push_back(bv2_side);

      } // loop over bases

} //addCVPointsAndBasis

void CVFEM_WorksetFactory::addCVPointsAndBasisBoundary(std::size_t num_cells,
                                                       const panzer::WorksetNeeds &needs,
                                                       panzer::WorksetDetails & details,
                                                       const Teuchos::RCP<panzer::WorksetDetails> other_details) const
{

    const panzer::CellData side_cell_data(num_cells,
                                          details.subcell_index,
                                          needs.cellData.getCellTopology());


    std::vector<Teuchos::RCP<const panzer::PureBasis> > bases;
    for(std::size_t i=0;i<needs.bases.size();i++)
      bases.push_back(Teuchos::rcp(new panzer::PureBasis(needs.bases[i]->type(),needs.bases[i]->order(),side_cell_data)));

    // create boundary integration rule
    std::string cvfem_type = "boundary";
    Teuchos::RCP<const panzer::IntegrationRule> cvfem_bc_rule = Teuchos::rcp(new panzer::IntegrationRule(side_cell_data,cvfem_type));
    details.ir_degrees->push_back(cvfem_bc_rule->order());

    // Check for an empty workset. This happens, for example, when a
    // processor doesn't have any cells from a physics block.
    bool empty_wkst = (details.numOwnedCells() <= 0);

    // loop over bases to set vector of basis_names
    for(std::size_t b=0;b<bases.size();b++) {

        // get combined basis integration rule layouts
        Teuchos::RCP<panzer::BasisIRLayout> b_bc_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_bc_rule));
        details.basis_names->push_back(b_bc_layout->name());

    }

    // get integration values
     Teuchos::RCP<panzer::IntegrationValues2<double> > cvfem_iv_bc =
               Teuchos::rcp(new panzer::IntegrationValues2<double>("",true));;
     cvfem_iv_bc->setupArrays(cvfem_bc_rule);

    // get corresponding integration rule index on other side
     if (Teuchos::nonnull(other_details))
     {
         std::size_t int_index =
             std::distance(other_details->ir_degrees->begin(),
                  std::find(other_details->ir_degrees->begin(),
                         other_details->ir_degrees->end(),
                         cvfem_bc_rule->order()));

        cvfem_iv_bc->evaluateValues(details.cell_vertex_coordinates, other_details->int_rules[int_index]->ip_coordinates);
     }
     else
        cvfem_iv_bc->evaluateValues(details.cell_vertex_coordinates);

     details.int_rules.push_back(cvfem_iv_bc);

     // loop over pure bases
     for(std::size_t b=0;b<bases.size();b++) {

        // get combined basis integration rule layouts
        Teuchos::RCP<panzer::BasisIRLayout> b_bc_layout = Teuchos::rcp(new panzer::BasisIRLayout(bases[b],*cvfem_bc_rule));

        // add basis values for cvfem boundary points
        std::size_t bc_index =
             std::distance(details.ir_degrees->begin(),
                       std::find(details.ir_degrees->begin(),
                                 details.ir_degrees->end(),
                                 cvfem_bc_rule->order()));

        Teuchos::RCP<panzer::BasisValues2<double> > bv2_bc =
                   Teuchos::rcp(new panzer::BasisValues2<double>("",true,false));

        bv2_bc->setupArrays(b_bc_layout);

        bv2_bc->evaluateValuesCV(details.int_rules[bc_index]->ref_ip_coordinates,
                                 details.int_rules[bc_index]->jac,
                                 details.int_rules[bc_index]->jac_det,
                                 // details.int_rules[bc_index]->jac_inv);  // cannot get coordinates if used
                                 details.int_rules[bc_index]->jac_inv,
                                 details.cell_vertex_coordinates,
                                 !empty_wkst,
                                 num_cells);

        details.bases.push_back(bv2_bc);

      } // loop over bases

} //addCVPointsAndBasisBoundary

panzer::WorksetNeeds CVFEM_WorksetFactory::
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


} // namespace Example
