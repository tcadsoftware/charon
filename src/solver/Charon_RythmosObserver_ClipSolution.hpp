
#ifndef CHARON_RYTHMOS_OBSERVER_CLIP_SOLUTION_HPP
#define CHARON_RYTHMOS_OBSERVER_CLIP_SOLUTION_HPP

#include "Charon_config.hpp"

#include "Rythmos_StepperBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_TimeRange.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_STK_Interface.hpp"

#include <vector>

namespace charon {

class RythmosObserver_ClipSolution : public Rythmos::IntegrationObserverBase<double> {
public:

   RythmosObserver_ClipSolution(const Teuchos::RCP<const panzer::GlobalIndexer>& dof_manager,
                                const Teuchos::RCP<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >& lof,
                                const std::vector<std::string> & clip_dof_names) :
      m_dof_manager(dof_manager),
      m_lof(lof),
      m_clip_dof_names(clip_dof_names)
   { }

   Teuchos::RCP<Rythmos::IntegrationObserverBase<double> >
   cloneIntegrationObserver() const
   { return Teuchos::rcp(new RythmosObserver_ClipSolution(m_dof_manager, m_lof, m_clip_dof_names)); }

   void resetIntegrationObserver(const Rythmos::TimeRange<double> & /* integrationTimeDomain */) { }

   void observeStartTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */) { }

   void observeCompletedTimeStep(const Rythmos::StepperBase<double> &stepper,
                                 const Rythmos::StepControlInfo<double> & /* stepCtrlInfo */,
                                 const int /* timeStepIter */)
   {
      // build field number, but only do it once!
      if(m_field_ids.size()!=m_clip_dof_names.size()) {
         m_field_ids.clear();
         for(std::size_t i=0;i<m_clip_dof_names.size();i++)
            m_field_ids.push_back(m_dof_manager->getFieldNum(m_clip_dof_names[i]));
      }

      // clip current solution
      clipSolution(stepper);
   }

   void observeEndTimeIntegration(const Rythmos::StepperBase<double> & /* stepper */) { }

private:

   void clipSolution(const Rythmos::StepperBase<double> &stepper)
   {
      std::vector<panzer::GlobalOrdinal> GIDs;
      std::vector<int> LIDs;

      // Grab solution vector
      Teuchos::RCP<Thyra::VectorBase<double> > solution = Teuchos::rcp_const_cast<Thyra::VectorBase<double> >(stepper.getStepStatus().solution);
      const Epetra_Map & ep_map = *m_lof->getMap(0);
      Teuchos::RCP<Epetra_Vector> ep_solution = Thyra::get_Epetra_Vector(ep_map, solution);

      // grab element blocks
      std::vector<std::string> eBlocks;
      m_dof_manager->getElementBlockIds(eBlocks);
      TEUCHOS_ASSERT(eBlocks.size()>0);

      // loop over element blocks clipping as necessary
      for(std::size_t blk=0;blk<eBlocks.size();blk++) {
         std::string blockId = eBlocks[blk];
         const std::vector<int> & localCellIds = m_dof_manager->getElementBlock(blockId);

         // NOTE: A reordering of these loops will likely improve performance
         //       The "getGIDFieldOffsets may be expensive.  However the
         //       "getElementGIDs" can be cheaper. However the lookup for LIDs
         //       may be more expensive!

         // gather operation for each cell in workset
         for(std::size_t cellIndex=0;cellIndex<localCellIds.size();++cellIndex) {
            std::size_t cellLocalId = localCellIds[cellIndex];

            m_dof_manager->getElementGIDs(cellLocalId,GIDs,blockId);

            // caculate the local IDs for this element
            LIDs.resize(GIDs.size());
            for(std::size_t i=0;i<GIDs.size();i++)
               LIDs[i] = ep_map.LID(GIDs[i]);

            // loop over the fields to be gathered
            for (std::size_t fieldIndex=0; fieldIndex<m_field_ids.size();fieldIndex++) {
               int fieldNum = m_field_ids[fieldIndex];
               const std::vector<int> & elmtOffset = m_dof_manager->getGIDFieldOffsets(blockId,fieldNum);

               // loop over basis functions and fill the fields
               for(std::size_t basis=0;basis<elmtOffset.size();basis++) {
                  int offset = elmtOffset[basis];
                  int lid = LIDs[offset];

                  // is the GID even on this processor (if not -- lid=-1, move on!)
                  if(lid<0) continue;

                  double value = (*ep_solution)[lid];
                  if(value<0.0)
                     (*ep_solution)[lid] = 0.0;
               }
            }
         }
      }
   }


   Teuchos::RCP<const panzer::GlobalIndexer> m_dof_manager;
   Teuchos::RCP<const panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> > m_lof;

   std::vector<std::string> m_clip_dof_names;
   std::vector<int> m_field_ids;
};

}

#endif
