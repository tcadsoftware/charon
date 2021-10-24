
#ifndef CHARON_EFFPG_WORKSETFACTORY_HPP
#define CHARON_EFFPG_WORKSETFACTORY_HPP

#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetNeeds.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"

#include "Panzer_STK_Interface.hpp"

namespace charon {

class EFFPG_WorksetFactory : public panzer_stk::WorksetFactory {

public:

   EFFPG_WorksetFactory() {}

   EFFPG_WorksetFactory(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh) : mesh_(mesh) {}

   virtual ~EFFPG_WorksetFactory() {}

   /** Set mesh
     */
   virtual
   void setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh);

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::BC & bc,
                   const panzer::PhysicsBlock & pb) const;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::BC & bc,
                   const panzer::WorksetNeeds & needs) const;

   /** Build sets of boundary condition worksets for the BCT_Interface case.
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::BC & bc,
                   const panzer::PhysicsBlock & pb_a,
                   const panzer::PhysicsBlock & pb_b) const;

   /** Build worksets specified by the workset descriptor.
     */
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getWorksets(const panzer::WorksetDescriptor & worksetDesc,
               const panzer::PhysicsBlock & pb) const;

   /////////////////////////////////////////////////////////////////////////////////////

   /** Build sets of boundary condition worksets for an interface case.
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::WorksetDescriptor & desc,
                   const panzer::WorksetNeeds & needs_a,
                   const panzer::WorksetNeeds & needs_b) const;

   /** Build sets of boundary condition worksets
     */
   virtual
   Teuchos::RCP<std::map<unsigned,panzer::Workset> >
   getSideWorksets(const panzer::WorksetDescriptor & desc,
                   const panzer::WorksetNeeds & needs) const;

   /** Build worksets specified by the workset descriptor.
     */
   virtual
   Teuchos::RCP<std::vector<panzer::Workset> >
   getWorksets(const panzer::WorksetDescriptor & worksetDesc,
               const panzer::WorksetNeeds & needs) const;

private:

   //! Temporary method extracting needs from the physics block
   panzer::WorksetNeeds getNeedsFromPhysicsBlock(const panzer::PhysicsBlock & pb) const;

   Teuchos::RCP<const panzer_stk::STK_Interface> mesh_;

};

}

#endif
