
#ifndef _CHARON_GatherScaledFields_IMPL_HPP_
#define _CHARON_GatherScaledFields_IMPL_HPP_

#include "Charon_Scaling_Parameters.hpp"

// **********************************************************************
template<typename EvalT, typename Traits>
charon::GatherScaledFields<EvalT, Traits>::
GatherScaledFields(Teuchos::RCP<panzer_stk::STK_Interface const> const& mesh,
                   Teuchos::ParameterList const& p)
{
  mesh_ = mesh;

  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::BasisIRLayout> basis =
    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis");

  Teuchos::RCP<charon::Scaling_Parameters> scale_params =
    p.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");

  scaleFactors_ = scale_params->varScaleFactors;

  gatherFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] =
      PHX::MDField<ScalarT,panzer::Cell,panzer::NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("Gather Scaled STK Fields");
}

// **********************************************************************
template<typename EvalT, typename Traits>
void charon::GatherScaledFields<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData /* d */,
                      PHX::FieldManager<Traits>& /* fm */)
{
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    std::string fieldName = gatherFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(mesh_->getNodeRank(), fieldName);

    if(stkFields_[fd]==0) {
      std::stringstream ss;
      mesh_->printMetaData(ss);
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                 "charon::GatherScaledFields: STK field " << "\"" << fieldName << "\" "
                                 "not found.\n STK meta data follows: \n\n" << ss.str());
    }
  }
}

// **********************************************************************
template<typename EvalT, typename Traits>
void charon::GatherScaledFields<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
   const std::vector<stk::mesh::Entity> & localElements = *mesh_->getElementsOrderedByLID();

   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;

   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
      stk::mesh::Entity const* relations = mesh_->getBulkData()->begin_nodes(localElements[cellLocalId]);

      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         VariableField * field = stkFields_[fieldIndex];

         std::size_t basisCnt = gatherFields_[fieldIndex].dimension(1);

         // If the variables were unscaled then they need to be scaled
         // prior to being used as the initial guess. This is controlled
         // by a flag in the input file for the initial conditions. Get
         // the scale factor for the variable and scale it
         double scale_factor = 1.0;
         std::map<std::string,double>::const_iterator sf_ent =
           (*scaleFactors_).find(gatherFields_[fieldIndex].fieldTag().name());
         if(sf_ent != (*scaleFactors_).end())
         {
           scale_factor = 1.0/(sf_ent->second);
         }

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<basisCnt;basis++) {
            stk::mesh::Entity node = relations[basis];

            // Scale the value and gather it to the local data structure.
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = *stk::mesh::field_data(*field,node) * scale_factor; // from STK
         }
      }
   }
}

#endif // _CHARON_GatherScaledFields_IMPL_HPP_
