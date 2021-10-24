#ifndef CHARON_CLOSUREMODEL_FACTORY_COMPOSITE_IMPL_HPP
#define CHARON_CLOSUREMODEL_FACTORY_COMPOSITE_IMPL_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include "Panzer_GlobalData.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Charon_Names.hpp"
#include "Charon_FreqDom_Parameters.hpp"

// ********************************************************************
template<typename EvalT>
charon::ClosureModelFactoryComposite<EvalT>::
ClosureModelFactoryComposite(const std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >& factories) :
  m_factories(factories)
{ }

// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
charon::ClosureModelFactoryComposite<EvalT>::
buildClosureModels(const std::string& model_id,
		   const Teuchos::ParameterList& models, 
		   const panzer::FieldLayoutLibrary& fl,
		   const Teuchos::RCP<panzer::IntegrationRule>& ir,
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{

  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

  if (!models.isSublist(model_id)) {
    std::stringstream msg;
    msg << "Falied to find requested model, \"" << model_id 
	<< "\" for equation set:\n" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
  }

  const ParameterList& my_model = models.sublist(model_id);

  // pull out the nonlist (not associated with any model) parameters
  // this will be used by each stored closure model
  Teuchos::ParameterList nonlist_params(models.name()); // make sure it maintains the models name
  for (ParameterList::ConstIterator model_it = models.begin(); 
       model_it != models.end(); ++model_it) {

    std::string key = model_it->first;
    if(!model_it->second.isList())
      nonlist_params.setEntry(key,model_it->second);
  }

  // build a copy of parameter list containing only the closure model of current relevance
  // with any supplemental non-list information contained in the parameter list
  ParameterList copy_of_my_model = nonlist_params;
  copy_of_my_model.sublist(model_id) = my_model; // copy my_model into copy of models

  // Get fd_suffixes if performing a frequency domain simulation
  std::vector<std::string> fd_suffixes = {};
  bool isFreqDom = default_params.get<std::string>("Type") == "Frequency Domain";

  //std::cout << "CompositeClosureModelFactory is in frequency domain mode: " << ( isFreqDom ? "true" : "false" ) << std::endl;
  const Teuchos::RCP<FreqDomParameters> freqDomParamsRCP = (isFreqDom ? 
							 default_params.sublist("Options").get<Teuchos::RCP<FreqDomParameters> >("Frequency Domain Parameters") 
                                                           : 
                                                         Teuchos::rcp(new FreqDomParameters()) ); 

  // TODO: take and use the suffixes from freqDomParams
  if(isFreqDom){
    // TEMP: hard coded indices                                                                                                                                                                                                                              
    for(int i = 0 ; i < freqDomParamsRCP->getNumTimeCollocationPoints() ; i++)
      fd_suffixes.push_back("_TP" + std::to_string(i) + "_");
    // TEMP: create CMF fields in frequency domain, too
    //for(int i = 0 ; i < freqDomParamsRCP->getNumTotalHarmonics() ; i++)
    //  fd_suffixes.push_back("_CosH" + std::to_string(i) + "_");
    //for(int i = 0 ; i < freqDomParamsRCP->getNumTotalHarmonics() ; i++)
    //  fd_suffixes.push_back("_SinH" + std::to_string(i) + "_");
  }
  else // !isFreqDom
    fd_suffixes = {""};

  // Loop over factories
  int factory_index = 0;
  for (std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >::const_iterator factory = m_factories.begin(); factory != m_factories.end(); ++factory) {

    //std::cout << "The current factory index is: " << std::to_string(factory_index) << std::endl;

    // HB mod: are these discontinuous fields and discontinuous suffixes necessary? 
    std::string discfields = ( default_params.isParameter("Discontinuous Fields") ? default_params.get<std::string>("Discontinuous Fields") : "");
    std::string discsuffix = ( default_params.isParameter("Discontinuous Suffix") ? default_params.get<std::string>("Discontinuous Suffix") : "");

    // HB mod: here, we set the suffixes for the frequency domain names, to initialize a charon::names with that suffix 
    ParameterList modifiedPL = default_params;
    modifiedPL.set<Teuchos::RCP<const charon::Names >>("Names", Teuchos::rcp(new charon::Names(1 /* equation_dimension */, 
                                                                                               "" /* prefix */, discfields, discsuffix, 
                                                                                               fd_suffixes[factory_index])));

    (*factory)->getAsObject<EvalT>()->setThrowOnModelNotFound(false);
    RCP< vector< RCP<Evaluator<panzer::Traits> > > > tmp_evaluators =
      (*factory)->getAsObject<EvalT>()->buildClosureModels(model_id, models /* copy_of_my_model was used in panzer::ClosureModel_Factory_Composite */, fl, ir,
                                                           (isFreqDom ? modifiedPL : default_params) /* use the unmodified PL if !isFreqDom */,
                                                           user_data,global_data,fm);

    if (tmp_evaluators->size() > 0) {
      for (vector< RCP<Evaluator<panzer::Traits> > >::const_iterator eval = tmp_evaluators->begin(); eval != tmp_evaluators->end(); ++eval)
        evaluators->push_back(*eval);

    }

    factory_index++;
  }

/*
  // for each model, try each factory until you get a nonnull
  // return, meaning you have built the evaluators for that model
  for (ParameterList::ConstIterator model_it = my_model.begin(); 
       model_it != my_model.end(); ++model_it) {
    
    std::string model_key = model_it->first;
    
    // Duplicate models sublist with just the particular model you
    // want to build
    ParameterList copy_of_models = nonlist_params;
    // Teuchos::ParameterList* tmp;
    // copy_of_models.sublist(model_id).sublist(model_key) = model_it->second.getValue(tmp);
    copy_of_models.sublist(model_id).setEntry(model_key,model_it->second);

    std::cout << "COPY OF MODELS = " << model_id << std::endl;
    copy_of_models.print(std::cout);
    
    // Loop over factories
    for (std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >::const_iterator factory = m_factories.begin(); factory != m_factories.end(); ++factory) {

      RCP< vector< RCP<Evaluator<panzer::Traits> > > > tmp_evaluators =
	(*factory)->getAsObject<EvalT>()->buildClosureModels(model_id,set,copy_of_models,default_params,user_data,global_data,fm);

      if (tmp_evaluators->size() > 0) {
	
	for (vector< RCP<Evaluator<panzer::Traits> > >::const_iterator eval = tmp_evaluators->begin(); eval != tmp_evaluators->end(); ++eval)
	  evaluators->push_back(*eval);

      }
      
    }

  }
  */

  return evaluators;
}

#endif
