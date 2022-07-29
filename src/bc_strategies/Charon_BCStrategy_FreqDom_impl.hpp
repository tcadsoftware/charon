
#ifndef CHARON_BCSTRATEGY_FREQDOM_IMPL_HPP
#define CHARON_BCSTRATEGY_FREQDOM_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Sum.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

#include "Charon_FreqDom_Parameters.hpp"
#include <math.h>

// BC Strategies
#include "Charon_BCStrategy_Dirichlet_OhmicContact.hpp"
#include "Charon_BCStrategy_Dirichlet_BJT1DBaseContact.hpp"
#include "Charon_BCStrategy_Dirichlet_ContactOnInsulator.hpp"

#include "Teuchos_ScalarTraits.hpp"

// ***********************************************************************
template <typename EvalT>
charon::BCStrategy_FreqDom<EvalT>::
BCStrategy_FreqDom(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) :
  panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>(bc, global_data)
{
  TEUCHOS_ASSERT( this->m_bc.strategy() == "Frequency Domain" );

  // default
  isLatTDof = false;
  isIonDof = false; 
  isFermiPin = false; 

  // need to modify the td_bc's parameter list for each time collocation point
  this->td_bc = Teuchos::rcp(&bc, false);
  this->td_global_data = global_data;

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_FreqDom<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& user_data)
{
  using Teuchos::RCP;
  using Teuchos::rcp; 
  using Teuchos::ParameterList;
  using std::vector;
  using std::string;
  using std::pair;


  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = side_pb.getParameterList();

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";
  int num_dim = 1;

  // get the freqDomParams pointer
  this->freqDomParamsRCP = eqSetPList.sublist("Options").get<Teuchos::RCP<FreqDomParameters> >("Frequency Domain Parameters");

  // set the fd DOF names required by the time domain BC strategy
  // first, create charon::Names object to handle each harmonic
  this->m_fd_names = std::vector<Teuchos::RCP<charon::Names> >();
  Teuchos::RCP<std::vector<double> > eta = this->freqDomParamsRCP->getRemappedHarmonics();
  for(int i = 0 ; i < this->freqDomParamsRCP->getNumTotalHarmonics() ; i++){
    (this->m_fd_names).emplace_back(Teuchos::rcp(new charon::Names(num_dim,prefix,discfields,discsuffix, "_CosH"+std::to_string((*eta)[i])+"_" )));
    if(i > 0)
      (this->m_fd_names).emplace_back(Teuchos::rcp(new charon::Names(num_dim,prefix,discfields,discsuffix, "_SinH"+std::to_string((*eta)[i])+"_" )));
  }

  // second, determine the time domain BC strategy
  RCP<const Teuchos::ParameterList> dataPList = (this->td_bc)->params();
  std::string time_domain_bc_strategy = dataPList->get<std::string>("Time Domain Strategy");
  bool solveElectrons = true;
  bool solveHoles = true;
  if(time_domain_bc_strategy == "BJT1D Base Contact"){
    std::string baseDopingType = dataPList->get<std::string>("Base Doping Type");
    solveElectrons = (baseDopingType != "Acceptor");
    solveHoles = (baseDopingType != "Donor");
  }
  if(time_domain_bc_strategy == "Contact On Insulator"){
    solveElectrons = false;
    solveHoles = false;
  }

  // finally, conditionally use those charon::Names objects to add the required DOFs
  // as required by the time domain BC strategy
  for(int i = 0 ; i < 2*((this->freqDomParamsRCP)->getNumTotalHarmonics() -1 ) +1 ; i++)
  {
    // get the Data parameter list
    RCP<const Teuchos::ParameterList> dataPList = this->m_bc.params();
    TEUCHOS_ASSERT(!Teuchos::is_null(dataPList));
  
    // check if Fermi Level Pinning is turned on
    if ( dataPList->isParameter("Fermi Level Pinning") )  
      isFermiPin = dataPList->template get<bool>("Fermi Level Pinning"); 

    // set the charon::Names object
    // this loops through both cos and sin dofs (even and odd, resp)
    const charon::Names& n = *(this->m_fd_names[i]); 

    // set dirichlet conditions on all dofs (potential and carrier density)
    const vector<pair<string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator dof_it =
      dofs.begin(); dof_it != dofs.end(); ++dof_it)
    {
// TODO: only add residuals for the BCs requiring them!
      // for the isothermal DD formulations
      if ( (dof_it->first == n.dof.phi) || 
           (dof_it->first == n.dof.edensity && solveElectrons) ||
           (dof_it->first == n.dof.hdensity && solveHoles) )
      {     

        // need the dof values to form the residual
        this->required_dof_names.push_back(dof_it->first);

        // unique residual name
        std::string residual_name = "Residual_" + dof_it->first;

        // map residual to dof
        this->residual_to_dof_names_map[residual_name] = dof_it->first;

        // map residual to target field
        this->residual_to_target_field_map[residual_name] = "Target_" + dof_it->first;

        // For now assume that potential and carrier density use the same basis
        this->basis = dof_it->second;
      }

      // for the DD+Lattice and DD+Lattice+Ion formulations, need to gather T (lattice temperature), since (phi,n,p) depend on T
      if (dof_it->first == n.dof.latt_temp)
      {
        isLatTDof = true; 
        this->required_dof_names.push_back(dof_it->first);
      }
    
      // for the DD+Ion and DD+Lattice+Ion formulations, need to gather ION_DENSITY
      if (dof_it->first == n.dof.iondensity)
      {
        isIonDof = true; 
        this->required_dof_names.push_back(dof_it->first);
      }
    }
  } // end for loop: over harmonics, to map target freqeuncy domain dof names

  //std::cout << "We are using the FreqDom BC Strategy." << std::endl;

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
                             "Error: the name \"" << this->m_bc.equationSetName()
                             << "\" is not a valid DOF for the boundary condition:\n"
                             << this->m_bc << "\n");

  // collect names of evaluated fd dof target fields
  for(int h = 0 ; h < this->freqDomParamsRCP->getNumTotalHarmonics() ; h++)
  {
    this->fd_phi_target_cos_names = Teuchos::rcp(new std::vector<std::string>() );
    this->fd_phi_target_sin_names = Teuchos::rcp(new std::vector<std::string>() );
    this->fd_elec_target_cos_names = Teuchos::rcp(new std::vector<std::string>() );
    this->fd_elec_target_sin_names = Teuchos::rcp(new std::vector<std::string>() );
    this->fd_hole_target_cos_names = Teuchos::rcp(new std::vector<std::string>() );
    this->fd_hole_target_sin_names = Teuchos::rcp(new std::vector<std::string>() );
  }

  for(int h = 0 ; h < this->freqDomParamsRCP->getNumTotalHarmonics() ; h++)
  {
    charon::Names fd_cos_names = charon::Names(num_dim,prefix,discfields,discsuffix, "_CosH"+std::to_string((*eta)[h])+"_");
    // store the names for later use
    this->fd_phi_target_cos_names->emplace_back("Target_"+fd_cos_names.dof.phi);
    if(solveElectrons)
      this->fd_elec_target_cos_names->emplace_back("Target_"+fd_cos_names.dof.edensity);
    if(solveHoles)
      this->fd_hole_target_cos_names->emplace_back("Target_"+fd_cos_names.dof.hdensity);
    charon::Names fd_sin_names = charon::Names(num_dim,prefix,discfields,discsuffix, "_SinH"+std::to_string((*eta)[h])+"_");
    // store the names for later use
    this->fd_phi_target_sin_names->emplace_back("Target_"+fd_sin_names.dof.phi);
    if(solveElectrons)
      this->fd_elec_target_sin_names->emplace_back("Target_"+fd_sin_names.dof.edensity);
    if(solveHoles)
      this->fd_hole_target_sin_names->emplace_back("Target_"+fd_sin_names.dof.hdensity);
  }

  // collect names of the time domain target fields
  this->td_phi_target_names = Teuchos::rcp(new std::vector<std::string>() );
  this->td_elec_target_names = Teuchos::rcp(new std::vector<std::string>() );
  this->td_hole_target_names = Teuchos::rcp(new std::vector<std::string>() );

  for(int tp = 0 ; tp < this->freqDomParamsRCP->getNumTimeCollocationPoints() ; tp++)
  {
    // collect their names, possibly for later use
    charon::Names td_names = charon::Names(num_dim,prefix,discfields,discsuffix, "_TP"+std::to_string(tp)+"_");
    this->td_phi_target_names->emplace_back("Target_"+td_names.dof.phi);
    if(solveElectrons)
      this->td_elec_target_names->emplace_back("Target_"+td_names.dof.edensity);
    if(solveHoles)
      this->td_hole_target_names->emplace_back("Target_"+td_names.dof.hdensity);
  }

}


// ***********************************************************************
template <typename EvalT>
void charon::BCStrategy_FreqDom<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& pb,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
			   const Teuchos::ParameterList& models,
			   const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;

  // build and register all closure models, creating the closure model fields
  // (note that this uses a composite closure model - so we don't need a loop here -
  //  with each component closure model corresponding to a time collocation point)
  pb.buildAndRegisterClosureModelEvaluators(fm,factory,models,user_data);

  // get the element block id and physics id for the given physics block
  string ebID_pb = pb.elementBlockID();
  string pbID = pb.physicsBlockID(); 

  // get the element block ID for the given boundary condition 
  string ebID_bc = this->m_bc.elementBlockID();
  
  if (ebID_pb != ebID_bc) 
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: " << pbID << " corresponds to " 
      << ebID_pb << ", while the BC corresponds to " << ebID_bc << "! \n");
    
  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = pb.getParameterList();
  
  // allow only one equation set per physics block
  if (pbParamList->numParams() > 1)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
       "The physics block " << pbParamList->name() << " has more than one equation sets ! ");

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  //const ParameterList& options = eqSetPList.sublist("Options");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<string>("Discontinuous Suffix") : "";

  Teuchos::RCP<PHX::DataLayout> data_layout = this->basis->functional;

  // get scaling parameters                 
  //Teuchos::RCP<charon::Scaling_Parameters> scaleParams = 
  //  user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");
  //double timeScaling = -1.0*scaleParams->scale_params.t0;
  double pi = 4.0*std::atan(1.0);

  std::string time_domain_bc_strategy;
  bool solveElectrons = true;
  bool solveHoles = true;

  // the tp=0 to tp=L+1 closure model fields have been built above
  // we now use those fields to apply the time domain bc strategy,
  // to evaluate the L+1 dof target fields in the time domain
  for(int tp = 0 ; tp < this->freqDomParamsRCP->getNumTimeCollocationPoints() ; tp++)
  {
    // convention: 
    // where a time domain strategy uses a parameter named Voltage, contained in parameter list 'Data',
    // the frequency domain bc strategy user_data has a parameter *list* named 'Voltage', 
    // with three parameters under it: 'Frequencies', 'Amplitudes', and 'Phase Shifts'

    // we should modify user_data_tp so that it evaluates the parameter Voltage 
    // at time collocation point tp and replaces the list Voltage with the parameter Voltage(tp), i.e., 
    // Voltage(tp) = sum_i amplitudes[i]*sin(2*pi*(remapped)frequencies[i] - phaseshifts[i])

    /* concretely, we have:

            <Parameter name="Strategy" type="string" value="Frequency Domain"/>
            <ParameterList name="Data">
                <Parameter name="Time Domain Strategy" type="string" value="Ohmic Contact"/>
                <ParameterList name="Voltage">
                    <Parameter name="Frequencies"  type="Array(double)" value="{0.0}"/>
                    <Parameter name="Amplitudes"   type="Array(double)" value="{1.0}"/>
                    <Parameter name="Phase Shifts" type="Array(double)" value="{0.0}"/>
                </ParameterList>
            </ParameterList>

    this should be transformed into

            <Parameter name="Strategy" type="string" value="Ohmic Contact"/>
            <ParameterList name="Data">
                <Parameter="Voltage"  type="double"  value=V(tp)>
            </ParameterList>

    where V(tp) is as described above */

    // formula for change:
    // 0) create a new panzer::BC object for this time collocation point,
    //    setting the "Strategy" to "Ohmic Contact" (or whatever is the value of "Time Domain Stragegy")
    // 1) loop over each parameter list L in 'Data'
    // 2) if L is validated to have three parameters 'Frequencies', 'Amplitudes', and 'Phase Shifts',
    //    each of type 'Array(TYPE)', calculate the value of L(tp) of type TYPE
    // 3) add a parameter named L to 'Data', with type TYPE, and value L(tp)
    // note: we retain all other parameters and parameter lists which fail the frequency domain validation
    // 4) use this parameter list at the current time collocation point

    // step 0
    RCP<const Teuchos::ParameterList> dataPList = (this->td_bc)->params();
    time_domain_bc_strategy = dataPList->get<std::string>("Time Domain Strategy");

    if(time_domain_bc_strategy == "BJT1D Base Contact"){
      std::string baseDopingType = dataPList->get<std::string>("Base Doping Type");
      solveElectrons = (baseDopingType == "Donor");
      solveHoles = (baseDopingType == "Acceptor");
    }
    if(time_domain_bc_strategy == "Contact On Insulator"){
      solveElectrons = false;
      solveHoles = false;
    }

    Teuchos::RCP<Teuchos::ParameterList> td_params_tp = Teuchos::parameterList();
    {
      //td_params_tp->set<panzer::BCType>("Type", td_bc->bcType());
      td_params_tp->set<std::string>("Type", "Dirichlet");
      td_params_tp->set<std::string>("Sideset ID", td_bc->sidesetID());
      td_params_tp->set<std::string>("Element Block ID", td_bc->elementBlockID());
      td_params_tp->set<std::string>("Equation Set Name", td_bc->equationSetName());
      td_params_tp->set<std::string>("Strategy", time_domain_bc_strategy);
      if (td_bc->bcType() == panzer::BCType::BCT_Interface) {
	td_params_tp->set<std::string>("Element Block ID2", td_bc->elementBlockID2());
	td_params_tp->set<std::string>("Equation Set Name2", td_bc->equationSetName2());
      }
    }

    // create a new "Data" parameter list, which will have L(tp)
    Teuchos::ParameterList dataPL;

    // step 1 : iterate through the remaining entries of dataPList
    for(ParameterList::ConstIterator dataPList_it = dataPList->begin();
         dataPList_it != dataPList->end(); ++dataPList_it) {
      std::string key = dataPList_it->first;
      // step 2 : check if it is a PARAMETER LIST with three entries 
      // ("Frequencies", "Amplitudes", and "Phase shifts")
      // which must themselves be lists
      if(dataPList_it->second.isList()){
        Teuchos::ParameterList sublist = dataPList->sublist(key);
        if(sublist.isType<Teuchos::Array<double>>("Frequencies")  &&
           sublist.isType<Teuchos::Array<double>>("Amplitudes")   &&
           sublist.isType<Teuchos::Array<double>>("Phase Shifts") ){
          Teuchos::Array<double> frequencies  = 
            sublist.get<Teuchos::Array<double>>("Frequencies");
          Teuchos::Array<double> amplitudes   = 
            sublist.get<Teuchos::Array<double>>("Amplitudes");
          Teuchos::Array<double> phase_shifts = 
            sublist.get<Teuchos::Array<double>>("Phase Shifts");
          // step 3 : calculate the proper value at this time collocation point
          double value = 0.0;
          double timepoint = (this->freqDomParamsRCP->getTimeCollocationPoints())[tp];
          for(int i = 0 ; i < frequencies.length() ; i++)
            value += (frequencies[i] == 0.0 ? amplitudes[i] : 
                      amplitudes[i]*std::sin(frequencies[i]*2.0*pi*timepoint - 2.0*pi*phase_shifts[i]));
            // note here that the argument for sine should be unitless
            // sin(2*pi*\eta_i*tp) is such that \eta_i is unitless and tp is unitless
            // tp = t/(t0 * f0) and tp ranges between 0 and 1 
          dataPL.set<double>(key, value);
          //std::cout << "At time " << std::to_string(tp) << " Voltage is " << value << " for " << ebID_pb << " " << pbID << " " << ebID_bc << std::endl;
        }
        // or, it can belong to a small signal perturbation
        // in which case the entries ("Frequency", "Amplitude", and "Phase Shift")
        // must be present, as doubles
        else if((key == "Small Signal Perturbation") &&
                sublist.isType<double>("Frequency")  && 
                sublist.isType<double>("Amplitude")  &&
                sublist.isType<double>("Phase Shift") ){
          double frequency = sublist.get<double>("Frequency");
          double amplitude = sublist.get<double>("Amplitude");
          double phase_shift = sublist.get<double>("Phase Shift");
          // step 3 : calculate the proper value at this time collocation point
          double timepoint = (this->freqDomParamsRCP->getTimeCollocationPoints())[tp];
          double value = amplitude*std::sin(frequency*2.0*pi*timepoint - 2.0*pi*phase_shift);
            // note here that the argument for sine should be unitless
            // sin(2*pi*\eta_i*tp) is such that \eta_i is unitless and tp is unitless
            // tp = t/(t0 * f0) and tp ranges between 0 and 1 
          dataPL.set<double>("Small Signal Perturbation", value);
          //std::cout << "At tp " << std::to_string(tp) << "/time " << std::to_string(timepoint) << " Small Signal Perturbation is " << value << " for " << ebID_pb << " " << pbID << " " << ebID_bc << std::endl;
        }

      }
      else // if it does not meet this criteria, simply copy it (parameter list or parameter) over
        dataPL.setEntry(key,dataPList_it->second);
    }
    // step 4 : set this parameter list for this time collocation point
    td_params_tp->set("Data", dataPL);

    panzer::BC td_bc_tp(td_bc->bcID(), // has the same bc ID 
                        *td_params_tp, // but a modified 'Data' sublist 
                        this->td_global_data);
    //std::cout << "The time domain strategy specified is: " << td_bc_tp.strategy() << std::endl;
    int num_dim = 1;
    if(time_domain_bc_strategy == "BJT1D Base Contact")
      num_dim = 3;
    // use this modified user_data to set up the time domain bc strategy at time tp
    Teuchos::RCP<panzer::BCStrategy<EvalT> > TD_BC_Strategy_RCP;

    Teuchos::RCP<charon::Names> td_names_tp = Teuchos::rcp(new charon::Names(num_dim,prefix,discfields,discsuffix,"_TP"+std::to_string(tp)+"_"));

    Teuchos::RCP<Teuchos::ParameterList> input_pl = Teuchos::rcp(new Teuchos::ParameterList());
    input_pl->set("Names",td_names_tp);
    input_pl->set("Basis",this->basis);

    if(time_domain_bc_strategy == "Ohmic Contact")
      TD_BC_Strategy_RCP = Teuchos::rcp(new charon::BCStrategy_Dirichlet_OhmicContact<EvalT>(td_bc_tp, this->td_global_data, input_pl) );  
    if(time_domain_bc_strategy == "BJT1D Base Contact")
      TD_BC_Strategy_RCP = Teuchos::rcp(new charon::BCStrategy_Dirichlet_BJT1DBaseContact<EvalT>(td_bc_tp, this->td_global_data, input_pl) );  
    if(time_domain_bc_strategy == "Contact On Insulator")
      TD_BC_Strategy_RCP = Teuchos::rcp(new charon::BCStrategy_Dirichlet_ContactOnInsulator<EvalT>(td_bc_tp, this->td_global_data, input_pl) );  

    // how to add more bc strategies:
    // if(time_domain_bc_strategy == "OTHER STRING")
    //   TD_BC_Strategy_RCP = Teuchos::rcp(new charon::BCStrategy_CORRESPONDING(*td_bc, td_global_data) );  

    // execute the time domain bc strategy's buldAndRegisterEvaluators at time point tp
    TD_BC_Strategy_RCP->buildAndRegisterEvaluators(fm, pb, factory, models, user_data);

  } // end for loop: creation of time domain bc strategy target dof fields

  // perform the discrete Fourier transform of the time domain bc strategy target dof fields
  Teuchos::RCP<std::vector<double> > eta = 
    this->freqDomParamsRCP->getRemappedHarmonics();
  // determine whether performing small or large signal analysis                                                              
  // grab the interwoven DFT coefficients to transform the td dof targets
  // note: for small signal runs, these equal delta_hk (h = harmonic number, k = time index)
  std::vector<std::vector<double> > cos_quad_weights = 
    this->freqDomParamsRCP->getCosQuadratureWeightsForHarmonicNumber();
  std::vector<std::vector<double> > sin_quad_weights = 
    this->freqDomParamsRCP->getSinQuadratureWeightsForHarmonicNumber();

  // take the Fourier transform of the time domain bc strategy target fields
  // (we actually use these to calculate the frequency domain bc's!)
// TODO: when Panzer_Sum supports panzer::EvaluatorStyle::CONTRIBUTES, rewrite these evaluators
  for(int h = 0 ; h < this->freqDomParamsRCP->getNumTotalHarmonics() ; h++)
  {
    // obtain the DFT coefficients to be used by all dofs
    const std::vector<double> const_cos_quad_weights = cos_quad_weights[h];
    Teuchos::RCP<const std::vector<double> > DFT_coeffs_cos_h = 
      Teuchos::rcp(&const_cos_quad_weights,false);
    const std::vector<double> const_sin_quad_weights = sin_quad_weights[h];
    Teuchos::RCP<const std::vector<double> > DFT_coeffs_sin_h = 
      Teuchos::rcp(&const_sin_quad_weights,false);

    // 'Target_ELECTRIC POTENTIAL_CosH'+'h_'
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_phi_target_cos_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_phi_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_cos_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // 'Target_ELECTRIC POTENTIAL_SinH'+'h_'
    if(h>0)
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_phi_target_sin_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_phi_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_sin_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    
    // 'Target_ELECTRON DENSITY_CosH'+'h_'
    if(solveElectrons)
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_elec_target_cos_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_elec_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_cos_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // 'Target_ELECTRON DENSITY_SinH'+'h_'
    if(solveElectrons && h > 0)
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_elec_target_sin_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_elec_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_sin_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // 'Target_HOLE DENSITY_CosH'+'h_'
    if(solveHoles)
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_hole_target_cos_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_hole_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_cos_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // 'Target_HOLE DENSITY_SinH'+'h_'
    if(solveHoles && h > 0)
    {
      Teuchos::ParameterList p;
      p.set("Sum Name",      (*this->fd_hole_target_sin_names)[h]); // this should be a frequency domain dof name
      p.set("Values Names",  this->td_hole_target_names); // this should be a vector of the time domain dof target names
      p.set("Scalars",       DFT_coeffs_sin_h);
      p.set("Data Layout",   data_layout);

      Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
        Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

  } // end of loop over harmonics h

}
#endif
