
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Charon_EquationSet_FreqDom_impl_hpp__
#define __Charon_EquationSet_FreqDom_impl_hpp__

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"

// begin HB mod
// equation set factory for time domain equation set
#include "Panzer_Sum.hpp"

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

#include "Charon_FreqDom_Parameters.hpp"
// end HB mod

// implementation outline
// INPUT DATA: the time domain equation set
//             fundamental frequencies
//             truncation scheme: box, diamond, alpha (for the \ell^\alpha ball)
//             options: transient assisted,
// 0) this equation set should include headers for all FEM time domain equation sets
// 1) call the appropriate time domain equation set, invoking it to create DOFs and fields
//    this may essentially involve running one steady state calculation, saving the fields
//    (note that this solution can be used for the linearized problem,
//    since the steady state is the mode 0 response)
// 2) create the appropriate frequency domain fields from this information,
//    including the HB residual, HB DOF's, and HB Jacobian
// 3) solve the resulting non-linear HB system of equations

// what I need to learn how to do in order to pull this off:
// 1) call an equation set from here, to set up all the relevant fields
// 2) learn how to work with the fields created from that set-up run
// 3) define the HB residual from the residuals defined above
// 3) define the HB Jacobian
// 4) set up the non-linear system and solve it
// 5) revert to the time domain?


// ***********************************************************************

template <typename EvalT>
charon::EquationSet_FreqDom<EvalT>::
EquationSet_FreqDom(const Teuchos::RCP<Teuchos::ParameterList>& params,
		    const int& default_integration_order,
		    const panzer::CellData& cell_data,
		    const Teuchos::RCP<panzer::GlobalData>& global_data,
		    const bool build_transient_support) :
charon::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, 
                        global_data, build_transient_support)
{
  // create vectors to store the names of the frequency domain dofs and their grads
  // this is used to calculate their Fourier transforms, required by panzer::Sum
  Teuchos::RCP<std::vector<std::string> > freqdom_dof_names = Teuchos::rcp(new std::vector<std::string>);
  Teuchos::RCP<std::vector<std::string> > freqdom_grad_dof_names = Teuchos::rcp(new std::vector<std::string>);

  this->td_eqnset_type = params->sublist("Options").get<std::string>("Time Domain Equation Set");

  // "Solve Electron" and "Solve Hole" should only be present in options for the following equation sets
  if((this->td_eqnset_type == "Drift Diffusion") ||
     (this->td_eqnset_type == "SGCVFEM Drift Diffusion")) {
    this->solveElectron = params->sublist("Options").get<std::string>("Solve Electron", "False");
    this->solveHole = params->sublist("Options").get<std::string>("Solve Hole", "False");
  }

  std::string srh_recomb;
  std::string rad_recomb;
  std::string auger_recomb;
  std::string ava_gen;
  std::string defect_cluster_recomb;
  std::string empirical_defect_recomb;
  std::string ionization_particle_strike;
  std::string trap_srh_recomb;
  std::string add_trap_charge;

  if((this->td_eqnset_type != "SGCVFEM Drift Diffusion") &&
     (this->td_eqnset_type != "SGCVFEM Laplace"))
  {
    srh_recomb = params->sublist("Options").get<std::string>("SRH", "Off");
    rad_recomb = params->sublist("Options").get<std::string>("Radiative", "Off");
    auger_recomb = params->sublist("Options").get<std::string>("Auger", "Off");
    ava_gen = params->sublist("Options").get<std::string>("Avalanche", "Off");
    defect_cluster_recomb = params->sublist("Options").get<std::string>("Defect Cluster", "Off");
    empirical_defect_recomb = params->sublist("Options").get<std::string>("Empirical Defect", "Off");
    ionization_particle_strike = params->sublist("Options").get<std::string>("Particle Strike", "Off");

    supg_stab = params->sublist("Options").get<std::string>("SUPG Stabilization", "Off");
    tau_e_type = params->sublist("Options").get<std::string>("Tau_E", "Tanh");
    tau_h_type = params->sublist("Options").get<std::string>("Tau_H", "Tanh");
    ls_type = params->sublist("Options").get<std::string>("Length Scale", "Stream");

    std::string source_stab = params->sublist("Options").get<std::string>("Add Source Term", "Off");
    if (source_stab == "On")
    {
      if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || 
          (ava_gen == "On") || (defect_cluster_recomb == "On") || (trap_srh_recomb == "On") )
        add_source_stab = true;
      else
        add_source_stab = false;
    }
    else
      add_source_stab = false;

    m_timedom_residual_pl.set<std::string>("supg_stab", this->supg_stab);
    m_timedom_residual_pl.set<std::string>("tau_e_type", this->tau_e_type);
    m_timedom_residual_pl.set<std::string>("tau_h_type", this->tau_h_type);
    m_timedom_residual_pl.set<std::string>("ls_type", this->ls_type);
  }

  // Must solve both continuity equations when any recombination model is on.
  if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") ||
      (ava_gen == "On") || (trap_srh_recomb == "On") ||
      (defect_cluster_recomb == "On") || (empirical_defect_recomb == "On") ||
      (ionization_particle_strike == "On"))
  {
    TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
    haveSource = true;
  }
  else
    haveSource = false;  // NO recomb. model

  if(this->td_eqnset_type == "SGCVFEM Drift Diffusion")
    m_timedom_residual_pl.get<std::string>("Driving Force", "EffectiveField");

  if((this->td_eqnset_type == "Drift Diffusion") ||
     (this->td_eqnset_type == "SGCVFEM Drift Diffusion")){
    trap_srh_recomb = params->sublist("Options").get<std::string>("Trap SRH", "Off");
    add_trap_charge = params->sublist("Options").get<std::string>("Add Trap Charge", "False");
    m_timedom_residual_pl.set<bool>("haveSource", haveSource);
    m_timedom_residual_pl.set<bool>("add_source_stab", this->add_source_stab);
    m_timedom_residual_pl.set<bool>("addTrapCharge", this->addTrapCharge);
  }
  // determine if want to add the trap charge to the Poisson equation
  addTrapCharge = (trap_srh_recomb == "On" && add_trap_charge == "True") ? true : false;

  std::string add_fix_charge = params->sublist("Options").get<std::string>("Fixed Charge", "False");
  m_timedom_residual_pl.set<bool>("addFixCharge", this->addFixCharge);
  addFixCharge = (add_fix_charge == "True") ? true : false;
   
  // determine OPTIONS required by the EquationSet_TimeDomain
  m_timedom_residual_pl.set<std::string>("solveElectron", this->solveElectron);
  m_timedom_residual_pl.set<std::string>("solveHole", this->solveHole);

  m_timedom_residual_pl.set<Teuchos::RCP<Teuchos::ParameterList>>("params", Teuchos::rcp(new Teuchos::ParameterList()) );
  m_timedom_residual_pl.set<int>("default_integration_order", default_integration_order);
  m_timedom_residual_pl.set<Teuchos::RCP<panzer::GlobalData>>("global_data", global_data);
  m_timedom_residual_pl.set<bool>("build_transient_support", build_transient_support);
  // "names" also needs to be set (and used by the time domain equation set)
  // this paramter is specified in the bARESE, and loops through the m_td_names vector
  //m_timedom_residual_pl.set<Teuchos::RCP<charon::Names>>("names", this->m_td_names[i]);

  // obtain the RCP<FreqDomParameters> from the main input deck 
  this->freqDomParamsRCP = params->sublist("Options").get<Teuchos::RCP<FreqDomParameters> >("Frequency Domain Parameters");
  //std::cout << "EqnSet_FreqDom was passed a FreqDomParameters pointer from Charon_main; the total number of time collocation points is: " << freqDomParamsRCP->getNumTimeCollocationPoints() << std::endl;
  //std::string fd_suffix = params->sublist("Options").get<std::string>("Frequency Domain Suffix");

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields", "");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix", "");
  std::string basis_type = params->get<std::string>("Basis Type", "");
  int basis_order = params->get<int>("Basis Order", 1);
  std::string model_id = params->get<std::string>("Model ID", "");
  int integration_order = params->get<int>("Integration Order", 1);

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

  // register the DOFs from the time domain equation set
  const Teuchos::RCP<Teuchos::ParameterList> constParams = params;

  if(this->td_eqnset_type == "Laplace")
    this->initializeEquationSet_Laplace(constParams,default_integration_order, 
                                        cell_data, global_data, build_transient_support);
  if(this->td_eqnset_type == "Drift Diffusion")
    this->initializeEquationSet_DriftDiffusion(constParams,default_integration_order, 
                                               cell_data, global_data, build_transient_support);
  if(this->td_eqnset_type == "SGCVFEM Drift Diffusion")
    this->initializeEquationSet_SGCVFEM_DriftDiffusion(constParams,default_integration_order, 
                                                       cell_data, global_data, build_transient_support);
  if(this->td_eqnset_type == "SGCVFEM Laplace")
    this->initializeEquationSet_SGCVFEM_Laplace(constParams,default_integration_order, 
                                                       cell_data, global_data, build_transient_support);

  // *************************************
  // Assemble DOF names and Residual names
  // *************************************

  if (this->buildTransientSupport())
    std::cout << "ERROR: we shouldn't be building with transient support!" << std::endl;

  this->m_td_names = std::vector<Teuchos::RCP<charon::Names> >();
  this->m_fd_names = std::vector<Teuchos::RCP<charon::Names> >();

  // initialize RCP<vector<string>> of FD DOF names, used to construct the TD DOF fields
  this->freqdom_dof_names_phi = Teuchos::rcp(new std::vector<std::string>() );
  this->freqdom_grad_dof_names_phi = Teuchos::rcp(new std::vector<std::string>() );

  if(this->solveElectron == "True")
  {
    this->freqdom_dof_names_elec = Teuchos::rcp(new std::vector<std::string>() );
    this->freqdom_grad_dof_names_elec = Teuchos::rcp(new std::vector<std::string>() );
  }
  if(this->solveHole == "True")
  {
    this->freqdom_dof_names_hole = Teuchos::rcp(new std::vector<std::string>() ); 
    this->freqdom_grad_dof_names_hole = Teuchos::rcp(new std::vector<std::string>() );
  }

// TODO: move this to the EqnSet_FreqDom constructor
// create the frequency domain DOFS
// this produces 2*this->freqDomParamsRCP->getNumTotalHarmonics() fields

Teuchos::RCP<std::vector<double> > eta = this->freqDomParamsRCP->getRemappedHarmonics();
for(int i =0 ; i < this->freqDomParamsRCP->getNumTotalHarmonics() ; i++)
{
  (this->m_fd_names).emplace_back(Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix, "_CosH"+std::to_string((*eta)[i])+"_" )));
  (this->m_fd_names).emplace_back(Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix, "_SinH"+std::to_string((*eta)[i])+"_" )));
  // TEMP: is the following necessary?
  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_fd_names[2*i]));

  // Always solve the Poisson equation (no transients required)
  for(int j = 0 ; j < 2 ; j++){  // j=0 for cosine, j=1 for sine
    if(i > 0 || j == 0){
      this->addDOF(m_fd_names[2*i+j]->dof.phi,basis_type,basis_order,integration_order,m_fd_names[2*i+j]->res.phi);
      this->addDOFGrad(m_fd_names[2*i+j]->dof.phi,m_fd_names[2*i+j]->grad_dof.phi);
      this->freqdom_dof_names_phi->emplace_back(m_fd_names[2*i+j]->dof.phi);
      this->freqdom_grad_dof_names_phi->emplace_back(m_fd_names[2*i+j]->grad_dof.phi);
    }
  }

  // Solve the Electron equation by default and users can disable it 
  if (solveElectron == "True"){
    for(int j = 0 ; j < 2 ; j++){  // j=0 for cosine, j=1 for sine
      if(i > 0 || j == 0){
        this->addDOF(m_fd_names[2*i+j]->dof.edensity,basis_type,basis_order,integration_order,m_fd_names[2*i+j]->res.edensity);
        this->addDOFGrad(m_fd_names[2*i+j]->dof.edensity,m_fd_names[2*i+j]->grad_dof.edensity);
        this->freqdom_dof_names_elec->emplace_back(m_fd_names[2*i+j]->dof.edensity);
        this->freqdom_grad_dof_names_elec->emplace_back(m_fd_names[2*i+j]->grad_dof.edensity);
      }
    }
  }
  
  // Solve the Hole equation by default and users can disable it 
  if (solveHole == "True"){
    for(int j = 0 ; j < 2 ; j++){  // j=0 for cosine, j=1 for sine
      if(i > 0 || j == 0){
        this->addDOF(m_fd_names[2*i+j]->dof.hdensity,basis_type,basis_order,integration_order,m_fd_names[2*i+j]->res.hdensity);
        this->addDOFGrad(m_fd_names[2*i+j]->dof.hdensity,m_fd_names[2*i+j]->grad_dof.hdensity);
        this->freqdom_dof_names_hole->emplace_back(m_fd_names[2*i+j]->dof.hdensity);
        this->freqdom_grad_dof_names_hole->emplace_back(m_fd_names[2*i+j]->grad_dof.hdensity);
      }
    }
  }

 } // end registering fd DOFs

  this->addClosureModel(model_id);
  this->setupDOFs();

  // create the vector of names used to create the TD DOFs
  this->m_td_names = std::vector<Teuchos::RCP<charon::Names> >();
  for(int time_coll_pt = 0 ; time_coll_pt < this->freqDomParamsRCP->getNumTimeCollocationPoints() ; time_coll_pt++)
      (this->m_td_names).emplace_back(Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix, "_TP"+std::to_string(time_coll_pt)+"_" )));

} // END constructor

// ********************************************************************

// BEGIN initializeEquationSet_Laplace

template <typename EvalT>
void charon::EquationSet_FreqDom<EvalT>::
initializeEquationSet_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
                              const int& default_integration_order,
                              const panzer::CellData& cell_data,
                              const Teuchos::RCP<panzer::GlobalData>& global_data,
                              const bool build_transient_support)
{
  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);
    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    // HB Mod: time domain equation set specification
    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");
    // this is separate from the FreqDomParameters input pl so that we can simply use that pl for instantiation of the FreqDomParameters object
    opt.set("Time Domain Equation Set","","Time domain equation set to be analyzed in the frequency domain");

    // HB Mod: validate the precense of a RCP<FreqDomParameters> parameter, created by Charon_Main
    opt.set("Frequency Domain Parameters", Teuchos::rcp(new FreqDomParameters()) ,"(For internal use: set by Charon_main.cpp)");

    // HB Mod: validate the the parameters specified to create the FreqDomParameters object from Charon_Main
    Teuchos::ParameterList& freqdom_opt = valid_parameters.sublist("Frequency Domain Options");  
    // the isSmallSignal boolean
    freqdom_opt.set("Enable Small Signal Analysis", false ,"Set to true for SS, or to false for LS.");
    // the truncation scheme
    Teuchos::setStringToIntegralParameter<int>(
      "Truncation Scheme",
      "Box", // default gives an invalid type
      "Choose the truncation scheme for the harmonic balance method.",
      Teuchos::tuple<std::string>("Box", "Diamond", "Hybrid"),
      &freqdom_opt
      );
    freqdom_opt.set("Hybrid Exponent",0.5,"Choose the hybrid truncation scheme exponent 0 < \alpha < 1.");
    // the truncation order
    freqdom_opt.set("Truncation Order",3,"Choose the truncation order of the harmonic balance method.");
    // the list of fundamental harmonics
    Teuchos::Array<double> default_harmonics;
    freqdom_opt.set("Fundamental Harmonics", default_harmonics,"Choose the fundamental harmonics of the system.");
    freqdom_opt.set("Remapped Fundamental Harmonics", default_harmonics,"Optionally, choose the remapped fundamental harmonics of the system.");
    // specify the number of time collocation points, defaulting to the Nyquist Sampling Theorem required minimum of 0
    freqdom_opt.set("Number of Time Collocation Points", 0, "Choose the integer number of time collocation points; choosing 0 defaults to the minimum required by the Nyquist Sampling Theorem.");

// these are replicated from the drift diffusion equation set, and are essentially unused
    // validate time domain equation set parameters
    Teuchos::setStringToIntegralParameter<int>(
      "Fermi Dirac",
      "False",
      "Determine if users want to use the Fermi-Dirac statistics",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "FD Formula",
      "Schroeder",  // default
      "Determine if users want to use Schroeder's or the Diffusion method for FD impl.",
      Teuchos::tuple<std::string>("Schroeder","Diffusion"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determine if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Trap SRH",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determine if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determine if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determine if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Defect Cluster",
      "Off",
      "Determine if users want to include Defect Cluster generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Particle Strike",
      "Off",
      "Determine if users want to include ionization due to a particle strike",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Empirical Defect",
      "Off",
      "Determine if users want to include Empirical Defect generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SUPG Stabilization",
      "Off",
      "Enable or disable SUPG stabilization in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_E",
      "None",
      "Determine the Tau_E model in the Electron equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_H",
      "None",
      "Determines the Tau_H model in the Hole equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Length Scale",
      "Stream",
      "Determine the Length Scale model in the SUPG stabilization",
      Teuchos::tuple<std::string>("Stream","Shakib"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Source Term",
      "Off",
      "Determine if users want to add the source term in the SUPG stabilization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determine if users want to turn on the BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Optical Generation",
      "Off",
      "Determine if users want to turn on the optical generation",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Acceptor Incomplete Ionization",
      "Off",
      "Determine if users want to include Acceptor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Donor Incomplete Ionization",
      "Off",
      "Determine if users want to include Donor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Trap SRH",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Trap Charge",
      "False",
      "Determine if users want to include the srh trap charge in the Poisson equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
// end borrowed parameters from drift diffusion

    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    // determine if a user wants to include bulk fixed charges
    addFixCharge = false; 
    if (params->sublist("Options").get<std::string>("Fixed Charge") == "True")
      addFixCharge = true; 
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  //int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  //int integration_order = params->get<int>("Integration Order");

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

/* TODO: Let FreqDom handle this!
  // ********************
  // Assemble DOF names and Residual names
  // ********************
  m_names = Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix));

  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_names));

  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);
  if (this->buildTransientSupport())
    this->addDOFTimeDerivative(m_names->dof.phi,m_names->dxdt.phi);

  this->addClosureModel(model_id);

  this->setupDOFs();
*/
}

// END initializeEquationSet_Laplace


// ********************************************************************

// BEGIN initializeEquationSet_DriftDiffusion
template <typename EvalT>
void charon::EquationSet_FreqDom<EvalT>::
initializeEquationSet_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
		    const int& default_integration_order,
		    const panzer::CellData& cell_data,
		    const Teuchos::RCP<panzer::GlobalData>& global_data,
		    const bool build_transient_support)
{
  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;

    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");

    // HB Mod: time domain equation set specification
    // this is separate from the FreqDomParameters input pl so that we can simply use that pl for instantiation of the FreqDomParameters object
    opt.set("Time Domain Equation Set","","Time domain equation set to be analyzed in the frequency domain");

    // HB Mod: validate the precense of a RCP<FreqDomParameters> parameter, created by Charon_Main
    opt.set("Frequency Domain Parameters", Teuchos::rcp(new FreqDomParameters()) ,"(For internal use: set by Charon_main.cpp)");

    // HB Mod: validate the the parameters specified to create the FreqDomParameters object from Charon_Main
    Teuchos::ParameterList& freqdom_opt = valid_parameters.sublist("Frequency Domain Options");  
    // the isSmallSignal boolean
    freqdom_opt.set("Enable Small Signal Analysis", false ,"Set to true for SS, or to false for LS.");
    // the truncation scheme
    Teuchos::setStringToIntegralParameter<int>(
      "Truncation Scheme",
      "Box", // default gives an invalid type
      "Choose the truncation scheme for the harmonic balance method.",
      Teuchos::tuple<std::string>("Box", "Diamond", "Hybrid"),
      &freqdom_opt
      );
    freqdom_opt.set("Hybrid Exponent",0.5,"Choose the hybrid truncation scheme exponent 0 < \alpha < 1.");
    // the truncation order
    freqdom_opt.set("Truncation Order",3,"Choose the truncation order of the harmonic balance method.");
    // the list of fundamental harmonics
    Teuchos::Array<double> default_harmonics;
    freqdom_opt.set("Fundamental Harmonics", default_harmonics,"Choose the fundamental harmonics of the system.");
    freqdom_opt.set("Remapped Fundamental Harmonics", default_harmonics,"Optionally, choose the remapped fundamental harmonics of the system.");
    // specify the number of time collocation points, defaulting to the Nyquist Sampling Theorem required minimum of 0
    freqdom_opt.set("Number of Time Collocation Points", 0, "Choose the integer number of time collocation points; choosing 0 defaults to the minimum required by the Nyquist Sampling Theorem.");

    // validate time domain equation set parameters
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Electron",
      "True",
      "Determine if users want to solve the Electron Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Hole",
      "True",
      "Determine if users want to solve the Hole Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Fermi Dirac",
      "False",
      "Determine if users want to use the Fermi-Dirac statistics",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "FD Formula",
      "Schroeder",  // default
      "Determine if users want to use Schroeder's or the Diffusion method for FD impl.",
      Teuchos::tuple<std::string>("Schroeder","Diffusion"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determine if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determine if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determine if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determine if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band2Band Tunneling",
      "Off",
      "Determines if users want to include Band2Band tunneling generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On", "Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Driving Force",
      "EffectiveField",
      "Determines driving force for different models in the continuity equation(s)",
      Teuchos::tuple<std::string>("EffectiveField","GradQuasiFermi","GradPotential"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Defect Cluster",
      "Off",
      "Determine if users want to include Defect Cluster generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Particle Strike",
      "Off",
      "Determine if users want to include ionization due to a particle strike",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Empirical Defect",
      "Off",
      "Determine if users want to include Empirical Defect generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SUPG Stabilization",
      "On",
      "Enable or disable SUPG stabilization in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_E",
      "Tanh",
      "Determine the Tau_E model in the Electron equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Tau_H",
      "Tanh",
      "Determines the Tau_H model in the Hole equation",
      Teuchos::tuple<std::string>("None","Linear","Tanh"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Length Scale",
      "Stream",
      "Determine the Length Scale model in the SUPG stabilization",
      Teuchos::tuple<std::string>("Stream","Shakib"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Source Term",
      "Off",
      "Determine if users want to add the source term in the SUPG stabilization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determine if users want to turn on the BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Optical Generation",
      "Off",
      "Determine if users want to turn on the optical generation",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Acceptor Incomplete Ionization",
      "Off",
      "Determine if users want to include Acceptor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Donor Incomplete Ionization",
      "Off",
      "Determine if users want to include Donor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Trap SRH",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Dynamic Traps",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Trap Charge",
      "False",
      "Determine if users want to include the srh trap charge in the Poisson equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");
    
    // Must solve at least one continuity equation
    //if ((solveElectron == "False") && (solveHole == "False"))
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: must solve at lease one"
    //       << " continuity equation !" << "\n");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH"); 
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");  
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");  
    const std::string& defect_cluster_recomb = params->sublist("Options").get<std::string>("Defect Cluster");  
    const std::string& empirical_defect_recomb = params->sublist("Options").get<std::string>("Empirical Defect");  
    const std::string& ionization_particle_strike = params->sublist("Options").get<std::string>("Particle Strike");  
    const std::string& trap_srh_recomb = params->sublist("Options").get<std::string>("Trap SRH");
    const std::string& add_trap_charge = params->sublist("Options").get<std::string>("Add Trap Charge");

    // determine if want to add the trap charge to the Poisson equation
    addTrapCharge = (trap_srh_recomb == "On" && add_trap_charge == "True") ? true : false;

    // Must solve both continuity equations when any recombination model is on. 
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || 
        (ava_gen == "On") || (trap_srh_recomb == "On") ||
        (defect_cluster_recomb == "On") || (empirical_defect_recomb == "On") ||
	(ionization_particle_strike == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
      haveSource = true; 
    }
    else
      haveSource = false;  // NO recomb. model  
    
    supg_stab = params->sublist("Options").get<std::string>("SUPG Stabilization", "Off");
    tau_e_type = params->sublist("Options").get<std::string>("Tau_E", "Tanh");
    tau_h_type = params->sublist("Options").get<std::string>("Tau_H", "Tanh");
    ls_type = params->sublist("Options").get<std::string>("Length Scale", "Stream");
    
    std::string source_stab = params->sublist("Options").get<std::string>("Add Source Term", "Off");
    if (source_stab == "On")
    {
      if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || 
          (ava_gen == "On") || (defect_cluster_recomb == "On") || (trap_srh_recomb == "On") )
        add_source_stab = true;
      else
        add_source_stab = false;
    }
    else
      add_source_stab = false;

    // error check for Incomplete Ionization models
    const std::string& acc_ioniz = 
      params->sublist("Options").get<std::string>("Acceptor Incomplete Ionization");
    const std::string& don_ioniz = 
      params->sublist("Options").get<std::string>("Donor Incomplete Ionization");
    if(acc_ioniz == "On" && solveHole == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
	  "Error: must solve hole continuity equation when " << 
	  "Acceptor Incomplete Ionization is On!" << "\n");
    if(don_ioniz == "On" && solveElectron == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
	  "Error: must solve electron continuity equation when " <<
	  "Donor Incomplete Ionization is On!" << "\n");
  }
  
}
// END initializeDriftDiffusion

// ***********************************************************************

// BEGIN initializeEquationSet_SGCVFEM_DriftDiffusion

template <typename EvalT>
void charon::EquationSet_FreqDom<EvalT>::
initializeEquationSet_SGCVFEM_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                             const int& default_integration_order,
                                             const panzer::CellData& cell_data,
                                             const Teuchos::RCP<panzer::GlobalData>& global_data,
                                             const bool build_transient_support)
{
  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;

    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");

    // HB block: time domain equation set specification
    // this is separate from the FreqDomParameters input pl so that we can simply use that pl for instantiation of the FreqDomParameters object
    opt.set("Time Domain Equation Set","","Time domain equation set to be analyzed in the frequency domain");
    // validate the precense of a RCP<FreqDomParameters> parameter, created by Charon_Main
    opt.set("Frequency Domain Parameters", Teuchos::rcp(new FreqDomParameters()) ,"(For internal use: set by Charon_main.cpp)");
    // validate the the parameters specified to create the FreqDomParameters object from Charon_Main
    Teuchos::ParameterList& freqdom_opt = valid_parameters.sublist("Frequency Domain Options");  
    freqdom_opt.set("Enable Small Signal Analysis", false ,"Set to true for SS, or to false for LS.");
    Teuchos::setStringToIntegralParameter<int>(
      "Truncation Scheme",
      "Box", // default gives an invalid type
      "Choose the truncation scheme for the harmonic balance method.",
      Teuchos::tuple<std::string>("Box", "Diamond", "Hybrid"),
      &freqdom_opt
      );
    freqdom_opt.set("Hybrid Exponent",0.5,"Choose the hybrid truncation scheme exponent 0 < \alpha < 1.");
    freqdom_opt.set("Truncation Order",3,"Choose the truncation order of the harmonic balance method.");
    Teuchos::Array<double> default_harmonics;
    freqdom_opt.set("Fundamental Harmonics", default_harmonics,"Choose the fundamental harmonics of the system.");
    freqdom_opt.set("Remapped Fundamental Harmonics", default_harmonics,"Optionally, choose the remapped fundamental harmonics of the system.");
    freqdom_opt.set("Number of Time Collocation Points", 0, "Choose the integer number of time collocation points; choosing 0 defaults to the minimum required by the Nyquist Sampling Theorem.");

    // validate time domain equation set parameters (copied from EquationSet_SGCVFEM_DriftDiffusion)
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Electron",
      "True",
      "Determines if users want to solve the Electron Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Hole",
      "True",
      "Determines if users want to solve the Hole Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Fermi Dirac",
      "False",
      "Determine if users want to use the Fermi-Dirac statistics",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "FD Formula",
      "Schroeder",  // default
      "Determine if users want to use Schroeder's or the Diffusion method for FD impl.",
      Teuchos::tuple<std::string>("Schroeder","Diffusion"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determines if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determines if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determines if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determines if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Driving Force",
      "EffectiveField",
      "Determines driving force for different models in the continuity equation(s)",
      Teuchos::tuple<std::string>("EffectiveField","GradQuasiFermi","GradPotential"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band2Band Tunneling",
      "Off",
      "Determines if users want to include Band2Band tunneling generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On", "Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Defect Cluster",
      "Off",
      "Determines if users want to include Defect Cluster generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Particle Strike",
      "Off",
      "Determine if users want to include ionization due to a particle strike",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Empirical Defect",
      "Off",
      "Determines if users want to include Empirical Defect generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determines if users want to turn on BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Optical Generation",
      "Off",
      "Determine if users want to turn on the optical generation",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Acceptor Incomplete Ionization",
      "Off",
      "Determine if users want to include Acceptor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Donor Incomplete Ionization",
      "Off",
      "Determine if users want to include Donor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Trap SRH",
      "Off",
      "Determine if users want to include Trap SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Add Trap Charge",
      "False",
      "Determine if users want to include the trap charge in the Poisson equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      ); // ACH: This is being validated, but not used...
    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");

    // Must solve at least one continuity equation
    //if ((solveElectron == "False") && (solveHole == "False"))
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: must solve at least one"
    //       << " continuity equation !" << "\n");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH");
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");
    const std::string& defect_cluster_recomb = params->sublist("Options").get<std::string>("Defect Cluster");
    const std::string& empirical_defect_recomb = params->sublist("Options").get<std::string>("Empirical Defect");
    const std::string& ionization_particle_strike = params->sublist("Options").get<std::string>("Particle Strike");
    const std::string& trap_srh_recomb = params->sublist("Options").get<std::string>("Trap SRH");
    const std::string& add_trap_charge = params->sublist("Options").get<std::string>("Add Trap Charge");
    const std::string& opt_gen = params->sublist("Options").get<std::string>("Optical Generation"); 

    // Must solve both continuity equations when any recomb./gen. model is on.
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") ||
        (trap_srh_recomb == "On") || (ava_gen == "On") || (opt_gen == "On") || 
        (defect_cluster_recomb == "On") || (empirical_defect_recomb == "On") ||
        (ionization_particle_strike == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
      haveSource = true;
    }
    else
      haveSource = false;  // NO recomb. model

    withAvaGen = (haveSource && ava_gen == "On") ? true : false;
    withTrapSRH = (haveSource && trap_srh_recomb == "On") ? true : false;
    addTrapCharge = (trap_srh_recomb == "On" && add_trap_charge == "True") ? true : false;
    drForce = params->sublist("Options").get<std::string>("Driving Force");

    // error check for Incomplete Ionization models
    const std::string& acc_ioniz =
      params->sublist("Options").get<std::string>("Acceptor Incomplete Ionization");
    const std::string& don_ioniz =
      params->sublist("Options").get<std::string>("Donor Incomplete Ionization");
    if(acc_ioniz == "On" && solveHole == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve hole continuity equation when " <<
          "Acceptor Incomplete Ionization is On!" << "\n");
    if(don_ioniz == "On" && solveElectron == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve electron continuity equation when " <<
          "Donor Incomplete Ionization is On!" << "\n");
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  std::string model_id = params->get<std::string>("Model ID");

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

}

// END initializeEquationSet_SGCVFEM_DriftDiffusion

// ***********************************************************************

// BEGIN initializeEquationSet_SGCVFEM_Laplace

template <typename EvalT>
void charon::EquationSet_FreqDom<EvalT>::
initializeEquationSet_SGCVFEM_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                      const int& default_integration_order,
                                      const panzer::CellData& cell_data,
                                      const Teuchos::RCP<panzer::GlobalData>& global_data,
                                      const bool build_transient_support)
{
  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;

    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    //valid_parameters.sublist("Options");
    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");

    // HB block: time domain equation set specification
    // this is separate from the FreqDomParameters input pl so that we can simply use that pl for instantiation of the FreqDomParameters object
    opt.set("Time Domain Equation Set","","Time domain equation set to be analyzed in the frequency domain");
    // validate the precense of a RCP<FreqDomParameters> parameter, created by Charon_Main
    opt.set("Frequency Domain Parameters", Teuchos::rcp(new FreqDomParameters()) ,"(For internal use: set by Charon_main.cpp)");

    // HB block: time domain equation set specification
    // this is separate from the FreqDomParameters input pl so that we can simply use that pl for instantiation of the FreqDomParameters object
    opt.set("Time Domain Equation Set","","Time domain equation set to be analyzed in the frequency domain");
    // validate the precense of a RCP<FreqDomParameters> parameter, created by Charon_Main
    opt.set("Frequency Domain Parameters", Teuchos::rcp(new FreqDomParameters()) ,"(For internal use: set by Charon_main.cpp)");
    // validate the the parameters specified to create the FreqDomParameters object from Charon_Main
    Teuchos::ParameterList& freqdom_opt = valid_parameters.sublist("Frequency Domain Options");  
    freqdom_opt.set("Enable Small Signal Analysis", false ,"Set to true for SS, or to false for LS.");
    Teuchos::setStringToIntegralParameter<int>(
      "Truncation Scheme",
      "Box", // default gives an invalid type
      "Choose the truncation scheme for the harmonic balance method.",
      Teuchos::tuple<std::string>("Box", "Diamond", "Hybrid"),
      &freqdom_opt
      );
    freqdom_opt.set("Hybrid Exponent",0.5,"Choose the hybrid truncation scheme exponent 0 < \alpha < 1.");
    freqdom_opt.set("Truncation Order",3,"Choose the truncation order of the harmonic balance method.");
    Teuchos::Array<double> default_harmonics;
    freqdom_opt.set("Fundamental Harmonics", default_harmonics,"Choose the fundamental harmonics of the system.");
    freqdom_opt.set("Remapped Fundamental Harmonics", default_harmonics,"Optionally, choose the remapped fundamental harmonics of the system.");
    freqdom_opt.set("Number of Time Collocation Points", 0, "Choose the integer number of time collocation points; choosing 0 defaults to the minimum required by the Nyquist Sampling Theorem.");

    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    // determine if a user wants to include bulk fixed charges
    addFixCharge = false; 
    if (params->sublist("Options").get<std::string>("Fixed Charge") == "True")
      addFixCharge = true; 
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  std::string model_id = params->get<std::string>("Model ID");

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

}

// END initializeEquationSet_SGCVFEM_Laplace

// ***********************************************************************

template <typename EvalT>
void charon::EquationSet_FreqDom<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& fl,
				      const Teuchos::ParameterList& user_data) const
{

  // For now, assume that all dofs use the same basis
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_fd_names[0]->dof.phi); 
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_fd_names[0]->dof.phi); 

  // take the Fourier transform of the frequency domain DOFs, 
  // producing this->freqDomParamsRCP->getNumTimeCollocationPoints() fields
  // these are DOF(t_i) for each t_i a time collocation point
  // and are used as the time domain equation sets' DOFs, but aren't *actual* DOFs!

  // get the required coefficients to calculate the time domain DOFs
  std::vector<std::vector<double> > interwoven_DFT_coeffs_at_time = this->freqDomParamsRCP->getDofInterwovenCoeffsAtTimeCollPt();
  // interwoven_DFT_coeffs_at_time[i], with i a time collocation point number, contains the coefficients of
  // u(t_i) = U_0 + (U_1^c) cos(\omega_1 t_i) + (U_1^s) sin(\omega_1 t_i) + ..
  // that is, interwoven_DFT_coeffs_at_time[i] = {0, cos(\omega_1 t_i) , sin(\omega_1 t_i) , ...}

  //std::cout << "We are now taking the Fourier transforms of the fd DOFs to produce the td DOFS" << std::endl;

for(int time_coll_pt = 0 ; time_coll_pt < this->freqDomParamsRCP->getNumTimeCollocationPoints() ; time_coll_pt++)
{
      // an RCP pointing to the vector<double> of this time collocation point's DFT coefficients
      const std::vector<double> const_interwoven_DFT_coeffs_at_tp = interwoven_DFT_coeffs_at_time[time_coll_pt];
      Teuchos::RCP<const std::vector<double> >  interwoven_DFT_coeffs = 
        Teuchos::rcp(&const_interwoven_DFT_coeffs_at_tp,false);

      //std::cout << "The interwoven DFT coefficients are: ";
      //for(auto item : *interwoven_DFT_coeffs)
      //  std::cout << std::to_string(item) << " ";
      //std:: cout << std::endl;

      // TODO: check whether cell_data.num_dim(), discfuelds, discsuffix

  // Always solve the Poisson equation (no transients required)
      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at nodes
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.phi);
        p.set("Values Names",  this->freqdom_dof_names_phi);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   basis->functional);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
	    Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }

      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.phi);
        p.set("Values Names",  this->freqdom_dof_names_phi);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_scalar);

       Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
           Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
      
      // Add the field named "GRAD"+dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->grad_dof.phi);
        p.set("Values Names",  this->freqdom_grad_dof_names_phi);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_vector);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
            Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }

  // Solve the Electron equation by default and users can disable it 
  if (solveElectron == "True") 
  {
      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at nodes
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.edensity);
        p.set("Values Names",  this->freqdom_dof_names_elec);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   basis->functional);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
            Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }

      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.edensity);
        p.set("Values Names",  this->freqdom_dof_names_elec);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_scalar);

       Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
           Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
      
      // Add the field named "GRAD"+dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->grad_dof.edensity);
        p.set("Values Names",  this->freqdom_grad_dof_names_elec);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_vector);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
            Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
  }
  
  // Solve the Hole equation by default and users can disable it 
  if (solveHole == "True") 
  {
      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at nodes
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.hdensity);
        p.set("Values Names",  this->freqdom_dof_names_hole);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   basis->functional);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
            Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }

      // Add the field named dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->dof.hdensity);
        p.set("Values Names",  this->freqdom_dof_names_hole);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_scalar);

       Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
           Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
      
      // Add the field named "GRAD"+dof_at_time_coll_pt to the field manager
      // Supported at IPs
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      this->m_td_names[time_coll_pt]->grad_dof.hdensity);
        p.set("Values Names",  this->freqdom_grad_dof_names_hole);
        p.set("Scalars",       interwoven_DFT_coeffs);
        p.set("Data Layout",   ir->dl_vector);

  	Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op = 
            Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
  }

}
// end taking DOFs of fd DOFs to produce td DOFs, with appropriate suffixes determined by m_td_names

  // execute the time domain equation set bARESE for each time collocation point
  // (NOT for each harmonic; e.g, note num harmonics < num time collocation points)
  std::vector<Teuchos::ParameterList> timeDomPL = {};
  std::vector<Teuchos::RCP<charon::EquationSet_DefaultImpl<EvalT> > > eqnSet_TimeDomain = {};
  std::vector<std::string> td_phi_residual_contributions = {};
  std::vector<std::string> td_elec_residual_contributions = {};
  std::vector<std::string> td_hole_residual_contributions = {};
  for(int i = 0 ; i < this->freqDomParamsRCP->getNumTimeCollocationPoints() ; i++)
  {
      timeDomPL.push_back(this->m_timedom_residual_pl);
      timeDomPL[i].set<Teuchos::RCP<charon::Names>>("Names", (this->m_td_names)[i]);

      if(this->td_eqnset_type == "Laplace")
        eqnSet_TimeDomain.push_back(Teuchos::rcp(new charon::EquationSet_Laplace<EvalT>(ir, basis, timeDomPL[i]) ) );
      if(this->td_eqnset_type == "Drift Diffusion")
        eqnSet_TimeDomain.push_back(Teuchos::rcp(new charon::EquationSet_DriftDiffusion<EvalT>(ir, basis, timeDomPL[i]) ) );
      if(this->td_eqnset_type == "SGCVFEM Drift Diffusion")
        eqnSet_TimeDomain.push_back(Teuchos::rcp(new charon::EquationSet_SGCVFEM_DriftDiffusion<EvalT>(ir, basis, timeDomPL[i]) ) );
      if(this->td_eqnset_type == "SGCVFEM Laplace")
        eqnSet_TimeDomain.push_back(Teuchos::rcp(new charon::EquationSet_SGCVFEM_Laplace<EvalT>(ir, basis, timeDomPL[i]) ) );

      eqnSet_TimeDomain[i]->buildAndRegisterEquationSetEvaluators(fm, fl, user_data);

      td_phi_residual_contributions.push_back(m_td_names[i]->res.phi);
      if((this->td_eqnset_type == "Drift Diffusion") || 
         (this->td_eqnset_type == "SGCVFEM Drift Diffusion"))
      {
        td_elec_residual_contributions.push_back(m_td_names[i]->res.edensity);
        td_hole_residual_contributions.push_back(m_td_names[i]->res.hdensity);
      }
  }

// build the HB residual equations (via the Discrete Fourier Transform)
  std::vector<double> truncated_harmonic_basis_unremapped =
    *((this->freqDomParamsRCP)->getUnRemappedHarmonics());
  std::vector<double> truncated_harmonic_basis =
    *((this->freqDomParamsRCP)->getRemappedHarmonics());
  // determine whether performing small or large signal analysis
  bool isSmallSignal = this->freqDomParamsRCP->queryIsSmallSignal();
  // get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = 
    user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");
  double timeScaling = -2.0*scaleParams->scale_params.t0;
  double pi = 4.0*std::atan(1.0);

// TEST HUGE FREQUENCY SCALING
//double huge_frequency_scaling = 1.0;///truncated_harmonic_basis_unremapped[1];

  std::vector<std::vector<double> > cos_quad_weights = 
    (this->freqDomParamsRCP)->getCosQuadratureWeightsForHarmonicNumber();
  std::vector<std::vector<double> > sin_quad_weights = 
    (this->freqDomParamsRCP)->getSinQuadratureWeightsForHarmonicNumber();
  // large signal: cos_quad_weights[i], with i a harmonic basis number, 
  //               contains the coefficients of \Sum_k w_k^i R(t_k)
  //               = (i==0 ? 1 : 2) * \int_[0,T] TimeDom_Residual(t) * cos(\omega_i t) dt
  // small signal: cos_quad_weights[i], with i a harmonic basis number, 
  //               contains the coefficients of \Sum_k w_k^i R(t_k)
  //               = (i==k ? R(t_k) : 0)

// TEST HUGE FREQUENCY SCALING
//for(int i = 1 ; i < cos_quad_weights.size() ; i++)
//  for(auto j : cos_quad_weights[i])
//    j *= huge_frequency_scaling;// /double(i);
//for(int i = 1 ; i < sin_quad_weights.size() ; i++)
//  for(auto j : sin_quad_weights[i])
//    j *= huge_frequency_scaling;// /double(i);

// phi cos(h) and sin(h) HB residual equations
  // phi cosine(h0)
  this->buildAndRegisterResidualSummationEvaluator(fm, 
    (*(this->freqdom_dof_names_phi))[0] /* residual equation's name, phi cos(h) */,
  	td_phi_residual_contributions,
    cos_quad_weights[0]);

  for(int h = 1 ; h < (this->freqDomParamsRCP)->getNumTotalHarmonics() ; h++){

      // phi cosine(harmonic h)
      this->buildAndRegisterResidualSummationEvaluator(fm, 
        (*(this->freqdom_dof_names_phi))[2*h-1] /* residual equation's name, phi cos(h) */,
	td_phi_residual_contributions,
        cos_quad_weights[h]);
      // phi sin(harmonic h)
      this->buildAndRegisterResidualSummationEvaluator(fm, 
        (*(this->freqdom_dof_names_phi))[2*h] /* residual equation's name, phi sin(h) */,
        td_phi_residual_contributions,
        sin_quad_weights[h]);
  } // end the harmonic loop
// end building the phi HB equations

// elec cos(h) and sin(h) HB residual equations
if(solveElectron == "True"){
  // elec cos(h0)
  this->buildAndRegisterResidualSummationEvaluator(fm, 
    (*(this->freqdom_dof_names_elec))[0] /* residual equation's name, elec cos(h) */,
    td_elec_residual_contributions,
    cos_quad_weights[0]);
  // elec cos(h>0) and sin(h>0)
  for(int h = 1 ; h < (this->freqDomParamsRCP)->getNumTotalHarmonics() ; h++){
      // elec cos(harmonic h)
      this->buildAndRegisterResidualSummationEvaluator(fm, 
        (*(this->freqdom_dof_names_elec))[2*h-1] /* residual equation's name, elec cos(h) */,
	td_elec_residual_contributions,
        cos_quad_weights[h]);
      // add the d/dt term for the h^th cosine harmonic
      if(h>0)
      {
        // the derivative of cos(w*t) is -w*sin(w*t), but d/dt is on the other side of the equation
        // small signal: no factor other than w is present
        // large signal: the 0.5 factor comes from \int cos(2*\pi*\omega*x)^2 from x=0 to x=1 equals 0.5
        double multiplier = ( isSmallSignal ?
	    2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h] :
	    0.5*2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h]);

// TEST HUGE FREQUENCY SCALING
//multiplier *= huge_frequency_scaling;// /double(h);

	Teuchos::RCP< PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
          Teuchos::rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(
	    panzer::EvaluatorStyle::CONTRIBUTES, 
            (this->m_fd_names[2*h+1])->res.edensity,
            (this->m_fd_names[2*h])->dof.edensity,
            *basis, *ir, multiplier));
        this->template registerEvaluator<EvalT>(fm, op);
      }

      // elec sin(harmonic h)
      if(h > 0)
        this->buildAndRegisterResidualSummationEvaluator(fm, 
          (*(this->freqdom_dof_names_elec))[2*h] /* residual equation's name, elec sin(h) */,
  	  td_elec_residual_contributions,
          sin_quad_weights[h]);

      // add the d/dt term for the h^th sine harmonic
      if(h>0)
      {
        // the derivative of sin(w*t) is w*cos(w*t), but d/dt is on the other side of the equation
        // small signal: no factor other than -w is present
        // large signal: the 0.5 factor comes from \int sin(2*\pi*\omega*x)^2 from x=0 to x=1 equals 0.5
        double multiplier = ( isSmallSignal ?
	    -2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h] :
	    -0.5*2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h]);

// TEST HUGE FREQUENCY SCALING
//multiplier *= huge_frequency_scaling;// /double(h);

        Teuchos::RCP< PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
          Teuchos::rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(
	    panzer::EvaluatorStyle::CONTRIBUTES, 
	    (this->m_fd_names[2*h])->res.edensity,
            (this->m_fd_names[2*h+1])->dof.edensity,
            *basis, *ir, multiplier));
        this->template registerEvaluator<EvalT>(fm, op);
      }
  } // end the harmonic loop
}
// end building the elec HB equations


// hole cos(h) and sin(h) HB residual equations
if(solveHole == "True"){
  this->buildAndRegisterResidualSummationEvaluator(fm, 
    (*(this->freqdom_dof_names_hole))[0] /* residual equation's name, hole cos(h) */,
    td_hole_residual_contributions,
    cos_quad_weights[0]);
  for(int h = 1 ; h < (this->freqDomParamsRCP)->getNumTotalHarmonics() ; h++){

      // hole cos(harmonic h)
      this->buildAndRegisterResidualSummationEvaluator(fm, 
        (*(this->freqdom_dof_names_hole))[2*h-1] /* residual equation's name, hole cos(h) */,
	td_hole_residual_contributions,
	cos_quad_weights[h]);

      // add the d/dt term for the h^th cosine harmonic
      if(h>0)
      {
        // the derivative of cos(w*t) is -w*sin(w*t), but d/dt is on the other side of the equation
        // small signal: no factor other than w is present
        // large signal: the 0.5 factor comes from \int cos(2*\pi*\omega*x)^2 from x=0 to x=1 equals 0.5
        double multiplier = ( isSmallSignal ?
	    2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h] :
	    0.5*2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h]);

// TEST HUGE FREQUENCY SCALING
//multiplier *= huge_frequency_scaling;// /double(h);

	Teuchos::RCP< PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
          Teuchos::rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(
	    panzer::EvaluatorStyle::CONTRIBUTES, 
            (this->m_fd_names[2*h+1])->res.hdensity,
            (this->m_fd_names[2*h])->dof.hdensity,
            *basis, *ir, multiplier));
        this->template registerEvaluator<EvalT>(fm, op);
      }

      // hole sin(harmonic h)
      if(h > 0)
        this->buildAndRegisterResidualSummationEvaluator(fm, 
          (*(this->freqdom_dof_names_hole))[2*h] /* residual equation's name, phi sin(h) */,
  	  td_hole_residual_contributions,
	  sin_quad_weights[h]);

      // add the d/dt term for the h^th sine harmonic
      if(h>0)
      {
        // the derivative of sin(w*t) is w*cos(w*t), but d/dt is on the other side of the equation
        // small signal: no factor other than -w is present
        // large signal: the 0.5 factor comes from \int sin(2*\pi*\omega*x)^2 from x=0 to x=1 equals 0.5
        double multiplier = ( isSmallSignal ?
	    -2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h] :
	    -0.5*2.0*pi*timeScaling*truncated_harmonic_basis_unremapped[h]);

// TEST HUGE FREQUENCY SCALING
//multiplier *= huge_frequency_scaling;// /double(h);

        Teuchos::RCP< PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
          Teuchos::rcp(new panzer::Integrator_BasisTimesScalar<EvalT,panzer::Traits>(
	    panzer::EvaluatorStyle::CONTRIBUTES, 
	    (this->m_fd_names[2*h])->res.hdensity,
            (this->m_fd_names[2*h+1])->dof.hdensity,
            *basis, *ir, multiplier));
        this->template registerEvaluator<EvalT>(fm, op);
      }
  } // end the harmonic loop
}
// end building the hole HB equations


/*
      // enforce temporal C^1 continuity 
      // candidate: (phi_tp0 - phi_tpL)^2 + ((tp1-tp0)/num_tp - (tp(L)-tp(L-1))/num_tp)^2
      //            =  (phi_tp0) * (phi_tp0) - 2 * (phi_tp0) * (phi_tpL) - (phi_tpL) * (phi_tpL)
      //             + 1/(num_tp)^2 * (phi_tp1 - phi_tp0 + phi_tp(L-1) - phi_tp(L)) ^2
      {
        Teuchos::ParameterList p;
        p.set("Sum Name",      (*this->fd_phi_target_sin_names)[h]);
        p.set("Values Names",  Teuchos::rcp(new std::vector<std::string>{
                               (*this->fd_phi_target_sin_names)[h],
                               (*(this->m_fd_names[2*h+0])).dof.phi} ));
        p.set("Scalars", Teuchos::rcp(new const std::vector<double>{1.0, -1.0}));
        p.set("Data Layout",   data_layout);

        Teuchos::RCP<PHX::EvaluatorWithBaseImpl<panzer::Traits> > op =
          Teuchos::rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, op);
      }
*/

}

// ***********************************************************************

#endif
