// @HEADER
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

#ifndef __Charon_EquationSet_FreqDom_decl_hpp__
#define __Charon_EquationSet_FreqDom_decl_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Charon_EquationSet_DefaultImpl.hpp"

#include "Panzer_Traits.hpp"
#include "Phalanx_FieldManager.hpp"
 
#ifndef __Charon_FreqDom_Parameters_cpp__
#include "Charon_FreqDom_Parameters.cpp"
#endif

namespace charon {

template <typename EvalT>
class EquationSet_FreqDom : public charon::EquationSet_DefaultImpl<EvalT> {

public:    

  EquationSet_FreqDom(const Teuchos::RCP<Teuchos::ParameterList>& params,
                      const int& default_integration_order,
                      const panzer::CellData& cell_data,
                      const Teuchos::RCP<panzer::GlobalData>& gd,
                      const bool build_transient_support);

  void initializeEquationSet_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                     const int& default_integration_order,
                                     const panzer::CellData& cell_data,
                                     const Teuchos::RCP<panzer::GlobalData>& global_data,
                                     const bool build_transient_support);

  void initializeEquationSet_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                            const int& default_integration_order,
                                            const panzer::CellData& cell_data,
                                            const Teuchos::RCP<panzer::GlobalData>& global_data,
                                            const bool build_transient_support);

  void initializeEquationSet_SGCVFEM_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                                    const int& default_integration_order,
                                                    const panzer::CellData& cell_data,
                                                    const Teuchos::RCP<panzer::GlobalData>& global_data,
                                                    const bool build_transient_support);

  void initializeEquationSet_SGCVFEM_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
                                             const int& default_integration_order,
                                             const panzer::CellData& cell_data,
                                             const Teuchos::RCP<panzer::GlobalData>& global_data,
                                             const bool build_transient_support);

  void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                             const panzer::FieldLibrary& field_library,
                                             const Teuchos::ParameterList& user_data) const;

  // begin HB mod
  // add evaluators from the Helmholtz equation set
  //  user_app::EquationSet_Helmholtz<EvalT> equationSet_Helmholtz(const bool build_transient_support, 
  //                                                           std::string input_dof_name_, 
  //                                                           Teuchos::RCP<panzer::IntegrationRule> input_ir, 
  //                                                           Teuchos::RCP<panzer::BasisIRLayout> input_basis);
  // end HB mod

private:

  std::string td_eqnset_type;

  int num_fundamental_harmonics;
  int num_total_harmonics;
  int num_time_collocation_points;

  Teuchos::RCP<Teuchos::ParameterList> freqdom_pl;

  std::string dof_name_;

  std::vector<Teuchos::RCP<charon::Names> >  m_td_names;
  std::vector<Teuchos::RCP<charon::Names> >  m_fd_names;

  Teuchos::RCP<std::vector<std::string> > freqdom_dof_names_phi;
  Teuchos::RCP<std::vector<std::string> > freqdom_grad_dof_names_phi;
  Teuchos::RCP<std::vector<std::string> > freqdom_dof_names_elec;
  Teuchos::RCP<std::vector<std::string> > freqdom_grad_dof_names_elec;
  Teuchos::RCP<std::vector<std::string> > freqdom_dof_names_hole;
  Teuchos::RCP<std::vector<std::string> > freqdom_grad_dof_names_hole;

  Teuchos::ParameterList timedom_residual_pl;

  // Options
  std::string solveElectron;
  std::string solveHole;

  std::string supg_stab;
  std::string tau_e_type;
  std::string tau_h_type;
  std::string ls_type;

  bool haveSource = false;
  bool add_source_stab = false; 
  bool addTrapCharge = false; 
  std::string drForce;

  bool withAvaGen = false;
  bool withTrapSRH = false;
  bool addFixCharge = false; 
  
  // for frequency domain simulations
  std::string m_dof_name_;
  Teuchos::RCP<panzer::IntegrationRule> m_ir;
  Teuchos::RCP<panzer::BasisIRLayout> m_basis;    
  Teuchos::ParameterList m_timedom_residual_pl;

  Teuchos::RCP<FreqDomParameters> freqDomParamsRCP;

};

}

#include "Charon_EquationSet_FreqDom_impl.hpp"

#endif
