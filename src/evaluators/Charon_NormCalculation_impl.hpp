#ifndef CHARON_NORM_CALCULATION_IMPL_HPP
#define CHARON_NORM_CALCULATION_IMPL_HPP

#include <string>
#include <cmath>

#include "Kokkos_ViewFactory.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IosAllSaver.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp" //For exo output

#include "Intrepid2_FunctionSpaceTools.hpp"



namespace charon {


//////////////////////////////////////////////////
//                                             ///
//                  L2 Error                   ///
//                                             ///
//////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Norm_L2Error<EvalT, Traits>::
Norm_L2Error(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using panzer::Dim;

  

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  //std::cout << "input params:" << std::endl;
  //std::cout << params << std::endl;
  //params.validateParameters(*valid_params);


  string name = params.get<string>("Name");
  comm = params.get<RCP<Teuchos::Comm<int> const> >("Comm");


  // Obtain BASIS information
  RCP<BasisIRLayout> basis = params.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;  
  numNodes = data_layout->dimension(1);
  basisName = basis->name();


  // Grab names of Fields
  //string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  string const errorPrefix = params.get< std::string >("Error Prefix");


  // Obtain IP information
  RCP<IntegrationRule const> ir = params.get<RCP<IntegrationRule> >("IR");
  quadOrder = ir->cubature_degree;
  numIPs = (ir->dl_scalar)->dimension(1);


  // Grab layouts for quadrature points on mesh
  //RCP<DataLayout> manufactured_layout_vec = ir->dl_vector; //NEW ADDITION
  RCP<DataLayout> manufactured_layout = ir->dl_scalar; //NEW ADDITION


  // Number of dimensions of gradient
  //num_dims = manufactured_layout_vec->dimension(2); //NEW ADDITION


  // Set values based on names of fields created in closure model factory
  simulation_solution = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);
  manufactured_solution = PHX::MDField<const ScalarT,Cell,Point>("analytic_"+name, manufactured_layout);
  

  // Add fields for DAG addDependentField lets Trilinos know to evaluate those fields before the addEvaluatedFields ones 
  this->addDependentField(simulation_solution);
  this->addDependentField(manufactured_solution);
  this->addEvaluatedField(error);

  outputNormName = simulation_solution.fieldTag().name()+"_L2_Error_Norm";
  outputNorm =
    panzer::createAndRegisterScalarParameter<EvalT>(outputNormName, *params.sublist("Norm ParameterList").get<RCP<panzer::ParamLib> >("ParamLib"));

  string n = "Norm Calculation: L2 Error " + name;
  this->setName(n);
}


template<typename EvalT, typename Traits>
void
Norm_L2Error<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd, 
  PHX::FieldManager<Traits>& /* fm */)
{
  quadIndex = panzer::getIntegrationRuleIndex(quadOrder,(*sd.worksets_)[0]);
  basisIndex = panzer::getBasisIndex(basisName,(*sd.worksets_)[0]);
  //CHANGED
  integral = Kokkos::createDynRankView(simulation_solution.get_static_view(), "integral", simulation_solution.dimension(0));
}



template<typename EvalT, typename Traits>
void
Norm_L2Error<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;


  // FieldContainer to temporarily hold L2 error at IPs
  //CHANGED
  Kokkos::DynRankView<ScalarT,PHX::Device> error_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "error_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "simulation_solution_ip", workset.num_cells, numIPs);
  //Kokkos::DynRankView<ScalarT,PHX::Device> manufactured_solution_para = Kokkos::createDynRankView(manufactured_solution.get_static_view(), "manufactured_solution_para", workset.num_cells, numIPs);
  Kokkos::DynRankView<const ScalarT,PHX::Device> manufactured_solution_para = manufactured_solution.get_view();
  
  // Zero-out arrays for 
  Kokkos::deep_copy(integral, ScalarT(0.0));
  Kokkos::deep_copy(error_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_ip, ScalarT(0.0));
  //Kokkos::deep_copy(manufactured_solution_para, ScalarT(0.0));

  // Convert the simulation solution at BASIS points to IPs and compute the pointwise error at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int ip = 0; ip < numIPs; ++ip)
    {   
      for (int node = 0; node < numNodes; ++node)
      {   
          simulation_solution_ip(cell,ip) += (workset.bases[basisIndex])->basis_scalar(cell,node,ip) * simulation_solution(cell,node);
      }   
      // Compute error at IPs and square it
      //manufactured_solution_para(cell,ip) = manufactured_solution(cell,ip); // MANUALLY COPYING manufactured solution to a parallelized object. This seems dumb.
      error_ip(cell,ip) = simulation_solution_ip(cell,ip) - manufactured_solution_para(cell,ip);
      error_ip(cell,ip) *= error_ip(cell,ip);
   
    }   
  }


  // Numerically integrate error_ip
  if (workset.num_cells > 0)
  {
     // Compute integral on each cell  
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::integrate(integral, error_ip, (workset.int_rules[quadIndex])->weighted_measure.get_view());
  
     // Combine integrals on cells
     for(index_t i = 0; i<workset.num_cells; i++)
     {
        l2error += integral(i);
     }
  }
} 
  

template <typename EvalT,typename Traits>
void Norm_L2Error<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
   // Initialize accumulated error
   l2error = Teuchos::ScalarTraits<ScalarT>::zero();
}

template <typename EvalT,typename Traits>
Teuchos::RCP<Teuchos::ParameterList> Norm_L2Error<EvalT,Traits>::getValidParameters() const
{
  
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Value", "string");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> integration_rule;
  p->set("IR", integration_rule);
  //Teuchos::RCP<charon::Names> n;
  p->set<std::string>("Name", "string");
  p->set<std::string>("Analytic Prefix", "string");
  p->set<std::string>("Error Prefix", "string");
  Teuchos::RCP<Teuchos::Comm<int> const> communicator;
  p->set("Comm", communicator);
  p->sublist("Norm ParameterList");
  p->sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));

  return p;
}

template <typename EvalT,typename Traits>
void Norm_L2Error<EvalT,Traits>::postEvaluate(typename Traits::PostEvalData /* d */)
{
  this->postprocess(std::cout);
  outputNorm->setValue(std::sqrt(globalL2Error));
}


template<typename EvalT, typename Traits>
void Norm_L2Error<EvalT, Traits>::postprocess(std::ostream& os) 
{
  // throw unless specialized for residual evaluations
  os << "WEIRD! It is a bad idea to use this evaluation type for Norm_L2Error!" << std::endl;
}


template<>
void Norm_L2Error<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{

  // Sum the error across processors and output it to the screen.
  globalL2Error = 0.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &l2error, &globalL2Error);

  if (comm->getRank() == 0) {
     panzer::ios_all_saver saver(os);

     const std::string outputStr = "L2 Error MMS "+simulation_solution.fieldTag().name();
     std::size_t precision = 8;
     std::size_t name_width = outputStr.size();
     std::size_t value_width = precision + 7;

     os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
     os << std::setw(name_width) << outputStr
        << " " << std::setw(value_width) << std::sqrt(globalL2Error)
        << std::endl;
  }
}


//////////////////////////////////////////////////
//                                             ///
//       L2 Norm of Simulation Solution        ///
//                                             ///
//////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Norm_L2<EvalT, Traits>::
Norm_L2(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using panzer::Dim;

  

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  //std::cout << "input params:" << std::endl;
  //std::cout << params << std::endl;
  //params.validateParameters(*valid_params);


  string name = params.get<string>("Name");
  comm = params.get<RCP<Teuchos::Comm<int> const> >("Comm");


  // Obtain BASIS information
  RCP<BasisIRLayout> basis = params.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;  
  numNodes = data_layout->dimension(1);
  basisName = basis->name();


  // Grab names of Fields
  //string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  string const errorPrefix = params.get< std::string >("Error Prefix");


  // Obtain IP information
  RCP<IntegrationRule const> ir = params.get<RCP<IntegrationRule> >("IR");
  quadOrder = ir->cubature_degree;
  numIPs = (ir->dl_scalar)->dimension(1);


  // Grab layouts for quadrature points on mesh
  //RCP<DataLayout> manufactured_layout_vec = ir->dl_vector; //NEW ADDITION
  RCP<DataLayout> manufactured_layout = ir->dl_scalar; //NEW ADDITION


  // Number of dimensions of gradient
  //num_dims = manufactured_layout_vec->dimension(2); //NEW ADDITION


  // Set values based on names of fields created in closure model factory
  simulation_solution = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);
  

  // Add fields for DAG addDependentField lets Trilinos know to evaluate those fields before the addEvaluatedFields ones 
  this->addDependentField(simulation_solution);
  this->addEvaluatedField(error);

  outputNormName = simulation_solution.fieldTag().name()+"_L2_Norm";
  outputNorm =
    panzer::createAndRegisterScalarParameter<EvalT>(outputNormName, *params.sublist("Norm ParameterList").get<RCP<panzer::ParamLib> >("ParamLib"));

  string n = "Norm Calculation: L2 " + name;
  this->setName(n);
}


template<typename EvalT, typename Traits>
void
Norm_L2<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd, 
  PHX::FieldManager<Traits>& /* fm */)
{
  quadIndex = panzer::getIntegrationRuleIndex(quadOrder,(*sd.worksets_)[0]);
  basisIndex = panzer::getBasisIndex(basisName,(*sd.worksets_)[0]);
  //CHANGED
  integral = Kokkos::createDynRankView(simulation_solution.get_static_view(), "integral", simulation_solution.dimension(0));
}



template<typename EvalT, typename Traits>
void
Norm_L2<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;


  // FieldContainer to temporarily hold L2 error at IPs
  //CHANGED
  Kokkos::DynRankView<ScalarT,PHX::Device> norm_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "norm_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "simulation_solution_ip", workset.num_cells, numIPs);

  
  // Zero-out arrays for 
  Kokkos::deep_copy(integral, ScalarT(0.0));
  Kokkos::deep_copy(norm_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_ip, ScalarT(0.0));
  

  // Convert the simulation solution at BASIS points to IPs and compute the pointwise error at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int ip = 0; ip < numIPs; ++ip)
    {   
      for (int node = 0; node < numNodes; ++node)
      {   
          simulation_solution_ip(cell,ip) += (workset.bases[basisIndex])->basis_scalar(cell,node,ip) * simulation_solution(cell,node);
      }   
      // Compute error at IPs and square it
      norm_ip(cell,ip) = simulation_solution_ip(cell,ip);
      norm_ip(cell,ip) *= norm_ip(cell,ip);
   
    }   
  }


  // Numerically integrate error_ip
  if (workset.num_cells > 0)
  {
     // Compute integral on each cell  
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::integrate(integral, norm_ip, (workset.int_rules[quadIndex])->weighted_measure.get_view());
  
     // Combine integrals on cells
     for(index_t i = 0; i<workset.num_cells; i++)
     {
        l2norm += integral(i);
     }
  }
} 
  

template <typename EvalT,typename Traits>
void Norm_L2<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
   // Initialize accumulated error
   l2norm = Teuchos::ScalarTraits<ScalarT>::zero();
}

template <typename EvalT,typename Traits>
Teuchos::RCP<Teuchos::ParameterList> Norm_L2<EvalT,Traits>::getValidParameters() const
{
  
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Value", "string");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> integration_rule;
  p->set("IR", integration_rule);
  //Teuchos::RCP<charon::Names> n;
  p->set<std::string>("Name", "string");
  p->set<std::string>("Analytic Prefix", "string");
  p->set<std::string>("Error Prefix", "string");
  Teuchos::RCP<Teuchos::Comm<int> const> communicator;
  p->set("Comm", communicator);
  p->sublist("Norm ParameterList");
  p->sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));

  return p;
}

template <typename EvalT,typename Traits>
void Norm_L2<EvalT,Traits>::postEvaluate(typename Traits::PostEvalData /* d */)
{
  this->postprocess(std::cout);
  outputNorm->setValue(std::sqrt(globalL2Norm));
}


template<typename EvalT, typename Traits>
void Norm_L2<EvalT, Traits>::postprocess(std::ostream& os) 
{
  // throw unless specialized for residual evaluations
  os << "WEIRD! It is a bad idea to use this evaluation type for Norm_L2Error!" << std::endl;
}


template<>
void Norm_L2<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{

  // Sum the error across processors and output it to the screen.
  globalL2Norm = 0.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &l2norm, &globalL2Norm);

  if (comm->getRank() == 0) {
     panzer::ios_all_saver saver(os);

     const std::string outputStr = "L2 Norm "+simulation_solution.fieldTag().name();
     std::size_t precision = 8;
     std::size_t name_width = outputStr.size();
     std::size_t value_width = precision + 7;

     os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
     os << std::setw(name_width) << outputStr
        << " " << std::setw(value_width) << std::sqrt(globalL2Norm)
        << std::endl;
  }
}

//////////////////////////////////////////////////
//                                             ///
//                  H1 Norm Error              ///
//                                             ///
//////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Norm_H1Error<EvalT, Traits>::
Norm_H1Error(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using panzer::Dim;

  

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  //std::cout << "input params:" << std::endl;
  //std::cout << params << std::endl;
  //params.validateParameters(*valid_params);


  string name = params.get<string>("Name");
  comm = params.get<RCP<Teuchos::Comm<int> const> >("Comm");


  // Obtain BASIS information
  RCP<BasisIRLayout> basis = params.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;  
  numNodes = data_layout->dimension(1);
  basisName = basis->name();


  // Grab names of Fields
  //string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  string const errorPrefix = params.get< std::string >("Error Prefix");


  // Obtain IP information
  RCP<IntegrationRule const> ir = params.get<RCP<IntegrationRule> >("IR");
  quadOrder = ir->cubature_degree;
  numIPs = (ir->dl_scalar)->dimension(1);


  // Grab layouts for quadrature points on mesh
  RCP<DataLayout> manufactured_layout_vec = ir->dl_vector; //NEW ADDITION
  RCP<DataLayout> manufactured_layout = ir->dl_scalar; //NEW ADDITION


  // Number of dimensions of gradient
  num_dims = manufactured_layout_vec->dimension(2); //NEW ADDITION
  

  // Set values based on names of fields created in closure model factory
  simulation_solution = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);
  manufactured_solution = PHX::MDField<const ScalarT,Cell,Point>("analytic_"+name, manufactured_layout);
  manufactured_gradient = PHX::MDField<const ScalarT,Cell,Point,Dim>("GRAD_analytic_"+name, manufactured_layout_vec);
   

  // Add fields for DAG addDependentField lets Trilinos know to evaluate those fields before the addEvaluatedFields ones 
  this->addDependentField(simulation_solution);
  this->addDependentField(manufactured_solution);
  this->addDependentField(manufactured_gradient);
  this->addEvaluatedField(error);

  outputNormName = simulation_solution.fieldTag().name()+"_H1_Error_Norm";
  outputNorm =
    panzer::createAndRegisterScalarParameter<EvalT>(outputNormName, *params.sublist("Norm ParameterList").get<RCP<panzer::ParamLib> >("ParamLib"));

  string n = "Norm Calculation: H1 Error " + name;
  this->setName(n);
}


template<typename EvalT, typename Traits>
void
Norm_H1Error<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd, 
  PHX::FieldManager<Traits>& /* fm */)
{
  quadIndex = panzer::getIntegrationRuleIndex(quadOrder,(*sd.worksets_)[0]);
  basisIndex = panzer::getBasisIndex(basisName,(*sd.worksets_)[0]);
  //CHANGED
  integral = Kokkos::createDynRankView(simulation_solution.get_static_view(), "integral", simulation_solution.dimension(0));
}



template<typename EvalT, typename Traits>
void
Norm_H1Error<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold L2 error at IPs
  Kokkos::DynRankView<ScalarT,PHX::Device> error_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "error_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "simulation_solution_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_grad_ip = Kokkos::createDynRankView(manufactured_gradient.get_static_view(), "simulation_solution_grad_ip", workset.num_cells, numIPs, num_dims);
  Kokkos::DynRankView<const ScalarT,PHX::Device> manufactured_solution_para = manufactured_solution.get_view();
  Kokkos::DynRankView<const ScalarT,PHX::Device> manufactured_solution_grad_para = manufactured_gradient.get_view();
  //Kokkos::DynRankView<ScalarT,PHX::Device> manufactured_solution_grad_para = Kokkos::createDynRankView(manufactured_gradient.get_static_view(), "simulation_solution_grad_ip", workset.num_cells, numIPs, num_dims);

  // Zero-out arrays for 
  Kokkos::deep_copy(integral, ScalarT(0.0));
  Kokkos::deep_copy(error_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_grad_ip, ScalarT(0.0));
  //Kokkos::deep_copy(manufactured_solution_grad_para, ScalarT(0.0));


  

  // Convert the simulation solution at BASIS points to IPs and compute the pointwise error at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int ip = 0; ip < numIPs; ++ip)
    {   
      for (int node = 0; node < numNodes; ++node)
      {   
          simulation_solution_ip(cell,ip) += (workset.bases[basisIndex])->basis_scalar(cell,node,ip) * simulation_solution(cell,node);
          for (int dim = 0; dim < num_dims; ++dim)
          {
            simulation_solution_grad_ip(cell,ip,dim) += (workset.bases[basisIndex])->grad_basis(cell,node,ip,dim) * simulation_solution(cell,node);
          }
      }   


      // Compute error at IPs and square it
      error_ip(cell,ip) = simulation_solution_ip(cell,ip) - manufactured_solution_para(cell,ip);
      error_ip(cell,ip) *= error_ip(cell,ip);
      // Add in gradient contributions to error
      for (int dim = 0; dim < num_dims; ++dim)
      {
        error_ip(cell,ip) += (simulation_solution_grad_ip(cell,ip,dim) - manufactured_solution_grad_para(cell,ip,dim)) * (simulation_solution_grad_ip(cell,ip,dim) - manufactured_solution_grad_para(cell,ip,dim));
      }
    }   
  }


  // Numerically integrate error_ip
  if (workset.num_cells > 0)
  {
     // Compute integral on each cell  
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::integrate(integral, error_ip, (workset.int_rules[quadIndex])->weighted_measure.get_view());
  
     // Combine integrals on cells
     for(index_t i = 0; i<workset.num_cells; i++)
     {
        h1error += integral(i);
     }
  }
} 
  

template <typename EvalT,typename Traits>
void Norm_H1Error<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
   // Initialize accumulated error
   h1error = Teuchos::ScalarTraits<ScalarT>::zero();
}

template <typename EvalT,typename Traits>
Teuchos::RCP<Teuchos::ParameterList> Norm_H1Error<EvalT,Traits>::getValidParameters() const
{
  
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Value", "string");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> integration_rule;
  p->set("IR", integration_rule);
  //Teuchos::RCP<charon::Names> n;
  p->set<std::string>("Name", "string");
  p->set<std::string>("Analytic Prefix", "string");
  p->set<std::string>("Error Prefix", "string");
  Teuchos::RCP<Teuchos::Comm<int> const> communicator;
  p->set("Comm", communicator);
  p->sublist("Norm ParameterList");
  p->sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));

  return p;
}

template <typename EvalT,typename Traits>
void Norm_H1Error<EvalT,Traits>::postEvaluate(typename Traits::PostEvalData /* d */)
{
  this->postprocess(std::cout);
  outputNorm->setValue(std::sqrt(globalH1Error));
}


template<typename EvalT, typename Traits>
void Norm_H1Error<EvalT, Traits>::postprocess(std::ostream& os) 
{
  // throw unless specialized for residual evaluations
  os << "WEIRD! It is a bad idea to use this evaluation type for Norm_H1Error!" << std::endl;
}


template<>
void Norm_H1Error<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{
  // Sum the error across processors and output it to the screen.
  globalH1Error = 0.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &h1error, &globalH1Error);

  if (comm->getRank() == 0) {
     panzer::ios_all_saver saver(os);

     const std::string outputStr = "H1 Error MMS "+simulation_solution.fieldTag().name();
     std::size_t precision = 8;
     std::size_t name_width = outputStr.size();
     std::size_t value_width = precision + 7;

     os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
     os << std::setw(name_width) << outputStr
        << " " << std::setw(value_width) << std::sqrt(globalH1Error)
        << std::endl;
  }
}



//////////////////////////////////////////////////
//                                             ///
//       H1 Norm of Simulation Solution        ///
//                                             ///
//////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Norm_H1<EvalT, Traits>::
Norm_H1(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using panzer::Dim;

  

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  //std::cout << "input params:" << std::endl;
  //std::cout << params << std::endl;
  //params.validateParameters(*valid_params);


  string name = params.get<string>("Name");
  comm = params.get<RCP<Teuchos::Comm<int> const> >("Comm");


  // Obtain BASIS information
  RCP<BasisIRLayout> basis = params.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;  
  numNodes = data_layout->dimension(1);
  basisName = basis->name();


  // Grab names of Fields
  //string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  string const errorPrefix = params.get< std::string >("Error Prefix");


  // Obtain IP information
  RCP<IntegrationRule const> ir = params.get<RCP<IntegrationRule> >("IR");
  quadOrder = ir->cubature_degree;
  numIPs = (ir->dl_scalar)->dimension(1);


  // Grab layouts for quadrature points on mesh
  RCP<DataLayout> manufactured_layout_vec = ir->dl_vector; //NEW ADDITION
  RCP<DataLayout> manufactured_layout = ir->dl_scalar; //NEW ADDITION


  // Number of dimensions of gradient
  num_dims = manufactured_layout_vec->dimension(2); //NEW ADDITION
  

  // Set values based on names of fields created in closure model factory
  simulation_solution = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  manufactured_gradient = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_analytic_solution", manufactured_layout_vec);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);
  

  // Add fields for DAG addDependentField lets Trilinos know to evaluate those fields before the addEvaluatedFields ones 
  this->addDependentField(simulation_solution);
  this->addEvaluatedField(error);

  outputNormName = simulation_solution.fieldTag().name()+"_H1_Norm";
  outputNorm =
    panzer::createAndRegisterScalarParameter<EvalT>(outputNormName, *params.sublist("Norm ParameterList").get<RCP<panzer::ParamLib> >("ParamLib"));

  string n = "Norm Calculation: H1 " + name;
  this->setName(n);
}


template<typename EvalT, typename Traits>
void
Norm_H1<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd, 
  PHX::FieldManager<Traits>& /* fm */)
{
  quadIndex = panzer::getIntegrationRuleIndex(quadOrder,(*sd.worksets_)[0]);
  basisIndex = panzer::getBasisIndex(basisName,(*sd.worksets_)[0]);
  //CHANGED
  integral = Kokkos::createDynRankView(simulation_solution.get_static_view(), "integral", simulation_solution.dimension(0));
}



template<typename EvalT, typename Traits>
void
Norm_H1<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold L2 error at IPs
  Kokkos::DynRankView<ScalarT,PHX::Device> norm_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "norm_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_ip = Kokkos::createDynRankView(simulation_solution.get_static_view(), "simulation_solution_ip", workset.num_cells, numIPs);
  Kokkos::DynRankView<ScalarT,PHX::Device> simulation_solution_grad_ip = Kokkos::createDynRankView(manufactured_gradient.get_static_view(), "simulation_solution_grad_ip", workset.num_cells, numIPs, num_dims);


  // Zero-out arrays for 
  Kokkos::deep_copy(integral, ScalarT(0.0));
  Kokkos::deep_copy(norm_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_ip, ScalarT(0.0));
  Kokkos::deep_copy(simulation_solution_grad_ip, ScalarT(0.0));
  

  // Convert the simulation solution at BASIS points to IPs and compute the pointwise error at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int ip = 0; ip < numIPs; ++ip)
    {   
      for (int node = 0; node < numNodes; ++node)
      {   
          simulation_solution_ip(cell,ip) += (workset.bases[basisIndex])->basis_scalar(cell,node,ip) * simulation_solution(cell,node);
          for (int dim = 0; dim < num_dims; ++dim)
          {
            simulation_solution_grad_ip(cell,ip,dim) += (workset.bases[basisIndex])->grad_basis(cell,node,ip,dim) * simulation_solution(cell,node);
          }
      }   


      // Compute norm at IPs and square it
      norm_ip(cell,ip) = simulation_solution_ip(cell,ip);
      norm_ip(cell,ip) *= norm_ip(cell,ip);
      // Add in gradient contributions to error
      for (int dim = 0; dim < num_dims; ++dim)
      {
        norm_ip(cell,ip) += simulation_solution_grad_ip(cell,ip,dim) * simulation_solution_grad_ip(cell,ip,dim);
      }
    }   
  }


  // Numerically integrate error_ip
  if (workset.num_cells > 0)
  {
     // Compute integral on each cell  
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::integrate(integral, norm_ip, (workset.int_rules[quadIndex])->weighted_measure.get_view());
  
     // Combine integrals on cells
     for(index_t i = 0; i<workset.num_cells; i++)
     {
        h1norm += integral(i);
     }
  }
} 
  

template <typename EvalT,typename Traits>
void Norm_H1<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
   // Initialize accumulated error
   h1norm = Teuchos::ScalarTraits<ScalarT>::zero();
}

template <typename EvalT,typename Traits>
Teuchos::RCP<Teuchos::ParameterList> Norm_H1<EvalT,Traits>::getValidParameters() const
{
  
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Value", "string");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> integration_rule;
  p->set("IR", integration_rule);
  //Teuchos::RCP<charon::Names> n;
  p->set<std::string>("Name", "string");
  p->set<std::string>("Analytic Prefix", "string");
  p->set<std::string>("Error Prefix", "string");
  Teuchos::RCP<Teuchos::Comm<int> const> communicator;
  p->set("Comm", communicator);
  p->sublist("Norm ParameterList");
  p->sublist("Norm ParameterList").set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", Teuchos::rcp(new panzer::ParamLib));

  return p;
}

template <typename EvalT,typename Traits>
void Norm_H1<EvalT,Traits>::postEvaluate(typename Traits::PostEvalData /* d */)
{
  this->postprocess(std::cout);
  outputNorm->setValue(std::sqrt(globalH1Norm));
}


template<typename EvalT, typename Traits>
void Norm_H1<EvalT, Traits>::postprocess(std::ostream& os) 
{
  // throw unless specialized for residual evaluations
  os << "WEIRD! It is a bad idea to use this evaluation type for Norm_H1Error!" << std::endl;
}


template<>
void Norm_H1<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{
  // Sum the norm across processors and output it to the screen.
  globalH1Norm = 0.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &h1norm, &globalH1Norm);

  if (comm->getRank() == 0) {
     panzer::ios_all_saver saver(os);

     const std::string outputStr = "H1 Norm "+simulation_solution.fieldTag().name();
     std::size_t precision = 8;
     std::size_t name_width = outputStr.size();
     std::size_t value_width = precision + 7;

     os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
     os << std::setw(name_width) << outputStr
        << " " << std::setw(value_width) << std::sqrt(globalH1Norm)
        << std::endl;
  }
}

}


#endif
