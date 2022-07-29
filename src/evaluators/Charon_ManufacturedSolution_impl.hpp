#ifndef CHARON_MANUFACTURED_SOLUTION_IMPL_HPP
#define CHARON_MANUFACTURED_SOLUTION_IMPL_HPP

#include <string>
#include <cmath>

#include "Kokkos_ViewFactory.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IosAllSaver.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace charon {



template <typename EvalT,typename Traits>
AnalyticSolution<EvalT,Traits>::AnalyticSolution(const std::string & name,  const Teuchos::RCP<panzer::IntegrationRule const> & ir) 
{
using Teuchos::RCP;
using panzer::Cell;
using panzer::Point;
using panzer::Dim;

Teuchos::RCP<PHX::DataLayout> data_layout_scalar = ir->dl_scalar;
Teuchos::RCP<PHX::DataLayout> data_layout_vector = ir->dl_vector;
ir_degree = ir->cubature_degree;

analytic_solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout_scalar);
analytic_solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_"+name, data_layout_vector);

this->addEvaluatedField(analytic_solution);
this->addEvaluatedField(analytic_solution_grad);
  
std::string n = "Analytic Solution";
this->setName(n);

}


template <typename EvalT,typename Traits>
void AnalyticSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd, PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0]);
}


template <typename EvalT,typename Traits>
void AnalyticSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim; 
  using panzer::index_t;
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < analytic_solution.extent_int(1); ++point) {

      const double & x = workset.int_rules[ir_index]->ip_coordinates(cell,point,0);
      const double & y = workset.int_rules[ir_index]->ip_coordinates(cell,point,1);

      analytic_solution(cell,point) = std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
      analytic_solution_grad(cell,point,0) = 2.0*M_PI*std::cos(2*M_PI*x)*std::sin(2*M_PI*y);
      analytic_solution_grad(cell,point,1) = 2.0*M_PI*std::sin(2*M_PI*x)*std::cos(2*M_PI*y);
    }   
  }
}


template <typename EvalT,typename Traits>
DD_RDH_1_AnalyticSolution<EvalT,Traits>::DD_RDH_1_AnalyticSolution(const Teuchos::RCP<panzer::IntegrationRule const> & ir) 
{
using Teuchos::RCP;
using panzer::Cell;
using panzer::Point;
using panzer::Dim;

Teuchos::RCP<PHX::DataLayout> data_layout_scalar = ir->dl_scalar;
Teuchos::RCP<PHX::DataLayout> data_layout_vector = ir->dl_vector;
ir_degree = ir->cubature_degree;
//std::cout<<"******NAME:"<<name<<std::endl;
//analytic_solution = PHX::MDField<ScalarT,Cell,Point>(name, data_layout_scalar);
//analytic_solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_"+name, data_layout_vector);
analytic_phi_solution = PHX::MDField<ScalarT,Cell,Point>("analytic_ELECTRIC_POTENTIAL", data_layout_scalar);
analytic_phi_solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_analytic_ELECTRIC_POTENTIAL", data_layout_vector);
analytic_edensity_solution = PHX::MDField<ScalarT,Cell,Point>("analytic_ELECTRON_DENSITY", data_layout_scalar);
analytic_edensity_solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_analytic_ELECTRON_DENSITY", data_layout_vector);
analytic_hdensity_solution = PHX::MDField<ScalarT,Cell,Point>("analytic_HOLE_DENSITY", data_layout_scalar);
analytic_hdensity_solution_grad = PHX::MDField<ScalarT,Cell,Point,Dim>("GRAD_analytic_HOLE_DENSITY", data_layout_vector);

//this->addEvaluatedField(analytic_solution);
//this->addEvaluatedField(analytic_solution_grad);
this->addEvaluatedField(analytic_phi_solution);
this->addEvaluatedField(analytic_phi_solution_grad);
this->addEvaluatedField(analytic_edensity_solution);
this->addEvaluatedField(analytic_edensity_solution_grad);
this->addEvaluatedField(analytic_hdensity_solution);
this->addEvaluatedField(analytic_hdensity_solution_grad);
  
std::string n = "Analytic Solution";
this->setName(n);

}


template <typename EvalT,typename Traits>
void DD_RDH_1_AnalyticSolution<EvalT,Traits>::postRegistrationSetup(typename Traits::SetupData sd, PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0]);
}


template <typename EvalT,typename Traits>
void DD_RDH_1_AnalyticSolution<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim; 
  using panzer::index_t;

  using std::atan;
  double vl = 0.4;
  double beta = 1e6;
  double length = 5e-5;
  double V0 = 2.5852029e-2; //From printout of Gary's V0
  double C0 = 1.0e10;
  double x_scale = 1.0e-4;
  double nhat = 3e3;
  double phat = 2e3;
  //-vl * (atan(beta*(x - length*0.5))/atan(beta*length*0.5));


  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < analytic_phi_solution.extent_int(1); ++point) {
    //for (int point = 0; point < analytic_phi_solution.extent_int(1); ++point) {

      const double & x = workset.int_rules[ir_index]->ip_coordinates(cell,point,0);
      //const double & y = workset.int_rules[ir_index]->ip_coordinates(cell,point,1);
      
      //analytic_solution(cell,point) = -vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0;
      //analytic_solution_grad(cell,point,0) = -vl*x_scale * ((beta/(1+(beta*(x*x_scale - length*0.5)*beta*(x*x_scale - length*0.5))))/atan(beta*length*0.5))/ V0;
      //analytic_solution_grad(cell,point,1) = 0.0;

      analytic_phi_solution(cell,point) = -vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0;
      analytic_phi_solution_grad(cell,point,0) = -vl*x_scale * ((beta/(1+(beta*(x*x_scale - length*0.5)*beta*(x*x_scale - length*0.5))))/atan(beta*length*0.5))/ V0;
      analytic_phi_solution_grad(cell,point,1) = 0.0;

      analytic_edensity_solution(cell,point) = nhat*std::exp(-vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0)/C0;
      analytic_edensity_solution_grad(cell,point,0) = nhat*std::exp(-vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0)*-vl*x_scale * ((beta/(1+(beta*(x*x_scale - length*0.5)*beta*(x*x_scale - length*0.5))))/atan(beta*length*0.5))/ V0/C0;
      analytic_edensity_solution_grad(cell,point,1) = 0.0;

      analytic_hdensity_solution(cell,point) = phat*std::exp(vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0)/C0;
      analytic_hdensity_solution_grad(cell,point,0) = -phat*std::exp(vl * (atan(beta*(x*x_scale - length*0.5))/atan(beta*length*0.5))/ V0)*-vl*x_scale * ((beta/(1+(beta*(x*x_scale - length*0.5)*beta*(x*x_scale - length*0.5))))/atan(beta*length*0.5))/ V0/C0;
      analytic_hdensity_solution_grad(cell,point,1) = 0.0;



    }   
  }
}


}
#endif
