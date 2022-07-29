#ifndef CHARON_MANUFACTURED_SOLUTION_HPP
#define CHARON_MANUFACTURED_SOLUTION_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Dimension.hpp"

#include "Kokkos_DynRankView.hpp"


namespace charon {

  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;

template<typename EvalT, typename Traits>
class AnalyticSolution : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

    //using panzer::Cell;
    //using panzer::Point;
    //using panzer::Dim;

    public:
        AnalyticSolution(const std::string & name, const Teuchos::RCP<panzer::IntegrationRule const> & ir);
        AnalyticSolution(){};
        void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& fm);

        void evaluateFields(typename Traits::EvalData d); 

        //typedef typename EvalT::ScalarT ScalarT;
        //PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_solution;
        //PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_solution_grad;

    private:
        typedef typename EvalT::ScalarT ScalarT;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_solution;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_solution_grad;
        int ir_degree, ir_index;
};


//template<typename EvalT, typename Traits>
//class DD_RDH_1_AnalyticSolution : public PHX::EvaluatorWithBaseImpl<Traits>,
//                        public PHX::EvaluatorDerived<EvalT, Traits>  {

//    public:
//        DD_RDH_1_AnalyticSolution(const std::string & name, const Teuchos::RCP<panzer::IntegrationRule const> & ir);
//        DD_RDH_1_AnalyticSolution(){};
//        void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& fm);

//        void evaluateFields(typename Traits::EvalData d); 


//    private:
//        typedef typename EvalT::ScalarT ScalarT;
//        PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_solution;
//        PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_solution_grad;
//        int ir_degree, ir_index;
//};


template<typename EvalT, typename Traits>
class DD_RDH_1_AnalyticSolution : public PHX::EvaluatorWithBaseImpl<Traits>,
                        public PHX::EvaluatorDerived<EvalT, Traits>  {

    public:
        //AnalyticSolution(const std::string & name, const Teuchos::RCP<panzer::IntegrationRule const> & ir);
        //AnalyticSolution(){};
        DD_RDH_1_AnalyticSolution(const Teuchos::RCP<panzer::IntegrationRule const> & ir);
        DD_RDH_1_AnalyticSolution(){};
        void postRegistrationSetup(typename Traits::SetupData d, PHX::FieldManager<Traits>& fm);

        void evaluateFields(typename Traits::EvalData d); 


    private:
        typedef typename EvalT::ScalarT ScalarT;
        //PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_solution;
        //PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_solution_grad;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_phi_solution;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_phi_solution_grad;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_edensity_solution;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_edensity_solution_grad;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point> analytic_hdensity_solution;
        PHX::MDField<ScalarT,panzer::Cell,panzer::Point,panzer::Dim> analytic_hdensity_solution_grad;
        int ir_degree, ir_index;
};
}
#endif
