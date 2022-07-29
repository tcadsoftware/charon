#ifndef CHARON_NORM_CALCULATION_HPP
#define CHARON_NORM_CALCULATION_HPP

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Panzer_Dimension.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Panzer_ScalarParameterEntry.hpp" //For norm output to exo
#include "Panzer_ParameterLibrary.hpp" //For norm output to exo

namespace charon {
 
  using panzer::Cell;
  using panzer::Point;
  using panzer::Dim;



template<typename EvalT, typename Traits>
class Norm_L2Error
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Norm_L2Error(
      const Teuchos::ParameterList& p); 

    Norm_L2Error(){};

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d); 

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  void postprocess(std::ostream& os);
  void preEvaluate(typename Traits::PreEvalData d); 
  void postEvaluate(typename Traits::PostEvalData d); 
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
  PHX::MDField<const ScalarT,Cell,Point> simulation_solution;
  PHX::MDField<ScalarT,Cell,Point> error;
  PHX::MDField<const ScalarT,Cell,Point> manufactured_solution;
  //PHX::MDField<ScalarT,Cell,Point,Dim> manufactured_gradient;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > outputNorm; // For exo output
  

  ScalarT l2error;
  Kokkos::DynRankView<ScalarT,PHX::Device> integral;


  int num_dims;
  int quadOrder, quadIndex;
  int numIPs, numNodes, basisIndex;
  std::string basisName;
  std::string outputNormName;
  ScalarT globalL2Error;


  Teuchos::RCP<Teuchos::Comm<int> const> comm;

}; // end of class Norm_L2Error


template<typename EvalT, typename Traits>
class Norm_H1Error
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Norm_H1Error(
      const Teuchos::ParameterList& p); 

    Norm_H1Error(){};

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d); 

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  void postprocess(std::ostream& os);
  void preEvaluate(typename Traits::PreEvalData d); 
  void postEvaluate(typename Traits::PostEvalData d); 
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
  PHX::MDField<const ScalarT,Cell,Point> simulation_solution;
  PHX::MDField<ScalarT,Cell,Point> error;
  PHX::MDField<const ScalarT,Cell,Point> manufactured_solution;
  PHX::MDField<const ScalarT,Cell,Point,Dim> manufactured_gradient;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > outputNorm; // For exo output
  

  ScalarT h1error;
  Kokkos::DynRankView<ScalarT,PHX::Device> integral;


  int num_dims;
  int quadOrder, quadIndex;
  int numIPs, numNodes, basisIndex;
  std::string basisName;
  std::string outputNormName;
  ScalarT globalH1Error;


  Teuchos::RCP<Teuchos::Comm<int> const> comm;

}; // end of class Norm_H1Error


template<typename EvalT, typename Traits>
class Norm_L2
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Norm_L2(
      const Teuchos::ParameterList& p); 

    Norm_L2(){};

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d); 

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  void postprocess(std::ostream& os);
  void preEvaluate(typename Traits::PreEvalData d); 
  void postEvaluate(typename Traits::PostEvalData d); 
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
  PHX::MDField<const ScalarT,Cell,Point> simulation_solution;
  PHX::MDField<ScalarT,Cell,Point> error;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > outputNorm; // For exo output
  

  ScalarT l2norm;
  Kokkos::DynRankView<ScalarT,PHX::Device> integral;


  int num_dims;
  int quadOrder, quadIndex;
  int numIPs, numNodes, basisIndex;
  std::string basisName;
  std::string outputNormName;
  ScalarT globalL2Norm;


  Teuchos::RCP<Teuchos::Comm<int> const> comm;

}; // end of class Norm_L2


template<typename EvalT, typename Traits>
class Norm_H1
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    Norm_H1(
      const Teuchos::ParameterList& p); 

    Norm_H1(){};

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d); 

  private:

    using ScalarT = typename EvalT::ScalarT;

private:
  void postprocess(std::ostream& os);
  void preEvaluate(typename Traits::PreEvalData d); 
  void postEvaluate(typename Traits::PostEvalData d); 
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;
  
  PHX::MDField<const ScalarT,Cell,Point> simulation_solution;
  PHX::MDField<ScalarT,Cell,Point,Dim> manufactured_gradient;
  PHX::MDField<ScalarT,Cell,Point> error;
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > outputNorm; // For exo output
  

  ScalarT h1norm;
  Kokkos::DynRankView<ScalarT,PHX::Device> integral;


  int num_dims;
  int quadOrder, quadIndex;
  int numIPs, numNodes, basisIndex;
  std::string basisName;
  std::string outputNormName;
  ScalarT globalH1Norm;


  Teuchos::RCP<Teuchos::Comm<int> const> comm;

}; // end of class Norm_H1


}

#endif
