
#ifndef CHARON_ANALYTIC_COMPARISON_IMPL_HPP
#define CHARON_ANALYTIC_COMPARISON_IMPL_HPP

#include <string>
#include <cmath>

#include "Kokkos_ViewFactory.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IosAllSaver.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace charon {

//**********************************************************************
// Straight difference of the analytic and computed solutions at the
// nodes.
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
AnalyticComparison<EvalT, Traits>::
AnalyticComparison(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;

  string name = params.get<string>("Name");

  Teuchos::RCP<PHX::DataLayout> data_layout = params.get< Teuchos::RCP<PHX::DataLayout> >("DataLayout");
  std::string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  std::string const errorPrefix = params.get< std::string >("Error Prefix");

  simulation_value = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  analytic_value = PHX::MDField<const ScalarT,Cell,Point>(analyticPrefix + name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);

  this->addDependentField(simulation_value);
  this->addDependentField(analytic_value);
  this->addEvaluatedField(error);

  string n = "Analytic Comparison: " + name;
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
AnalyticComparison<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;
  using panzer::index_t;


  // Simple difference of the simulation and analytic solutions at the
  // nodes.
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    for (size_type point = 0; point < error.dimension(1); ++point)
      error(cell,point) = simulation_value(cell,point) - analytic_value(cell,point);

}

//**********************************************************************
// A relative error at the nodes. The relative breaks down if the
// solution is at, or very close to, zero. In that case the error
// reverts to a straight difference of the analytic and computed
// solution.
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
AnalyticComparison_RelError<EvalT, Traits>::
AnalyticComparison_RelError(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;

  string name = params.get<string>("Name");

  useAbs = params.get<bool>("Use Absolute");

  Teuchos::RCP<PHX::DataLayout> data_layout = params.get< Teuchos::RCP<PHX::DataLayout> >("DataLayout");
  std::string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  std::string const errorPrefix = params.get< std::string >("Error Prefix");

  simulation_value = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  analytic_value = PHX::MDField<const ScalarT,Cell,Point>(analyticPrefix + name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);

  this->addDependentField(simulation_value);
  this->addDependentField(analytic_value);
  this->addEvaluatedField(error);

  string n = "Analytic Comparison: " + name;
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
AnalyticComparison_RelError<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;
  using panzer::index_t;


  // For efficiency it'd probably be better to brink the "useAbs"
  // outside the inner loops, but this makes for less code
  // duplication. Since it's just comparing values speed isn't really
  // critical anyway.
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (size_type point = 0; point < error.dimension(1); ++point)
    {
      ScalarT diffVals = simulation_value(cell,point) - analytic_value(cell,point);
      if (std::abs(diffVals) < std::numeric_limits<ScalarT>::min() ||
          std::abs(analytic_value(cell,point)) < std::numeric_limits<ScalarT>::min())
      {
        if (useAbs)
        {
          error(cell,point) = std::abs(diffVals);
        }
        else
        {
          error(cell,point) = diffVals;
        }
      }
      else
      {
        if (useAbs)
        {
          error(cell,point) = std::abs(diffVals / analytic_value(cell,point));
        }
        else
        {
          error(cell,point) = diffVals / std::abs(analytic_value(cell,point));
        }
      }

    }
  }

}

//**********************************************************************
// This computes an L2 error of the computed and analytic solution. This
// includes the numerical integral over the problem domain for a global
// measure of the overall error.
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
AnalyticComparison_L2Error<EvalT, Traits>::
AnalyticComparison_L2Error(
  const Teuchos::ParameterList&  params)
{
  using Teuchos::RCP;
  using std::string;
  using panzer::BASIS;
  using PHX::DataLayout;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  string name = params.get<string>("Name");
  comm = params.get<RCP<Teuchos::Comm<int> const> >("Comm");

  // Teuchos::RCP<PHX::DataLayout> data_layout = params.get< Teuchos::RCP<PHX::DataLayout> >("DataLayout");

  // Obtain BASIS information
  RCP<BasisIRLayout> basis = params.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  numNodes = data_layout->dimension(1);
  basisName = basis->name();

  string const analyticPrefix = params.get< std::string >("Analytic Prefix");
  string const errorPrefix = params.get< std::string >("Error Prefix");

  // Obtain IP information
  RCP<IntegrationRule const> ir = params.get<RCP<IntegrationRule> >("IR");
  quadOrder = ir->cubature_degree;
  numIPs = (ir->dl_scalar)->dimension(1);

  simulation_value = PHX::MDField<const ScalarT,Cell,Point>(name, data_layout);
  analytic_value = PHX::MDField<const ScalarT,Cell,Point>(analyticPrefix + name, data_layout);
  error = PHX::MDField<ScalarT,Cell,Point>(errorPrefix + name, data_layout);

  this->addDependentField(simulation_value);
  this->addDependentField(analytic_value);
  this->addEvaluatedField(error);

  string n = "Analytic Comparison: L2 Error" + name;
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
AnalyticComparison_L2Error<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  quadIndex = panzer::getIntegrationRuleIndex(quadOrder,(*sd.worksets_)[0]);
  basisIndex = panzer::getBasisIndex(basisName,(*sd.worksets_)[0]);

  integral = Kokkos::createDynRankView(error.get_static_view(), "integral", error.dimension(0));
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
AnalyticComparison_L2Error<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData  workset)
{
  using panzer::index_t;

  // FieldContainer to temporarily hold L2 error at IPs
  Kokkos::DynRankView<ScalarT,PHX::Device> error_ip = Kokkos::createDynRankView(error.get_static_view(), "error_ip", workset.num_cells, numIPs);

  // Zero-out arrays
  Kokkos::deep_copy(integral, ScalarT(0.0));
  Kokkos::deep_copy(error_ip, ScalarT(0.0));

  // Compute the error at BASIS points
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < numNodes; ++point)
    {
      error(cell,point) = simulation_value(cell,point) - analytic_value(cell,point);
      error(cell,point) *= error(cell,point);
    }
  }

  // Convert the error at BASIS points to IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int ip = 0; ip < numIPs; ++ip)
    {
      for (int node = 0; node < numNodes; ++node)
      {
        error_ip(cell,ip) += (workset.bases[basisIndex])->basis_scalar(cell,node,ip) * error(cell,node);
        // if (ip == 0) std::cout << "cell=" << cell << ", node=" << node << ", error=" << error(cell,node) << std::endl;
      }
      // std::cout << "cell=" << cell << ", ip=" << ip << ", error_ip=" << error_ip(cell,ip) << std::endl;
    }
  }

  // Numerically integrate the L2 error over the problem domain
  if (workset.num_cells > 0)
  {
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::
       integrate(integral, error_ip, (workset.int_rules[quadIndex])->weighted_measure.get_view());

     //Intrepid2::FunctionSpaceTools<PHX::exec_space>::
     //  integrate(integral, error.get_view(), (workset.int_rules[quadIndex])->weighted_measure.get_view());

     for(index_t i = 0; i < workset.num_cells; i++)
     {
        l2error += integral(i);
        // std::cout << "cell=" << i << ", integral=" << integral(i) << std::endl;
     }
     // std::cout << "l2error=" << l2error << std::endl;
  }

}

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT,typename Traits>
void AnalyticComparison_L2Error<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
  // Initialize accumulated error
  l2error = Teuchos::ScalarTraits<ScalarT>::zero();
}

///////////////////////////////////////////////////////////////////////////////
//
//  postEvaluate()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT,typename Traits>
void AnalyticComparison_L2Error<EvalT,Traits>::postEvaluate(typename Traits::PostEvalData /* d */)
{
  this->postprocess(std::cout);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postprocess()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void AnalyticComparison_L2Error<EvalT, Traits>::postprocess(std::ostream& os)
{
  // throw unless specialized for residual evaluations
  os << "WEIRD! It is a bad idea to use this evaluation type for AnalyticComparison_L2Error!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//
//  postprocess() (Residual specialization)
//
///////////////////////////////////////////////////////////////////////////////
template<>
void AnalyticComparison_L2Error<panzer::Traits::Residual, panzer::Traits>::postprocess(std::ostream& os)
{

  // Sum the error across processors and output it to the screen.
  ScalarT globalL2Error = 0.0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, static_cast<int>(1), &l2error, &globalL2Error);

  if (comm->getRank() == 0) {
     panzer::ios_all_saver saver(os);

     const std::string outputStr = "L2 Error "+simulation_value.fieldTag().name();
     std::size_t precision = 8;
     std::size_t name_width = outputStr.size();
     std::size_t value_width = precision + 7;

     os << std::scientific << std::showpoint << std::setprecision(precision) << std::left;
     os << std::setw(name_width) << outputStr
        << " " << std::setw(value_width) << std::sqrt(globalL2Error)
        << std::endl;
  }
}

}

#endif
