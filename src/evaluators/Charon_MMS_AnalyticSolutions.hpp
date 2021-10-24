
#ifndef CHARON_MMS_ANALYTICSOLUTIONS_HPP
#define CHARON_MMS_ANALYTICSOLUTIONS_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_GlobalData.hpp"

#include "Charon_config.hpp"
#include "Charon_Names.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Charon_Physical_Constants.hpp"
#include "Charon_MMS_AnalyticFunctions.hpp"

/**
 * This file contains the analytic solutions of MMS problems used in
 * Charon.
 */

namespace charon {

template<typename EvalT, typename Traits>
class MMS_AnalyticSolutionAlg :
    public PHX::EvaluatorWithBaseImpl<Traits>,
    public PHX::EvaluatorDerived<EvalT, Traits>  {

public:

  virtual ~MMS_AnalyticSolutionAlg() = 0;


protected:
  MMS_AnalyticSolutionAlg() {;}

  typedef typename EvalT::ScalarT ScalarT;

  template <typename T> struct Type {}; // used for member specialization


  PHX::MDField<const ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> getCoordField(Type<panzer::BASIS>,
                                                                 std::string const& fieldName,
                                                                 panzer::FieldLibraryBase const& fieldLayoutLibrary,
                                                                 Teuchos::RCP<panzer::IntegrationRule const> const& ir);

  Teuchos::RCP<PHX::DataLayout> getFieldDataLayout(Type<panzer::BASIS>,
                                                   const std::string & fieldName,
                                                   const panzer::FieldLibraryBase & fieldLayoutLibrary,
                                                   const Teuchos::RCP<const panzer::IntegrationRule> & /* ir */)
    { return fieldLayoutLibrary.lookupBasis(fieldName)->functional; }

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters();

  // Required fields
  PHX::MDField<const ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> coordinates;
  double V0;

};


/**
 * \brief Class for analytic solution of an NLP MMS.
 *
 * This class is used to compute the analytic solution to a manufactured
 * solution for the nonlinear Poisson equation. The solution is given by
 * \f$V(x) = 0.3* \text{erfc}\left[ 20,000 * (x - 0.00025)\right] -
 * 0.3\f$. There is no source term associated with this MMS, instead
 * there is an analytic function for the doping that serves the same
 * purpose. See the included documentation for a description of that.
 */
template<typename EvalT, typename Traits>
class MMS_NLP_GLH_1_AnalyticSolution :
  public MMS_AnalyticSolutionAlg<EvalT, Traits> {

using MMS_AnalyticSolutionAlg<EvalT, Traits>::coordinates;
using MMS_AnalyticSolutionAlg<EvalT, Traits>::V0;

public:


  MMS_NLP_GLH_1_AnalyticSolution(std::string const& prefix,
                                 charon::Names const& names,
                                 Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                                 Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                                 Teuchos::ParameterList const& options);

  void evaluateFields(typename Traits::EvalData d);

private:


  typedef typename EvalT::ScalarT ScalarT;

  template <typename T> struct TypeD : public MMS_AnalyticSolutionAlg<EvalT, Traits>::template Type<T> {};

  // Simulation solution
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> analytic_phi;

};

/**
 * @brief The analytic solution for a semiconductor drift-diffusion MMS.
 *
 * This class contains the analytic solution for a 1D, steady-state
 * drift-diffusion MMS. The MMS is complete in that it represents all
 * three drift-diffusion equations with the appropriate coupling between
 * equations. The solutions for the unknowns are:
 *
 * - @f$\psi\left(x\right) = \text{VH} \left[\text{erfc}\left(
 *    \left(x - \text{VXC}\right) \right) - 1 \right]@f$
 * - @f$n\left(x\right) = 10^{\left(\text{NDMAX}\
 *    \text{erfc}\left(\text{NW}\left(x - \text{NXC}\right)\right) +
 *    \text{NSHIFT}\right)}@f$
 * - @f$p\left(x\right) = 10^{\left(\text{NAMAX}\
 *    \text{erfc}\left(\text{PW}\left(\text{PXC} - x\right)\right) +
 *    \text{PSHIFT}\right)}@f$
 *
 * The parameters \em VH, \em VXC, \em NDMAX, \em NW, \em NXC, \em
 * NSHIFT, \em NAMAX, \em PW, \em PXC and \em PSHIFT are specified via
 * \em DEBUG \em PARAMETER statements in the input file.
 *
 * The problem should be run on a rectangle 0.5um in length in the x
 * direction. The y direction is not significant but is necessary for
 * %Charon. The material parameters must currently be constant and
 * should be set to:
 *
 *   - relative permittivity - 11.8
 *   - electron mobility - 1000.0
 *   - electron diffusion coefficient - 25.85215
 *   - hole mobility - 500.0
 *   - hole diffusion coefficient - 12.92608
 */
template<typename EvalT, typename Traits>
class MMS_DD_RDH_AnalyticSolution :
  public MMS_AnalyticSolutionAlg<EvalT, Traits> {

//using MMS_AnalyticSolutionAlg<EvalT, Traits>::coordinates;
//using MMS_AnalyticSolutionAlg<EvalT, Traits>::V0;

public:

  virtual ~MMS_DD_RDH_AnalyticSolution() = 0;

  void evaluateFields(typename Traits::EvalData d);


protected:

  MMS_DD_RDH_AnalyticSolution() {;}

  Teuchos::RCP<charon::MMS_DD_RDH_AnalyticFunction> mms;

  typedef typename EvalT::ScalarT ScalarT;

  template <typename T> struct TypeD : public MMS_AnalyticSolutionAlg<EvalT, Traits>::template Type<T> {};

  // Required fields
  PHX::MDField<const ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim> coordinates;

  // scaling parameters
  double C0;
  double V0;

  // Simulation solution
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> analytic_phi;
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> analytic_edensity;
  PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS> analytic_hdensity;

};

/**
 * @brief The analytic solution for a semiconductor drift-diffusion MMS.
 *
 * This class contains the mobility independent version
 *
 * */

template<typename EvalT, typename Traits>
class MMS_DD_RDH_1_AnalyticSolution :
  public MMS_DD_RDH_AnalyticSolution<EvalT, Traits> {

using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::coordinates;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::V0;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::C0;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_phi;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_edensity;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_hdensity;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::mms;

public:

  MMS_DD_RDH_1_AnalyticSolution(std::string const& prefix,
                                 charon::Names const& names,
                                 Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                                 Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                                 Teuchos::ParameterList const& options);

private:

  typedef typename EvalT::ScalarT ScalarT;

  template <typename T> struct TypeD : public MMS_AnalyticSolutionAlg<EvalT, Traits>::template Type<T> {};
};

/**
 * @brief The analytic solution for a semiconductor drift-diffusion MMS.
 *
 * This class contains the mobility dependent version
 *
 * */

template<typename EvalT, typename Traits>
class MMS_DD_RDH_2_AnalyticSolution :
  public MMS_DD_RDH_AnalyticSolution<EvalT, Traits> {

using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::coordinates;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::V0;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::C0;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_phi;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_edensity;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::analytic_hdensity;
using MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::mms;

public:

  MMS_DD_RDH_2_AnalyticSolution(std::string const& prefix,
                                 charon::Names const& names,
                                 Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                                 Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                                 Teuchos::ParameterList const& options);

private:

  typedef typename EvalT::ScalarT ScalarT;

  template <typename T> struct TypeD : public MMS_AnalyticSolutionAlg<EvalT, Traits>::template Type<T> {};
};
}

#endif
