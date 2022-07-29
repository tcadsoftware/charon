
#ifndef CHARON_MMS_ANALYTICSOLUTIONS_IMPL_HPP
#define CHARON_MMS_ANALYTICSOLUTIONS_IMPL_HPP

#include <cmath>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_GatherIntegrationCoordinates.hpp"

#include "Charon_Names.hpp"
#include "Charon_Scaling_Parameters.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Destructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
MMS_AnalyticSolutionAlg<EvalT, Traits>::~MMS_AnalyticSolutionAlg()
{
}


///////////////////////////////////////////////////////////////////////////////
//
//  getCoordField()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
PHX::MDField<const typename MMS_AnalyticSolutionAlg<EvalT, Traits>::ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim>
MMS_AnalyticSolutionAlg<EvalT, Traits>::getCoordField(
  Type<panzer::BASIS>,
  std::string const& fieldName,
  panzer::FieldLibraryBase const& fieldLayoutLibrary,
  Teuchos::RCP<panzer::IntegrationRule const> const& /* ir */)
{
  Teuchos::RCP<panzer::PureBasis const> basis = fieldLayoutLibrary.lookupBasis(fieldName);
  std::string const coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis->name());
  return PHX::MDField<const ScalarT,panzer::Cell,panzer::BASIS,panzer::Dim>(coordName, basis->coordinates);

}

//**********************************************************************
// MMS_NLP_GLH_1
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
MMS_NLP_GLH_1_AnalyticSolution<EvalT, Traits>::
MMS_NLP_GLH_1_AnalyticSolution(std::string const& prefix,
                               charon::Names const& names,
                               Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                               Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                               Teuchos::ParameterList const& options)
{

  Teuchos::RCP<PHX::DataLayout> data_layout = this->getFieldDataLayout(TypeD<panzer::BASIS>(),
                                                                 names.dof.phi,
                                                                 *fieldLayoutLibrary, ir);
  analytic_phi = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.phi, data_layout);
  this->addEvaluatedField(analytic_phi);

  // Coordinates are required for the analytic solution and scale factor
  // for the electric potential
  coordinates = this->getCoordField(TypeD<panzer::BASIS>(),names.dof.phi,*fieldLayoutLibrary,ir);
  this->addDependentField(coordinates);

  // scaling parameter
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = options.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;

  // Name for this evaluator
  this->setName("MMS_NLP_GLH_1_AnalyticSolution");
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
void MMS_NLP_GLH_1_AnalyticSolution<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;
  typedef typename PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>::size_type size_type;

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (size_type point = 0; point < analytic_phi.dimension(1); ++point) {

      double const x = Sacado::ScalarValue<ScalarT>::eval(coordinates(cell,point,0));

      // This assumes the x coordinate is given in microns. The function
      // itself was originally formulated in CGS units thus the
      // multiplication by 1.0e-4.
      analytic_phi(cell, point) = (0.3 * std::erfc(20000*(x*1.0e-4 - 2.5e-4)) - 0.3) / V0;
    }
  }
}

//**********************************************************************
// MMS_DD_RDH Base
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Destructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::~MMS_DD_RDH_AnalyticSolution()
{
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
void MMS_DD_RDH_AnalyticSolution<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  using std::atan;
  using std::pow;
  using panzer::index_t;

  typedef typename PHX::MDField<ScalarT,panzer::Cell,panzer::BASIS>::size_type size_type;

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (size_type point = 0; point < analytic_phi.dimension(1); ++point) {

      // This assumes the x coordinate is given in microns. The function
      // itself was originally formulated in CGS units thus the
      // multiplication by 1.0e-4.
      double const x = Sacado::ScalarValue<ScalarT>::eval(coordinates(cell,point,0))*1.0e-4;

      analytic_phi(cell,point) = mms->potential(x) / V0;
      analytic_edensity(cell,point) = mms->edensity(x) / C0;
      analytic_hdensity(cell,point) = mms->hdensity(x) / C0;
    }
  }
}

//**********************************************************************
// MMS_DD_RDH_1
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
MMS_DD_RDH_1_AnalyticSolution<EvalT, Traits>::
MMS_DD_RDH_1_AnalyticSolution(std::string const& prefix,
                               charon::Names const& names,
                               Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                               Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                               Teuchos::ParameterList const& options)
{

  Teuchos::RCP<PHX::DataLayout> data_layout = this->getFieldDataLayout(TypeD<panzer::BASIS>(),
                                                                 names.dof.phi,
                                                                 *fieldLayoutLibrary, ir);

  analytic_phi = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.phi, data_layout);
  analytic_edensity = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.edensity, data_layout);
  analytic_hdensity = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.hdensity, data_layout);
  this->addEvaluatedField(analytic_phi);
  this->addEvaluatedField(analytic_edensity);
  this->addEvaluatedField(analytic_hdensity);

  // Coordinates are required for the analytic solution and scale factor
  // for the electric potential
  coordinates = this->getCoordField(TypeD<panzer::BASIS>(),names.dof.phi,*fieldLayoutLibrary,ir);
  this->addDependentField(coordinates);

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = options.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  C0 = scaleParams->scale_params.C0;
  std::cout << std::scientific << std::showpoint << std::setprecision(12) << std::left;

  // Name for this evaluator
  this->setName("MMS_DD_RDH_1_AnalyticSolution");

  mms = Teuchos::rcp(new charon::MMS_DD_RDH_1_AnalyticFunction(options));

}


//**********************************************************************
// MMS_DD_RDH_2
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits>
MMS_DD_RDH_2_AnalyticSolution<EvalT, Traits>::
MMS_DD_RDH_2_AnalyticSolution(std::string const& prefix,
                               charon::Names const& names,
                               Teuchos::RCP<panzer::FieldLibraryBase const> const& fieldLayoutLibrary,
                               Teuchos::RCP<panzer::IntegrationRule const> const& ir,
                               Teuchos::ParameterList const& options)
{

  Teuchos::RCP<PHX::DataLayout> data_layout = this->getFieldDataLayout(TypeD<panzer::BASIS>(),
                                                                 names.dof.phi,
                                                                 *fieldLayoutLibrary, ir);

  analytic_phi = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.phi, data_layout);
  analytic_edensity = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.edensity, data_layout);
  analytic_hdensity = PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS>(prefix+names.dof.hdensity, data_layout);
  this->addEvaluatedField(analytic_phi);
  this->addEvaluatedField(analytic_edensity);
  this->addEvaluatedField(analytic_hdensity);

  // Coordinates are required for the analytic solution and scale factor
  // for the electric potential
  coordinates = this->getCoordField(TypeD<panzer::BASIS>(),names.dof.phi,*fieldLayoutLibrary,ir);
  this->addDependentField(coordinates);

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = options.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  C0 = scaleParams->scale_params.C0;

  // Name for this evaluator
  this->setName("MMS_DD_RDH_2_AnalyticSolution");

  mms = Teuchos::rcp(new charon::MMS_DD_RDH_2_AnalyticFunction(options));

}

}

#endif
