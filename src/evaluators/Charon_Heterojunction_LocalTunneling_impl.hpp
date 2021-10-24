
#ifndef CHARON_HETEROJUNCTION_LOCALTUNNELING_IMPL_HPP
#define CHARON_HETEROJUNCTION_LOCALTUNNELING_IMPL_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

namespace charon {

/**
 * @brief This evaluator computes the local tunneling (LT) factor (1+tunnel), which
 * multiplies the net thermionic emission (TE) current density to produce the
 * normal TE+LT flux boundary condition at a heterojunction (HJ).
 * The LT model is implemented for both FEM-SUPG and CVFEM-SG discretization
 * schemes. For FEM-SUGP, the LT factor is computed at the FEM integration points
 * (IPs) of a HJ; while for CVFEM-SG, the LT factor is computed at the CVFEM IPs
 * of a HJ.
*/

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Heterojunction_LocalTunneling<EvalT, Traits>::
Heterojunction_LocalTunneling(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve the data layout
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  // Retrieve parameters
  flux_name = p.get<string>("Flux Name");
  bandOffset = p.get<double>("Band Offset");
  tunnelMass = p.get<double>("Tunneling Effective Mass");
  detailIndex = p.get<int>("Details Index");

  if (p.get<string>("Primary Side") == "Left")
    isPrimarySideLeft = true;
  else
    isPrimarySideLeft = false;

  // Evaluated field
  flux_local_tunnel = MDField<ScalarT,Cell,Point>(flux_name, scalar);
  this->addEvaluatedField(flux_local_tunnel);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  E0 = scaleParams->scale_params.E0;

  // Dependent field
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);
  normal_dot_grad = MDField<const ScalarT,Cell,Point>(p.get<string>("Normal Dot Gradient"), scalar);

  this->addDependentField(latt_temp);
  this->addDependentField(normal_dot_grad);

  std::string name = "Heterojunction_LocalTunneling";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Heterojunction_LocalTunneling<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Heterojunction_LocalTunneling<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double h = cpc.h;       // [J.s]
  double q = cpc.q;       // [C]
  double m0 = cpc.m0;     // [kg]
  double pi = cpc.pi;
  double kb = cpc.kb*q;   // [J/K]

  // 100.0 is used to convert from [meter] to [centimeter]
  double prefactor = std::pow(3.0*h*q/(8.0*pi*std::sqrt(2.0*m0*tunnelMass))*100.0, 2.0/3.0);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain lattice temperature in [K]
      const ScalarT& latT = latt_temp(cell,point) * T0;
      ScalarT kbT = kb * latT;  // [J]

      ScalarT intUpLimit = std::fabs(bandOffset) / (kbT/q);   // [unitless], note bandOffset can be negative

      ScalarT eField = std::fabs(normal_dot_grad(cell,point)) * E0;   // [V/cm]

      ScalarT fieldParam = 1.0/kbT * prefactor * std::pow(eField, 2.0/3.0);

      ScalarT integralValue = evaluateIntegration(intUpLimit, fieldParam);

      flux_local_tunnel(cell,point) = 1.0 + integralValue;

      // std::cout << "intUpLimit=" << intUpLimit << ", eField=" << eField <<", fieldParam="
      // std::cout << fieldParam << ", flux_local_tunnel=" << flux_local_tunnel(cell,point) << std::endl;

    }  // end of loop over points
  }  // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Heterojunction_LocalTunneling<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set<std::string>("Flux Name", "");
  p->set<std::string>("Primary Side", "", "Determine if Primary Side is Left or Right");
  p->set<std::string>("Normal Dot Gradient", "");

  p->set<double>("Band Offset", 0., "Either conduction or valence band offset in units of [eV], can be a positive or negative value");
  p->set<double>("Tunneling Effective Mass", 0., "Effective mass used for the local tunneling calculation");

  p->set<int>("Details Index", 0, "Determine if on side 1 or on side 2");

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}


// Perform the integration: int_{0}^{intUpLimit} exp(u-(u/fieldParam)^1.5) du

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateIntegration()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
typename Heterojunction_LocalTunneling<EvalT,Traits>::ScalarT
Heterojunction_LocalTunneling<EvalT, Traits>::evaluateIntegration
  (const ScalarT & intUpLimit, const ScalarT & fieldParam)
{
  int num_points = 1000;

  ScalarT du = intUpLimit / double(num_points);
  ScalarT tunnel = 0.0;

  // perform the summation
  for (int i = 0; i < num_points; i++)
  {
    ScalarT u = (i + 0.5) * du;
    tunnel += std::exp(u - std::pow(u/fieldParam, 1.5));
  }

  // multiply the constant step size
  tunnel *= du;

  return tunnel;
}

}

#endif
