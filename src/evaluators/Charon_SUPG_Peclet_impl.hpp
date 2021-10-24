
#ifndef CHARON_SUPG_PECLET_IMPL_HPP
#define CHARON_SUPG_PECLET_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SUPG_Peclet<EvalT, Traits>::
SUPG_Peclet(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& names =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  ir_degree = ir->cubature_degree;
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  ls_type = p.get<string>("Length Scale");
  carrType = p.get<string>("Carrier Type");

  bool includeSoret = false;
  if (p.isParameter("Include Soret Effect"))
    includeSoret = p.get<bool>("Include Soret Effect");

  if (carrType == "Electron")
  {
    peclet = MDField<ScalarT,Cell,Point>(names.field.elec_peclet,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.elec_diff_coeff,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.elec_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.elec_velocity,vector);
  }

  else if (carrType == "Hole")
  {
    peclet = MDField<ScalarT,Cell,Point>(names.field.hole_peclet,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.hole_diff_coeff,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.hole_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.hole_velocity,vector);
  }

  else if (carrType == "Ion")
  {
    peclet = MDField<ScalarT,Cell,Point>(names.field.ion_peclet,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.ion_diff_coeff,scalar);
    if (includeSoret)
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.ion_eff_velocity,vector);
    else
      velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.ion_velocity,vector);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole, or Ion !");

  this->addEvaluatedField(peclet);

  this->addDependentField(diffcoeff);
  this->addDependentField(velocity);

  std::string n = "SUPG_Carrier_Peclet";
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SUPG_Peclet<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  ir_index = panzer::getIntegrationRuleIndex(ir_degree,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SUPG_Peclet<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  typedef typename PHX::MDField<ScalarT,Cell,Point>::size_type size_type;

  using panzer::Cell;
  using panzer::IP;
  using panzer::Dim;
  using panzer::index_t;

  PHX::MDField<double,Cell,IP,Dim,Dim> gc =
    (workset.int_rules[ir_index])->contravarient;

  PHX::MDField<double,Cell,IP> norm_gc =
    (workset.int_rules[ir_index])->norm_contravarient;

  PHX::MDField<double,Cell,IP,Dim,Dim> gC =
    (workset.int_rules[ir_index])->covarient;

  std::vector<ScalarT> gc_vt(gc.dimension(2));  // gc.dimension(2) gives number of dimensions
  std::vector<ScalarT> gC_vt(gc.dimension(2));

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      for (size_type i = 0; i < gc.dimension(2); ++i)
      {
        gc_vt[i] = 0.0;
        gC_vt[i] = 0.0;
        for (size_type j = 0; j < gc.dimension(3); ++j)
        {
          const ScalarT vel = velocity(cell,point,j);
          gc_vt[i] += gc(cell,point,i,j) * vel;
          gC_vt[i] += gC(cell,point,i,j) * vel;
        }
      }

      ScalarT v_gc_vt = 0.0;  // contravarient
      ScalarT v_gC_vt = 0.0;  // covarient
      ScalarT veloc_mag2 = 0.0;
      for (size_type i = 0; i < gc.dimension(2); ++i)
      {
        const ScalarT vel = velocity(cell,point,i);
        v_gc_vt += vel * gc_vt[i];
        v_gC_vt += vel * gC_vt[i];
        veloc_mag2 += vel*vel;
      }

      // Calculate size of the element in the direction of velocity.
      ScalarT length_scale = 0.0;
      if (ls_type == "Stream")
      {
        if (v_gC_vt != 0.0)
          length_scale = 2.0 * std::sqrt(v_gC_vt/veloc_mag2);
      }

      if (ls_type == "Shakib")
      {
        ScalarT norm_gc2 = norm_gc(cell,point) * norm_gc(cell,point);
        if (v_gc_vt != 0.0)
          length_scale = 2.0 * std::sqrt(num_dims*v_gc_vt/(norm_gc2*veloc_mag2));
      }

      // Carrier peclet number, this expression is also used in Charon 1.0.
      // All quantities in the peclet number expression should be unitless (scaled quantities)
      ScalarT zero(0.0);
      peclet(cell,point) = (veloc_mag2 == 0.0 ? zero : std::sqrt(veloc_mag2) * length_scale / (2.0 * diffcoeff(cell,point) ) );
      //peclet(cell,point) = std::sqrt(veloc_mag2) * length_scale / (2.0 * diffcoeff(cell,point) );
    }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SUPG_Peclet<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Length Scale", "?");
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  p->set<bool>("Include Soret Effect", false);

  return p;
}

}

#endif
