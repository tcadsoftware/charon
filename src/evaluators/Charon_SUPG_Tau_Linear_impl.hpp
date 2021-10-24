
#ifndef CHARON_SUPG_TAU_LINEAR_IMPL_HPP
#define CHARON_SUPG_TAU_LINEAR_IMPL_HPP

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
SUPG_Tau_Linear<EvalT, Traits>::
SUPG_Tau_Linear(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  ir_degree = ir->cubature_degree;
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  add_source_stab = p.get<bool>("Add Source Stabilization");
  carrType = p.get<string>("Carrier Type");

  if (carrType == "Electron")
  {
    tau = MDField<ScalarT,Cell,Point>(names.field.tau_stab_e,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.elec_diff_coeff,scalar);
    velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.elec_velocity,vector);
  }
  else if (carrType == "Hole")
  {
    tau = MDField<ScalarT,Cell,Point>(names.field.tau_stab_h,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.hole_diff_coeff,scalar);
    velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.hole_velocity,vector);
  }
  else if (carrType == "Ion")
  {
    tau = MDField<ScalarT,Cell,Point>(names.field.tau_stab_ion,scalar);
    diffcoeff = MDField<const ScalarT,Cell,Point>(names.field.ion_diff_coeff,scalar);
    velocity = MDField<const ScalarT,Cell,Point,Dim>(names.field.ion_velocity,vector);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron, Hole, or Ion !");

  this->addEvaluatedField(tau);

  this->addDependentField(diffcoeff);
  this->addDependentField(velocity);

  if (add_source_stab)
  {
    if (carrType == "Electron")
      recomb_deriv = MDField<const ScalarT,Cell,Point>(names.field.recomb_deriv_e,scalar);
    else if (carrType == "Hole")
      recomb_deriv = MDField<const ScalarT,Cell,Point>(names.field.recomb_deriv_h,scalar);

    this->addDependentField(recomb_deriv);
  }

  std::string n = "SUPG_Tau_Linear";
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SUPG_Tau_Linear<EvalT, Traits>::
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
SUPG_Tau_Linear<EvalT, Traits>::
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

  std::vector<ScalarT> gc_vt(gc.dimension(2));  // gc.dimension(2) gives number of dimensions

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      for (size_type i = 0; i < gc.dimension(2); ++i)
      {
        gc_vt[i] = 0.0;
        for (size_type j = 0; j < gc.dimension(3); ++j)
          gc_vt[i] += gc(cell,point,i,j) * velocity(cell,point,j);
      }

      ScalarT v_gc_vt = 0.0;
      for (size_type i = 0; i < gc.dimension(2); ++i)
        v_gc_vt += velocity(cell,point,i) * gc_vt[i];

      ScalarT tmp = 0.0;

      // convection stab.
      if (std::fabs(v_gc_vt) > 0.)
        tmp += sqrt(v_gc_vt);

      // diffusion stab.
      tmp += sqrt(3.0) * diffcoeff(cell,point) * norm_gc(cell,point);

      // source stab.
      if (add_source_stab)
        tmp += std::fabs( recomb_deriv(cell,point) );

      // transient stab.
      if (workset.evaluate_transient_terms)
      {
        // TODO:
      }

      tau(cell,point) = 1.0 / tmp;

    }  // end of loop over QPs

  }  // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
SUPG_Tau_Linear<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set<std::string>("Carrier Type", "?");
  p->set<bool>("Add Source Stabilization", false);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  return p;
}

}

#endif
