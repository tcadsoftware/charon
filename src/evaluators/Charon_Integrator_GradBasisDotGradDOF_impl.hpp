
#ifndef CHARON_GRADBASISDOTGRADDOF_IMPL_HPP
#define CHARON_GRADBASISDOTGRADDOF_IMPL_HPP

#include "Kokkos_ViewFactory.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Integrator_GradBasisDotGradDOF<EvalT, Traits>::
Integrator_GradBasisDotGradDOF(
  const Teuchos::ParameterList& p) :
  residual( p.get<std::string>("Residual Name"),
            p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  flux( p.get<std::string>("Flux Name"),
        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  this->addEvaluatedField(residual);
  this->addDependentField(flux);

  multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"))
  {
    const std::vector<std::string>& field_multiplier_names =
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin();
      name != field_multiplier_names.end(); ++name)
    {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n =
    "Integrator_GradBasisDotGradDOF: " + residual.fieldTag().name();
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_GradBasisDotGradDOF<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  num_nodes = residual.dimension(1);
  num_qp = flux.dimension(1);
  num_dim = flux.dimension(2);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  tmp = Kokkos::createDynRankView(residual.get_static_view(),"tmp",flux.dimension(0), num_qp, num_dim);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Integrator_GradBasisDotGradDOF<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // for (int i=0; i < residual.size(); ++i)
  //   residual[i] = 0.0;
  residual.deep_copy(ScalarT(0.0));

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t qp = 0; qp < num_qp; ++qp)
    {
      ScalarT tmpVar = 1.0;
      for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
           field != field_multipliers.end(); ++field)
        tmpVar = tmpVar * (*field)(cell,qp);

      for (std::size_t dim = 0; dim < num_dim; ++dim)
        tmp(cell,qp,dim) = multiplier * tmpVar * flux(cell,qp,dim);
    }
  }

  if(workset.num_cells>0)
     Intrepid2::FunctionSpaceTools<PHX::exec_space>::
       integrate(residual.get_view(), tmp,
                 (workset.bases[basis_index])->weighted_grad_basis.get_view());
}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Integrator_GradBasisDotGradDOF<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Flux Name", "?");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<double>("Multiplier", 1.0);
  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

}

#endif
