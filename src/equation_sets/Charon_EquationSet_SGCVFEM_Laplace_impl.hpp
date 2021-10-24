
#ifndef CHARON_EQUATIONSET_SGCVFEM_LAPLACE_IMPL_HPP
#define CHARON_EQUATIONSET_SGCVFEM_LAPLACE_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_Sum.hpp"

#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PotentialFlux.hpp"
#include "Charon_SGCVFEM_PotentialFlux.hpp"
#include "Charon_Integrator_SubCVFluxDotNorm.hpp"
#include "Charon_Integrator_SubCVNodeScalar.hpp"


// ***********************************************************************
template <typename EvalT>
charon::EquationSet_SGCVFEM_Laplace<EvalT>::
EquationSet_SGCVFEM_Laplace(const Teuchos::RCP<Teuchos::ParameterList>& params,
                            const int& default_integration_order,
                            const panzer::CellData& cell_data,
                            const Teuchos::RCP<panzer::GlobalData>& global_data,
                            const bool build_transient_support) :
  charon::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support)
{
  // ********************
  // Options
  // ********************
  {
    Teuchos::ParameterList valid_parameters;

    this->setDefaultValidParameters(valid_parameters);

    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("Prefix","","Prefix for using multiple instantiations of the equation set");
    valid_parameters.set("Discontinuous Fields","","List of fields which are discontinuous");
    valid_parameters.set("Discontinuous Suffix","","Suffix for enabling discontinuous fields");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    // valid_parameters.sublist("Options");
    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");
    Teuchos::setStringToIntegralParameter<int>(
      "Fixed Charge",
      "False",
      "Determine if users want to add fixed charges in an insulator region",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "TID",
      "Off",
      "Determine if users want to add TID models in an insulator region",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    // determine if a user wants to include bulk fixed charges
    addFixCharge = false; 
    if (params->sublist("Options").get<std::string>("Fixed Charge") == "True")
      addFixCharge = true; 
    addTID = false; 
    if (params->sublist("Options").get<std::string>("TID") == "On")
      addTID = true; 
  }

  std::string prefix = params->get<std::string>("Prefix");
  std::string discfields = params->get<std::string>("Discontinuous Fields");
  std::string discsuffix = params->get<std::string>("Discontinuous Suffix");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  std::string model_id = params->get<std::string>("Model ID");
  int integration_order = params->get<int>("Integration Order");

  this->getEvaluatorParameterList()->sublist("Options") = params->sublist("Options");
  this->getEvaluatorParameterList()->set("Type",params->get<std::string>("Type"));

  // ***************************************************************************
  // Assemble DOF names and Residual names
  // ***************************************************************************
  m_names = Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix));

  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_names));

  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);
  if (this->buildTransientSupport())
    this->addDOFTimeDerivative(m_names->dof.phi,m_names->dxdt.phi);

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_SGCVFEM_Laplace<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<charon::Names> m_names = (this->isFreqDom ? this->base_names : this->m_names);
  const charon::Names& n = (this->isFreqDom ? *(this->base_names) : *(this->m_names));

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  Teuchos::RCP<panzer::IntegrationRule> ir = (this->isFreqDom ? this->base_ir : this->getIntRuleForDOF(m_names->dof.phi));
  Teuchos::RCP<panzer::BasisIRLayout> basis = (this->isFreqDom ? this->base_basis : this->getBasisIRLayoutForDOF(m_names->dof.phi));

  // Get basis and IR layout for CVFEM
  panzer::CellData celldata(basis->numCells(),ir->topology);
  std::string cvfem_type = "side";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_side_ir = Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));
  Teuchos::RCP<const panzer::PureBasis> Hgradbasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),ir->topology));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_side_cvfem = Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_side_ir));

  cvfem_type = "volume";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_vol_ir =
    Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_vol_cvfem =
    Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_vol_ir));


  // ***************************************************************************
  // Construct the Laplace equation based on the CVFEM formulation
  // ***************************************************************************

  // Laplacian Operator: \int_cv_boundary (-\lambda2 \epsilon_r \grad_phi) \cdot {\bf n} dS
  {
    // Compute the potential flux: lambda2*epsilon_r*grad_phi @ the midpoints
    // of subcv edges(2D)/faces(3D) using nodal basis functions (scalar)
    {
      ParameterList p("CVFEM-SG Potential Flux");
      p.set("Flux Name", n.field.phi_flux);
      p.set("DOF Name", n.dof.phi);
      p.set("Basis", hgrad_side_cvfem);
      p.set("Scaling Parameters", scaleParams);
      p.set< RCP<const charon::Names> >("Names", m_names);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::SGCVFEM_PotentialFlux<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the CVFEM-SG potential flux over the subcv boundary
    {
      ParameterList p("Laplacian Residual");
      p.set("Residual Name", n.res.phi+n.op.laplacian);
      p.set("Flux Name", n.field.phi_flux);
      p.set<RCP<const charon::Names> >("Names", m_names);
      p.set("Basis", basis);
      p.set("IR", cvfem_side_ir);
      p.set("Multiplier", -1.0);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Add the fixed charges to the equation residual
  if (this->isFreqDom ? this->base_addFixCharge : this->addFixCharge)
  {
    ParameterList p("FixCharge Source Residual");
    p.set("Residual Name", n.res.phi+"_FIXCHARGE_SOURCE_OP");
    p.set("Value Name", n.field.fixed_charge);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set("Multiplier", -1.0);
    p.set("WithInterpolation", false);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // Add the TID trapped hole charges to the equation residual
  if (this->addTID && !this->isFreqDom)
  {
    {
      ParameterList p("Ins Trapped Charge Residual");
      p.set("Residual Name", n.res.phi+"_INSTRAPPEDCHARGE_OP");
      p.set("Value Name", n.field.ins_htrappedcharge);
      p.set("Basis", hgrad_vol_cvfem);
      p.set("IR", cvfem_vol_ir);
      p.set("Multiplier", -1.0);
      p.set("WithInterpolation", false);

      RCP< PHX::Evaluator<panzer::Traits> > op =
	rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
  }

  // Use a sum operator to form the overall residual for the equation
  // - this way we avoid loading each operator separately into the
  // global residual and Jacobian
  {
    ParameterList p;
    p.set("Sum Name", n.res.phi);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);
    sum_names->push_back(n.res.phi+n.op.laplacian);

    if (this->isFreqDom ? this->base_addFixCharge : this->addFixCharge) 
      sum_names->push_back(n.res.phi + "_FIXCHARGE_SOURCE_OP");

    if (this->addTID && !this->isFreqDom)
      sum_names->push_back(n.res.phi + "_INSTRAPPEDCHARGE_OP");

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

}

// ***********************************************************************

#endif
