
#ifndef CHARON_EQUATIONSET_SGCHARON1_DRIFTDIFFUSION_IMPL_HPP
#define CHARON_EQUATIONSET_SGCHARON1_DRIFTDIFFUSION_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_DOF.hpp"
#include "Panzer_DOFGradient.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Panzer_Sum.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_PoissonSource.hpp"
#include "Charon_SGCVFEM_EdgeCurrDens.hpp"
#include "Charon_SGCharon1_SubCVCurrDens.hpp"
#include "Charon_SGCharon1_PotentialFlux.hpp"
#include "Charon_Integrator_SubCVNodeScalar.hpp"
#include "Charon_Integrator_SubCVFluxDotNorm.hpp"


// CVFEM-SG formulation of the Poisson + Continuity equations

// ***********************************************************************
template <typename EvalT>
charon::EquationSet_SGCharon1_DriftDiffusion<EvalT>::
EquationSet_SGCharon1_DriftDiffusion(const Teuchos::RCP<Teuchos::ParameterList>& params,
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

    Teuchos::ParameterList& opt = valid_parameters.sublist("Options");
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Electron",
      "True",
      "Determines if users want to solve the Electron Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Solve Hole",
      "True",
      "Determines if users want to solve the Hole Continuity equation",
      Teuchos::tuple<std::string>("True","False"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "SRH",
      "Off",
      "Determines if users want to include SRH recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Radiative",
      "Off",
      "Determines if users want to include Radiative recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Auger",
      "Off",
      "Determines if users want to include Auger recombination in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Avalanche",
      "Off",
      "Determines if users want to include Avalanche generation in the continuity equation(s)",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Band Gap Narrowing",
      "Off",
      "Determines if users want to turn on BGN effect due to heavy doping",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Acceptor Incomplete Ionization",
      "Off",
      "Determine if users want to include Acceptor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );
    Teuchos::setStringToIntegralParameter<int>(
      "Donor Incomplete Ionization",
      "Off",
      "Determine if users want to include Donor Incomplete Ionization",
      Teuchos::tuple<std::string>("On","Off"),
      &opt
      );

    params->validateParametersAndSetDefaults(valid_parameters);

    solveElectron = params->sublist("Options").get<std::string>("Solve Electron");
    solveHole = params->sublist("Options").get<std::string>("Solve Hole");

    // Must solve at least one continuity equation
    //if ((solveElectron == "False") && (solveHole == "False"))
    //  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error: must solve at least one"
    //       << " continuity equation !" << "\n");

    const std::string& srh_recomb = params->sublist("Options").get<std::string>("SRH");
    const std::string& rad_recomb = params->sublist("Options").get<std::string>("Radiative");
    const std::string& auger_recomb = params->sublist("Options").get<std::string>("Auger");
    const std::string& ava_gen = params->sublist("Options").get<std::string>("Avalanche");

    // Must solve both continuity equations when any recomb./gen. model is on.
    if ((srh_recomb == "On") || (rad_recomb == "On") || (auger_recomb == "On") || (ava_gen == "On"))
    {
      TEUCHOS_ASSERT((solveElectron == "True") && (solveHole == "True"));
      haveSource = true;
    }
    else
      haveSource = false;  // NO recomb. model

    // error check for Incomplete Ionization models
    const std::string& acc_ioniz =
      params->sublist("Options").get<std::string>("Acceptor Incomplete Ionization");
    const std::string& don_ioniz =
      params->sublist("Options").get<std::string>("Donor Incomplete Ionization");
    if(acc_ioniz == "On" && solveHole == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve hole continuity equation when " <<
          "Acceptor Incomplete Ionization is On!" << "\n");
    if(don_ioniz == "On" && solveElectron == "False")
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "Error: must solve electron continuity equation when " <<
          "Donor Incomplete Ionization is On!" << "\n");
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

  // *************************************
  // Assemble DOF names and Residual names
  // *************************************
  m_names = Teuchos::rcp(new charon::Names(cell_data.baseCellDimension(),prefix,discfields,discsuffix));

  this->getEvaluatorParameterList()->set("Names",Teuchos::RCP<const charon::Names>(m_names));

  // Always solve the Poisson equation
  this->addDOF(m_names->dof.phi,basis_type,basis_order,integration_order,m_names->res.phi);
  this->addDOFGrad(m_names->dof.phi,m_names->grad_dof.phi);

  // Solve the Electron equation by default and users can disable it
  if (solveElectron == "True")
  {
    this->addDOF(m_names->dof.edensity,basis_type,basis_order,integration_order,m_names->res.edensity);
    this->addDOFGrad(m_names->dof.edensity,m_names->grad_dof.edensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.edensity,m_names->dxdt.edensity);
  }

  // Solve the Hole equation by default and users can disable it
  if (solveHole == "True")
  {
    this->addDOF(m_names->dof.hdensity,basis_type,basis_order,integration_order,m_names->res.hdensity);
    this->addDOFGrad(m_names->dof.hdensity,m_names->grad_dof.hdensity);
    if (this->buildTransientSupport())
      this->addDOFTimeDerivative(m_names->dof.hdensity,m_names->dxdt.hdensity);
  }

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void charon::EquationSet_SGCharon1_DriftDiffusion<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                      const panzer::FieldLibrary& /* fl */,
                                      const Teuchos::ParameterList& user_data) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const charon::Names& n = *m_names;

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // For now, assume that all dofs use the same basis
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(m_names->dof.phi);
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(m_names->dof.phi);

  // Get basis and IR layout for CVFEM
  panzer::CellData celldata(basis->numCells(),ir->topology);
  std::string cvfem_type = "volume";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_vol_ir = Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));
  cvfem_type = "side";
  Teuchos::RCP<panzer::IntegrationRule> cvfem_side_ir = Teuchos::rcp(new panzer::IntegrationRule(celldata,cvfem_type));

  Teuchos::RCP<const panzer::PureBasis> Hgradbasis = Teuchos::rcp(new panzer::PureBasis("HGrad",1,basis->numCells(),ir->topology));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_vol_cvfem = Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_vol_ir));
  Teuchos::RCP<panzer::BasisIRLayout> hgrad_side_cvfem = Teuchos::rcp(new panzer::BasisIRLayout(Hgradbasis,*cvfem_side_ir));


  // ***************************************************************************
  // Construct Poisson Equation (solved with continuity equations) based on the
  // Charon1-SG formulation implemented in Charon1
  // ***************************************************************************

  // Laplacian Operator: \int_cv_boundary (-\lambda2 \epsilon_r \grad_phi) \cdot {\bf n} dS
  {
    // Compute the potential flux: lambda2*epsilon_r*grad_phi @ the midpoints
    // of primary edges using finite difference
    {
      ParameterList p("Charon1-SG Potential Flux");
      p.set("Flux Name", n.field.phi_flux);
      p.set("DOF Name", n.dof.phi);
      p.set("Basis", basis);
      p.set("Scaling Parameters", scaleParams);
      p.set< RCP<const charon::Names> >("Names", m_names);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::SGCharon1_PotentialFlux<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Integrate the Charon1-SG potential flux over the subcv boundary
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

  // Source Operator: Assembles \int_cv g dV
  {
    // Compute the Poisson source term: (p-n+Nd-Na) @ BASIS
    {
      ParameterList p("Poisson Source");
      p.set("Source Name", n.field.psrc);
      p.set("Data Layout", basis->functional);
      p.set< RCP<const charon::Names> >("Names", m_names);
      p.set("Solve Electron", solveElectron);
      p.set("Solve Hole", solveHole);
      p.set("Scaling Parameters", scaleParams);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::PoissonSource<EvalT,panzer::Traits>(p));

      fm.template registerEvaluator<EvalT>(op);
    }

    // Interpolate the source term at primary nodes to the center point of a
    // subcontrol volume and integrate over subcv
    {
      ParameterList p("Source Residual");
      p.set("Residual Name", n.res.phi+n.op.src);
      p.set("Value Name", n.field.psrc);
      p.set("Basis", hgrad_vol_cvfem);
      p.set("IR", cvfem_vol_ir);
      p.set<RCP<const charon::Names> >("Names", m_names);
      p.set("Multiplier", -1.0);
      p.set("WithInterpolation", true);

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
    sum_names->push_back(n.res.phi+n.op.src);

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }


  // ***************************************************************************
  // Construct Electron Continuity Equation based on the Charon1-SG formulation
  // ***************************************************************************

 if (solveElectron == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    // integrate dn_dt over subcv, \int_cv dn_dt dV, where dn_dt = \sum_BFi dni_dt BFi
    // interpolate dn_dt at primary nodes to the center point of subcv and integrate over subcv
    ParameterList p("Electron Transient Residual");
    p.set("Residual Name", n.res.edensity+n.op.trans);
    p.set("Value Name", n.dxdt.edensity);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute electron CVFEM-SG edge current density (scalar) at the midpoints of primary edges
   {
    ParameterList p("Electron CVFEM-SG Edge Current Density");
    p.set("Carrier Type", "Electron");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_EdgeCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // compute electron Charon1-SG current density (vector) at the midpoints of primary edges
   {
    ParameterList p("Electron Charon1-SG SubCV Current Density");
    p.set("Carrier Type", "Electron");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCharon1_SubCVCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // integrate electron Charon1-SG current density over the subcv boundary
   {
    ParameterList p("Electron Convection Diffusion Residual");
    p.set("Residual Name", n.res.edensity+n.op.conv_diff);
    p.set("Flux Name", n.field.elec_curr_dens_cvedge);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", cvfem_side_ir);
    p.set("Multiplier", -1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   // interpolate the total source term (exclude the Avalanche gen.) at primary nodes
   // to the center point of subcontrol volume and integrate over subcv
   if (haveSource)
   {
    ParameterList p("Electron Source Residual");
    p.set("Residual Name", n.res.edensity+n.op.src);
    p.set("Value Name", n.field.total_recomb);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Use a sum operator to form the overall residual for the equation
  {
    ParameterList p;
    p.set("Sum Name", n.res.edensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.edensity + n.op.conv_diff);

    if (this->buildTransientSupport())
      sum_names->push_back(n.res.edensity + n.op.trans);

    if (haveSource)
      sum_names->push_back(n.res.edensity + n.op.src);

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveElectron == "True")


  // ***************************************************************************
  // Construct Hole Continuity Equation based on the Charon1-SG formulation
  // ***************************************************************************

 if (solveHole == "True")
 {
  // Transient Residual
  if (this->buildTransientSupport())
  {
    // integrate dp_dt over subcv, \int_cv dp_dt dV, where dp_dt = \sum_BFi dpi_dt BFi
    // interpolate dp_dt at primary nodes to the center point of subcv and integrate over subcv
    ParameterList p("Hole Transient Residual");
    p.set("Residual Name", n.res.hdensity+n.op.trans);
    p.set("Value Name", n.dxdt.hdensity);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Convection-Diffusion Residual
  {
   // compute hole CVFEM-SG edge current density (scalar) at the midpoints of primary edges
   {
    ParameterList p("Hole CVFEM-SG Edge Current Density");
    p.set("Carrier Type", "Hole");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCVFEM_EdgeCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // compute hole Charon1-SG current density (vector) at the midpoints of primary edges
   {
    ParameterList p("Hole Charon1-SG SubCV Current Density");
    p.set("Carrier Type", "Hole");
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::SGCharon1_SubCVCurrDens<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }

   // integrate hole Charon-SG current density over the subcv boundary
   {
    ParameterList p("Hole Convection Diffusion Residual");
    p.set("Residual Name", n.res.hdensity+n.op.conv_diff);
    p.set("Flux Name", n.field.hole_curr_dens_cvedge);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Basis", basis);
    p.set("IR", cvfem_side_ir);
    p.set("Multiplier", 1.0);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVFluxDotNorm<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Source Residual
  {
   // interpolate the total source term (exclude the Avalanche gen.) at primary nodes
   // to the center point of subcontrol volume and integrate over subcv
   if (haveSource)
   {
    ParameterList p("Hole Source Residual");
    p.set("Residual Name", n.res.hdensity+n.op.src);
    p.set("Value Name", n.field.total_recomb);
    p.set("Basis", hgrad_vol_cvfem);
    p.set("IR", cvfem_vol_ir);
    p.set<RCP<const charon::Names> >("Names", m_names);
    p.set("Multiplier", 1.0);
    p.set("WithInterpolation", true);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::Integrator_SubCVNodeScalar<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
   }
  }

  // Use a sum operator to form the overall residual for the equation
  {
    ParameterList p;
    p.set("Sum Name", n.res.hdensity);

    RCP<std::vector<std::string> > sum_names =
      rcp(new std::vector<std::string>);

    sum_names->push_back(n.res.hdensity + n.op.conv_diff);

    if (this->buildTransientSupport())
      sum_names->push_back(n.res.hdensity + n.op.trans);

    if (haveSource)
      sum_names->push_back(n.res.hdensity+ n.op.src);

    p.set("Values Names", sum_names);
    p.set("Data Layout", basis->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

 }  // end of if (solveHole == "True")


}

// ***********************************************************************

#endif
