
#ifndef __Charon_ResponseEvaluatorFactory_Current_impl_hpp__
#define __Charon_ResponseEvaluatorFactory_Current_impl_hpp__

#include "Panzer_Normals.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Integrator_Scalar.hpp"

#include "Charon_Scaling_Parameters.hpp"
#include "Charon_FEM_ElectricField.hpp"
#include "Charon_FEM_CurrentDensity.hpp"


namespace charon {

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_Current<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  const panzer::CellData & cd = physicsBlock.cellData();

  TEUCHOS_TEST_FOR_EXCEPTION(!cd.isSide(),std::logic_error,
                     "Charon Error!: Calculating current on a side is required, volume integration is not implemented.");

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = physicsBlock.getParameterList();
  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  // get any prefix or suffix parameters
  std::string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<std::string>("Prefix") : "";
  std::string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<std::string>("Discontinuous Fields") : "";
  std::string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<std::string>("Discontinuous Suffix") : "";
  // append suffix
  names->applySuffixes(discfields, discsuffix);

  // create ir
  RCP<panzer::IntegrationRule> ir = rcp(new panzer::IntegrationRule(this->getCubatureDegree(),physicsBlock.cellData()));

  // create basis
  RCP<const panzer::FieldLibrary> fieldLib = physicsBlock.getFieldLibrary();
  Teuchos::RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(1, prefix, discfields, discsuffix, "_CosH0.000000_"));
  RCP<const panzer::PureBasis> pureBasis = fieldLib->lookupBasis(isFreqDom_ ? fd_names->dof.phi : names->dof.phi);

  RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(pureBasis, *ir));

  // Get scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams = user_data.get<Teuchos::RCP<charon::Scaling_Parameters> >("Scaling Parameter Object");

  // Normals evaluator
  {
    Teuchos::ParameterList p;
    p.set<std::string>("Name","Side Normal");
    p.set<int>("Side ID",cd.side());
    p.set< Teuchos::RCP<panzer::IntegrationRule> >("IR", ir);
    p.set<bool>("Normalize",true);

    RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Computing FEM electron electric field: Fn = grad(Ei/V0-0.5*dEg/V0) at contact IPs,
  // Ei = intrinsic Fermi energy, valid for BGN = On and Off (dEg=0 when BGN = Off).
  if (enableElectrons_)
  {
    ParameterList p("Electron Electric Field");
    p.set("Carrier Type", "Electron");
    p.set< RCP<const charon::Names> >("Names", names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_ElectricField<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Computing electron current density: Jn = n*mun*Fn + Dn*grad(n) at IPs of Ohmic contact
  if (enableElectrons_)
  {
    ParameterList p("Electron Current Density");
    p.set("Carrier Type", "Electron");
    p.set("Current Name", names->field.elec_contact_currdens);
    p.set< RCP<const charon::Names> >("Names", names);
    p.set("IR", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_CurrentDensity<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Computing hole electric field: Fp = grad(Ei/V0 + 0.5*dEg/V0) at contact IPs,
  // Ei = intrinsic Fermi energy, valid for BGN = On and Off (dEg=0 when BGN = Off).
  if (enableHoles_)
  {
    ParameterList p("Hole Electric Field");
    p.set("Carrier Type", "Hole");
    p.set< RCP<const charon::Names> >("Names", names);
    p.set("Basis", basis);
    p.set("IR", ir);
    p.set("Scaling Parameters", scaleParams);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_ElectricField<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // Computing hole current density: Jp = p*mup*Fp - Dp*grad(p) at IPs of Ohmic contact
  if (enableHoles_)
  {
    ParameterList p("Hole Current Density");
    p.set("Carrier Type", "Hole");
    p.set("Current Name", names->field.hole_contact_currdens);
    p.set< RCP<const charon::Names> >("Names", names);
    p.set("IR", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new charon::FEM_CurrentDensity<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  // force a sum over all fields containing the current vector
  std::string currentVectorName = "CurrentVector_" + responseName;
  if(currents_.size()>1)
  {
    Teuchos::ParameterList p("MultiCurrentSum");;
    p.set("Sum Name", currentVectorName);
    p.set("Values Names",Teuchos::rcp_const_cast<std::vector<std::string> >(Teuchos::rcpFromRef(currents_)));
    p.set("Data Layout",ir->dl_vector);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::Sum<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }
  else
    currentVectorName = currents_[0]; // use whatever current vector is needed

  // dot product evaluator
  {
    Teuchos::ParameterList p;
    p.set("Result Name", responseName);
    p.set("Vector A Name","Side Normal");
    p.set("Vector B Name",currentVectorName);
    p.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);

    RCP< PHX::Evaluator<panzer::Traits> > op =
      rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));

    fm.template registerEvaluator<EvalT>(op);
  }

  {
    // obtain the dimensionality of the simulation domain
    int dims = cd.baseCellDimension();

    // scaling factor
    double scaling = scaleParams->scale_params.J0;  // for 1D, contact current in [A/cm^2]
    if (dims == 2)
      scaling *= scaleParams->scale_params.X0;  // for 2D, contact current in [A/cm]
    else if (dims == 3)
      scaling *= scaleParams->scale_params.X0 * scaleParams->scale_params.X0;  // for 3D, contact current in [A]

    Teuchos::ParameterList pl;
    pl.set("Integral Name",responseName);
    pl.set("Integrand Name",responseName);
    pl.set("IR",ir);
    pl.set("Multiplier",scaling);

    Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
        = Teuchos::rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(pl));

    fm.template registerEvaluator<EvalT>(eval);
  }

  panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::buildAndRegisterEvaluators(responseName,fm,physicsBlock,user_data);
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_Current<EvalT,LO,GO>::
typeSupported() const
{
#if 0
  if(PHX::TypeString<EvalT>::value==PHX::TypeString<panzer::Traits::Residual>::value)
    return true;
  if(PHX::TypeString<EvalT>::value==PHX::TypeString<panzer::Traits::Jacobian>::value)
    return true;

  return false;
#else
  return (typeid(EvalT) == typeid(panzer::Traits::Residual) ||
          typeid(EvalT) == typeid(panzer::Traits::Jacobian));
#endif
}

}

#endif
