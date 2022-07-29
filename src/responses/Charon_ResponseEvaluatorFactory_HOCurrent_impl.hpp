
#ifndef __Charon_ResponseEvaluatorFactory_HOCurrent_impl_hpp__
#define __Charon_ResponseEvaluatorFactory_HOCurrent_impl_hpp__

#include "Panzer_Normals.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Panzer_SubcellSum.hpp"

#include "Charon_Subtract.hpp"


namespace charon {

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_HOCurrent<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::string;

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = physicsBlock.getParameterList();

  // allow only one equation set per physics block
  if (pbParamList->numParams() > 1)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
       "The physics block " << pbParamList->name() << " has more than one equation sets ! ");

  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  const ParameterList& options = eqSetPList.sublist("Options");

  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<string>("Prefix") : "";
  std::string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<std::string>("Discontinuous Fields") : "";
  std::string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<std::string>("Discontinuous Suffix") : "";
  // append suffix
  names->applySuffixes(discfields, discsuffix);

  // determine if the electron and/or hole equations are solved
  bool solveElectron = false;
  bool solveHole = false;

  // "Solve Electron" and "Solve Hole" do not exist for the "Lattice" equation set.
  // So options.get<string>("Solve Electron") directly will produce an error.
  if (options.isParameter("Solve Electron"))
    if (options.get<string>("Solve Electron") == "True") solveElectron = true;
  if (options.isParameter("Solve Hole"))
    if (options.get<string>("Solve Hole") == "True") solveHole = true;

  // obtain the residual name
  string eResidualName = "";
  string hResidualName = "";
  if (solveElectron) eResidualName = names->res.edensity;
  if (solveHole) hResidualName = names->res.hdensity;

  // create basis 
  RCP<const panzer::FieldLibrary> fieldLib = physicsBlock.getFieldLibrary();
  Teuchos::RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(1, prefix, discfields, discsuffix, "_CosH0.000000_"));

  // obtain the dimensionality (1 for 1D, 2 for 2D, 3 for 3D)
  const panzer::CellData & cdata = physicsBlock.cellData();
  int dims = cdata.baseCellDimension();

  // obtain the scaling factor
  double scaling = scaleParams_->scale_params.J0;  // for 1D, contact current in [A/cm^2]
  if (dims == 2)
    scaling *= scaleParams_->scale_params.X0;  // for 2D, contact current in [A/cm]
  else if (dims == 3)
    scaling *= scaleParams_->scale_params.X0 * scaleParams_->scale_params.X0;  // for 3D, contact current in [A]

  // Limit the following calculation to the case of solveElectron or solveHole = true.
  // If none of them is true, should not compute the terminal current.
  // An example is the "Lattice" equation set, where solveElectron = solveHole = false.
  
  if (solveElectron || solveHole)
  {
    // create basis
    RCP<const panzer::FieldLibrary> fieldLib = physicsBlock.getFieldLibrary();
    RCP<const panzer::PureBasis> basis;
    if (solveElectron) basis = fieldLib->lookupBasis(isFreqDom_ ? fd_names->dof.edensity : names->dof.edensity);
    if (solveHole) basis = fieldLib->lookupBasis(isFreqDom_ ? fd_names->dof.hdensity : names->dof.hdensity);

    // Note eCurrent = \int (res.) dV, while hCurrent = -\int (res.) dV, hence
    // we need (eResidual - hResidual) to get the total current
    {
      Teuchos::ParameterList p;
      p.set("Difference Name", responseName+"_residual_diff");
      p.set("Value A", eResidualName);
      p.set("Value B", hResidualName);
      // p.set("Data Layout",eBasis->functional);
      p.set("Data Layout",basis->functional);
      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new charon::Subtract<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }

    // sum over the sub cells on the edge (use subcell_dim and subcell_id)
    {
      Teuchos::ParameterList pl;
      pl.set("Sum Name",responseName);
      pl.set("Field Name",responseName+"_residual_diff");
      // pl.set("Basis",eBasis);
      pl.set("Basis",basis);
      pl.set("Multiplier",scaling);
      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer::SubcellSum<EvalT,panzer::Traits>(pl));
      fm.template registerEvaluator<EvalT>(eval);
    }

    panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::buildAndRegisterEvaluators(responseName,fm,physicsBlock,user_data);
  }
  
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_HOCurrent<EvalT,LO,GO>::
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
