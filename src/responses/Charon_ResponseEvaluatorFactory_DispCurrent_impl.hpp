
#ifndef __Charon_ResponseEvaluatorFactory_DispCurrent_impl_hpp__
#define __Charon_ResponseEvaluatorFactory_DispCurrent_impl_hpp__

#include "Panzer_Normals.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_DotProduct.hpp"
#include "Panzer_Integrator_Scalar.hpp"

#include "Charon_FEM_ElectricField.hpp"
#include "Charon_FEM_CurrentDensity.hpp"
#include "Charon_Physical_Constants.hpp"

#include "Panzer_DOFGradient.hpp"
#include "Panzer_Constant.hpp"

#include "Charon_DisplacementCurrentOnContact.hpp"


namespace charon {

template <typename EvalT,typename LO,typename GO>
void ResponseEvaluatorFactory_DispCurrent<EvalT,LO,GO>::
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::string;
  using std::vector;  

  const panzer::CellData & cd = physicsBlock.cellData();
  TEUCHOS_TEST_FOR_EXCEPTION(!cd.isSide(),std::logic_error,
                     "Charon Error!: Calculating current on a side is required, volume integration is not implemented.");

  // obtain the dimensionality (1 for 1D, 2 for 2D, 3 for 3D)
  int dims = cd.baseCellDimension();

  // get the physics block parameter list
  RCP<const ParameterList> pbParamList = physicsBlock.getParameterList();

  // allow only one equation set per physics block
  if (pbParamList->numParams() > 1)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
       "The physics block " << pbParamList->name() << " has more than one equation sets ! ");
  
  // get the equation set parameter list
  const ParameterList& eqSetPList = pbParamList->sublist("child0");
  
  // get any prefix or suffix parameters
  string prefix = eqSetPList.isParameter("Prefix") ? eqSetPList.get<std::string>("Prefix") : "";
  string discfields = eqSetPList.isParameter("Discontinuous Fields") ? eqSetPList.get<std::string>("Discontinuous Fields") : "";
  string discsuffix = eqSetPList.isParameter("Discontinuous Suffix") ? eqSetPList.get<std::string>("Discontinuous Suffix") : "";
  RCP<charon::Names> names = Teuchos::rcp(new charon::Names(1, prefix, discfields, discsuffix));

  // append suffix
  names->applySuffixes(discfields, discsuffix);

  // create ir
  RCP<panzer::IntegrationRule> ir = rcp(new panzer::IntegrationRule(this->getCubatureDegree(),physicsBlock.cellData()));

  // create basis
  RCP<const panzer::FieldLibrary> fieldLib = physicsBlock.getFieldLibrary();
  RCP<charon::Names> fd_names = Teuchos::rcp(new charon::Names(1, prefix, discfields, discsuffix, "_CosH"+std::to_string(0.0)+"_"));
  RCP<const panzer::PureBasis> pureBasis = fieldLib->lookupBasis(isFreqDom_ ? fd_names->dof.phi : names->dof.phi);
  RCP<panzer::BasisIRLayout> basis = rcp(new panzer::BasisIRLayout(pureBasis, *ir));

  // obtain the scaling factor
  double scaling = scaleParams_->scale_params.E0;  // Electric field scaling factor in [V/cm]
  if (dims == 2)
    scaling *= scaleParams_->scale_params.X0;  // for 2D, scaling factor = E0*X0
  else if (dims == 3)
    scaling *= scaleParams_->scale_params.X0 * scaleParams_->scale_params.X0;  // for 3D, scaling factor = E0*X0*X0

  // get physical constants
  const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  double eps0 = phyConst.eps0;  // Vacuum permittivity in [C/(V.cm)]
  
  // define string names
  string side_normal = "Side Normal";
  string cos_suffix = "_CosH"+std::to_string(1.0)+"_";  //1.0 for small signal harmonic
  string sin_suffix = "_SinH"+std::to_string(1.0)+"_";
  string grad_phi_cos = names->grad_dof.phi + cos_suffix;  
  string grad_phi_sin = names->grad_dof.phi + sin_suffix;
  string normal_gradphi_cos = "Normal Dot GradPhi CosH";
  string normal_gradphi_sin = "Normal Dot GradPhi SinH";
  
  // analyze fd_suffix_ and set the right calculation flag
  bool isDispCurrCosH = false;
  bool isDispCurrSinH = false; 
  if (cos_suffix.compare(fd_suffix_) == 0) isDispCurrCosH = true; 
  if (sin_suffix.compare(fd_suffix_) == 0) isDispCurrSinH = true; 
  
  // get the normal vector of a sideset at IP from the "Normals" evaluator
  {
    Teuchos::ParameterList p;
    p.set<std::string>("Name", side_normal);
    p.set<int>("Side ID",cd.side());
    p.set< Teuchos::RCP<panzer::IntegrationRule> >("IR", ir);
    p.set<bool>("Normalize",true);

    RCP< PHX::Evaluator<panzer::Traits> > op = rcp(new panzer::Normals<EvalT,panzer::Traits>(p));
    fm.template registerEvaluator<EvalT>(op);
  }

  // compute displacement current for single-frequency FD analysis
  if (isSingleFreq_) 
  {
   // add eps0 and omega (radian frequency) to the "scaling" factor
   scaling *= eps0 * omega_;   // in [A/cm] for 2D and [A] for 3D 

   // compute displacement current for the sin(wt) component
   if (isDispCurrSinH) 
   { 
    // compute the side normal dot the potential gradient cos-component
    {
      ParameterList p(normal_gradphi_cos);
      p.set("Result Name", normal_gradphi_cos);
      p.set("Vector A Name", grad_phi_cos);
      p.set("Vector B Name", side_normal);
      p.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);

      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }

    // integrate normal_gradphi_cos over a contact with the proper scaling factor
    {
      scaling *= 1.0;  
      ParameterList p;
      p.set("Integral Name",responseName);
      p.set("Integrand Name",normal_gradphi_cos);
      p.set("IR",ir);
      p.set("Multiplier",scaling);

      RCP<vector<string> > fms = Teuchos::rcp(new vector<string>);
      fms->push_back(names->field.rel_perm+"_TP0_");
      p.set< RCP<const vector<string> > >("Field Multipliers",fms);

      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(eval);
    }
   }  // end of if (isDispCurrSinH) 
   
   // compute displacement current for the cos(wt) component
   else if (isDispCurrCosH)
   { 
    // compute the side normal dot the potential gradient sin-component
    {
      ParameterList p(normal_gradphi_sin);
      p.set("Result Name", normal_gradphi_sin);
      p.set("Vector A Name", grad_phi_sin);
      p.set("Vector B Name", side_normal);
      p.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);

      const RCP< PHX::Evaluator<panzer::Traits> >
        op = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, op);
    }

    // integrate normal_gradphi_sin over a contact with the proper scaling factor
    {
      scaling *= -1.0;  
      ParameterList p;
      p.set("Integral Name",responseName);
      p.set("Integrand Name",normal_gradphi_sin);
      p.set("IR",ir);
      p.set("Multiplier",scaling);

      RCP<vector<string> > fms = Teuchos::rcp(new vector<string>);
      fms->push_back(names->field.rel_perm+"_TP0_");
   
      p.set< RCP<const vector<string> > >("Field Multipliers",fms);

      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(eval);
    }
   }  // end of else if (isDispCurrCosH)
     
  }  // end of if (isSingleFreq_) 

  // compute displacement current for multiple-frequency FD analysis
  else if (isFreqDom_ && !isSingleFreq_) 
  {
    // TODO:
  }

  // compute displacement current for transient simulation
  else if (isTransient_) 
  {
    string disp_curr = names->field.cont_disp_curr_density;
    string normal_disp_curr = "Normal Dot DispCurrent";
    RCP<charon::Scaling_Parameters> scaleParams = 
	user_data.get<RCP<charon::Scaling_Parameters>>("Scaling Parameter Object");

    {
      ParameterList p("Cont Disp Current Density");
      p.set("Current Name", disp_curr);
      p.set<Teuchos::RCP<const charon::Names> >("Names", names);
      p.set("Scaling Parameters", scaleParams);
      p.set("IR", ir);
      Teuchos::RCP< PHX::Evaluator<panzer::Traits> > eval = 
	Teuchos::rcp(new charon::DisplacementCurrentOnContact<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm,eval);
    }

    {
      ParameterList p(normal_disp_curr);
      p.set("Result Name", normal_disp_curr);
      p.set("Vector A Name", disp_curr);
      p.set("Vector B Name", side_normal);
      p.set<Teuchos::RCP<const panzer::PointRule> >("Point Rule", ir);

      const RCP< PHX::Evaluator<panzer::Traits> >
	eval = rcp(new panzer::DotProduct<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, eval);
    }

    // integrate normal_disp_curr over a contact with the proper scaling factor
    {
        
      // scaling factor
      double sc = scaleParams->scale_params.J0;  // for 1D, contact current in [A/cm^2]
      if (dims == 2)
	sc *= scaleParams->scale_params.X0;  // for 2D, contact current in [A/cm]
      else if (dims == 3)
	sc *= scaleParams->scale_params.X0 * scaleParams->scale_params.X0;  // for 3D, contact current in [A]
      
      ParameterList p;
      p.set("Integral Name",responseName);
      p.set("Integrand Name",normal_disp_curr);
      p.set("IR",ir);
      p.set("Multiplier",sc);

      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
          = Teuchos::rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(p));
      fm.template registerEvaluator<EvalT>(eval);
    }
  }

  panzer::ResponseEvaluatorFactory_Functional<EvalT,LO,GO>::buildAndRegisterEvaluators(responseName,fm,physicsBlock,user_data);
}

template <typename EvalT,typename LO,typename GO>
bool ResponseEvaluatorFactory_DispCurrent<EvalT,LO,GO>::
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
