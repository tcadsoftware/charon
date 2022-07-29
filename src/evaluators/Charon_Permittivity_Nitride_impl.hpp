#ifndef CHARON_PERMITTIVITY_NITRIDE_IMPL_HPP
#define CHARON_PERMITTIVITY_NITRIDE_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
 * Computes permittivity of nitrides as a function of mole fraction
 * O.Ambacher et al. Appl.Phys. Vol87 No1 2000
 */

namespace charon {

//**********************************************************************
PHX_EVALUATOR_CTOR(Permittivity_Nitride,p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  auto scalar = p.get< RCP<DataLayout> >("Data Layout"); 
  num_points = scalar->dimension(1);

  // Material name
  materialName = p.get<string>("Material Name");
  TEUCHOS_TEST_FOR_EXCEPTION(((materialName !="AlGaN") and 
    (materialName !="InGaN")),
      std::logic_error, "Invalid ternary material: "<< materialName << "!" 
        << std::endl);
  
  // Evaluated fields
  rel_perm = MDField<ScalarT,Cell,Point>(n.field.rel_perm,scalar);

  this->addEvaluatedField(rel_perm);
  
  // Dependent fields
  molefrac = MDField<const ScalarT,Cell,Point>(n.field.mole_frac,scalar);

  this->addDependentField(molefrac);

  std::string name = "Permittivity_Nitride";
  this->setName(name);
}


//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Permittivity_Nitride,sd,fm)
{
  this->utils.setFieldData(rel_perm,fm);
  
  this->utils.setFieldData(molefrac,fm);
}


//**********************************************************************
PHX_EVALUATE_FIELDS(Permittivity_Nitride,workset)
{
  using panzer::index_t;

  ScalarT perm{0.0};

  // loop over the cells 
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT x = molefrac(cell,point); 

      if (materialName == "AlGaN")
        perm = 8.5*x + 8.9*(1 - x);
      if (materialName == "InGaN")
        perm = 15.3*x + 8.9*(1 - x);

      rel_perm(cell,point) = perm; 
    }
  }  
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Permittivity_Nitride<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Relative Permittivity ParameterList", false, "");
  p->sublist("Relative Permittivity ParameterList").set<std::string>("Value", "Nitride", "Specify nitride permittivity");

  return p;
}

//**********************************************************************

}

#endif
