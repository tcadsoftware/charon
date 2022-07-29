#ifndef CHARON_RELATIVE_PERMITTIVITY_IMPL_HPP
#define CHARON_RELATIVE_PERMITTIVITY_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"



namespace charon {


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Relative_Permittivity<EvalT, Traits>::
Relative_Permittivity(const Teuchos::ParameterList& p) 
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

  // obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  withMoleFrac = matProperty.hasMoleFracDependence(materialName);
  if (withMoleFrac)
    comp_mat = matProperty.getMoleFracMaterial(materialName);
  
  // Evaluated fields
  rel_perm = MDField<ScalarT,Cell,Point>(n.field.rel_perm,scalar);
  this->addEvaluatedField(rel_perm);
  
  // Dependent fields
  if (withMoleFrac) {
    xMoleFrac = MDField<const ScalarT,Cell,Point>(n.field.xMoleFrac,scalar);
    this->addDependentField(xMoleFrac);
    if (matProperty.getArityType(materialName) == "Quaternary") {
      yMoleFrac = MDField<const ScalarT,Cell,Point>(n.field.yMoleFrac,scalar);
      this->addDependentField(yMoleFrac);
    }
  }

  std::string name = "Relative_Permittivity";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Relative_Permittivity<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();
  using panzer::index_t;

  // loop over the cells 
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (int point = 0; point < num_points; ++point) {
      if (withMoleFrac) {
	const std::string& arity = matProperty.getArityType(materialName);
	ScalarT T = 300.0;
	if (arity == "Binary" or arity == "Ternary") {
	  ScalarT xFrac = xMoleFrac(cell,point);
	  rel_perm(cell,point) = comp_mat->compute_Eps<EvalT>(
		    Sacado::ScalarValue<ScalarT>::eval(xFrac),0.0,T);
	} else { // "Quaternary"
	  ScalarT xFrac = xMoleFrac(cell,point);
	  ScalarT yFrac = yMoleFrac(cell,point);
	  rel_perm(cell,point) = comp_mat->compute_Eps<EvalT>(
		    Sacado::ScalarValue<ScalarT>::eval(xFrac),
		    Sacado::ScalarValue<ScalarT>::eval(yFrac),
		    T);
	}
      } else {
	;
      }
    } // points
  } // cells 
}



///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Relative_Permittivity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Relative Permittivity ParameterList", false, "");
  p->sublist("Relative Permittivity ParameterList").set<double>(
			     "Value", 0.0, "Relative Permittivity");
  Teuchos::ParameterList& moleFracParams = 
    p->sublist("Relative Permittivity ParameterList").sublist(
			    "Mole Fraction Parameters", false, "");
  Teuchos::ParameterList& Value = moleFracParams.sublist("Value", false, "");
  Value.set<double>("b", 0., "Relative Permittivity mole fraction 'b' interpolation coefficient");
  Value.set<double>("c", 0., "Relative Permittivity mole fraction 'c' interpolation coefficient");
  moleFracParams.set<double>("Value(x=0)", 0., "Relative Permittivity for x=0");
  moleFracParams.set<double>("Value(x=1)", 0., "Relative Permittivity for x=1");

  return p;
}

//**********************************************************************

}

#endif
