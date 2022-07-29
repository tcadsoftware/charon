#ifndef CHARON_BANDGAP_NITRIDE_IMPL_HPP
#define CHARON_BANDGAP_NITRIDE_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
 *Computes Bandgap for nitrides as with respect to mole fraction
 *R.Passler Phys.Stat.Solidi(b) 216 (1999):975
 */

namespace charon {

//**********************************************************************
PHX_EVALUATOR_CTOR(BandGap_Nitride,p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout"); 
  num_points = scalar->dimension(1);

  // material name
  materialName = p.get<string>("Material Name");
  isBinary = false;
  isTernary = false;

  // bandgap parameterlist
  //const ParameterList& bgParamList = p.sublist("Bandgap ParameterList");

  // obtain the instance of charon::Material_Properties.
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();
  auto arity = matProperty.getArityType(materialName);
  if(arity == "Binary")
    isBinary = true;
  if(arity == "Ternary")
    isTernary = true;

  // scaling parameter
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  
  // Evaluated fields
  band_gap = MDField<ScalarT,Cell,Point>(n.field.band_gap,scalar);

  // Dependent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);
  if(isTernary)
    molefrac = MDField<const ScalarT,Cell,Point>(n.field.mole_frac,scalar);

  // evaluated fields
  this->addEvaluatedField(band_gap);
  
  // dependent fields
  this->addDependentField(latt_temp);
  if(isTernary)
    this->addDependentField(molefrac);

  std::string name = "BandGap_Nitride";
  this->setName(name);
}


//**********************************************************************
PHX_POST_REGISTRATION_SETUP(BandGap_Nitride,sd,fm)
{
  this->utils.setFieldData(band_gap,fm);

  this->utils.setFieldData(latt_temp,fm);
  if(isTernary)
    this->utils.setFieldData(molefrac,fm);
}


//**********************************************************************
PHX_EVALUATE_FIELDS(BandGap_Nitride,workset)
{
  using panzer::index_t;

  // Note all energy-related fields saved in the FM are in the units of eV, 
  // and all other fields are scaled !  
  
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain temperature in [K]
      ScalarT lattT = latt_temp(cell,point) * T0;

      // compute band gap in [eV]
      if (isBinary)
        band_gap(cell,point) = binaryBandgap(lattT, materialName) ;   
      if (isTernary)
      {
        ScalarT mole_frac = molefrac(cell,point); 
        band_gap(cell,point) = ternaryBandgap(lattT, materialName, mole_frac) ;
      }   
    }
  }

}

//**********************************************************************
template<typename EvalT, typename Traits>
auto
BandGap_Nitride<EvalT, Traits>::binaryBandgap(const ScalarT& lattT, 
                                              const std::string& materialName)
  ->typename BandGap_Nitride<EvalT,Traits>::ScalarT
{
  using std::pow;
  ScalarT result;

  if (materialName == "GaN")
    result = 3.507 - (0.909e-3*pow(lattT,2)/(lattT + 830.0));
  else if (materialName == "AlN")
    result = 6.23 - (1.799e-3*pow(lattT,2)/(lattT + 1462.0));
  else if (materialName == "InN")
    result = 1.994 - (0.245e-3*pow(lattT,2)/(lattT + 624.0));
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid binary bandgap material: "
      << materialName << "!" << std::endl);

  return result;
}
//**********************************************************************
template<typename EvalT, typename Traits>
auto
BandGap_Nitride<EvalT, Traits>::ternaryBandgap(const ScalarT& lattT, 
                                              const std::string& materialName, 
                                              const ScalarT& mole_frac)
  ->typename BandGap_Nitride<EvalT,Traits>::ScalarT
{
  ScalarT eggan;
  ScalarT egaln;
  ScalarT eginn;
  ScalarT result;

  if (materialName == "AlGaN")
  {
    egaln = binaryBandgap(lattT, "AlN");
    eggan = binaryBandgap(lattT, "GaN");
    result = egaln*mole_frac + eggan*(1 - mole_frac) - 
                                   1.3*mole_frac*(1 - mole_frac);
  }

  else if (materialName == "InGaN")
  {
    eginn = binaryBandgap(lattT, "InN");
    eggan = binaryBandgap(lattT, "GaN");
    result = eginn*mole_frac + eggan*(1 - mole_frac) - 
                                    3.8*mole_frac*(1 - mole_frac);
  }
  else if (materialName == "AlInN")
  {
    egaln = binaryBandgap(lattT, "AlN");
    eginn = binaryBandgap(lattT, "InN");
    result = egaln*mole_frac + eginn*(1 - mole_frac);
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid ternary bandgap material: "
      << materialName << "!" << std::endl);
  
  return result;
}

//**********************************************************************

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
BandGap_Nitride<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Material Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Bandgap ParameterList", false, "");
  p->sublist("Bandgap ParameterList").set<std::string>("Value", "Nitride", "Specify nitride band gap");  

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);
  
  return p;
}

//**********************************************************************

}

#endif
