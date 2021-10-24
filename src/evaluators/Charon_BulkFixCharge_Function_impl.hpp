
#ifndef CHARON_BULKFIXCHARGE_FUNCTION_IMPL_HPP
#define CHARON_BULKFIXCHARGE_FUNCTION_IMPL_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Panzer_GlobalData.hpp"
#include "Charon_Names.hpp"

/** 
 * @brief "Fixed Charge" specification is similar to that of Doping, except that 
 * Fixed Charge currently supports only Uniform profiles.
 * An example of specifying Fixed Charge is given below:
 *
 * <ParameterList name="Fixed Charge">
 *    <Parameter name="Value" type="string" value="Function"/>
 *    <ParameterList name="Function 1">
 *       <Parameter name="Function Type" type="string" value="Uniform"/>
 *       <Parameter name="Charge Density" type="double" value="1e18"/>
 *    </ParameterList>
 *    <ParameterList name="Function 2">
 *       <Parameter name="Function Type" type="string" value="Uniform"/>
 *       <Parameter name="Charge Density" type="double" value="-5e17"/>
 *       <Parameter name="Xmin" type="double" value="0.6"/>
 *       <Parameter name="Xmax" type="double" value="1.0"/>
 *       <Parameter name="Ymin" type="double" value="0.0"/>
 *       <Parameter name="Ymax" type="double" value="0.1"/>
 *    </ParameterList>
 * </ParameterList>
 */

namespace charon {

void uniformBulkFixQParams::parseUniform (const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;

  value = 0.0;
  varyingFixCharge = false;
  if(plist.isParameter("Charge Density"))
    {
      value = plist.get<double>("Charge Density");  // in unit of cm^{-3}
    }
  if(plist.isParameter("Varying Charge Density"))
    {
      varyingFixCharge = true;
    }

  xmin = -1e100, ymin = -1e100, zmin = -1e100;
  xmax =  1e100, ymax =  1e100, zmax =  1e100;

  if (plist.isParameter("Xmin"))  xmin = plist.get<double>("Xmin");
  if (plist.isParameter("Xmax"))  xmax = plist.get<double>("Xmax");
  if (plist.isParameter("Ymin"))  ymin = plist.get<double>("Ymin");
  if (plist.isParameter("Ymax"))  ymax = plist.get<double>("Ymax");
  if (plist.isParameter("Zmin"))  zmin = plist.get<double>("Zmin");
  if (plist.isParameter("Zmax"))  zmax = plist.get<double>("Zmax");
  
  // std::cout << "Xmin = " << xmin << ", Xmax = " << xmax << std::endl; 
}


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
BulkFixCharge_Function<EvalT, Traits>::
BulkFixCharge_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using Teuchos::rcp;

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  RCP<DataLayout> ip_vector = ir->dl_vector;
  int_rule_degree = ir->cubature_degree;
  num_ips = ip_vector->dimension(1);
  num_dims = ip_vector->dimension(2);

  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();
  num_basis = basis_scalar->dimension(1);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;

  // bulk fixed charge parameterlist
  bulkFixQParamList = p.sublist("Bulk FixCharge ParameterList");


  // evaluated fields
  fixed_charge = MDField<ScalarT,Cell,Point>(n.field.fixed_charge,ip_scalar);
  fixed_charge_basis = MDField<ScalarT,Cell,Point>(n.field.fixed_charge,basis_scalar);
  this->addEvaluatedField(fixed_charge);
  this->addEvaluatedField(fixed_charge_basis);

  std::string name = "BulkFixCharge_Function";
  this->setName(name);

  for (ParameterList::ConstIterator model_it = bulkFixQParamList.begin();
       model_it != bulkFixQParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      if (funcType == "Uniform")
      {
        uniformBulkFixQParams umfp_;
        umfp_.parseUniform(funcParamList);
        umfp_vec.push_back(umfp_);
	varyingFixCharge = false;
	if(umfp_.varyingFixCharge)
	  {
	    varyingFixCharge = true;
	  }
      }
      
      // to be done for other function types

    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop


  if(varyingFixCharge)
    {
      // varying bulk fixed charge parameters
      user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
      user_value->setRealValue(1);
      //      bulkFixQParamList.set<std::string>("Varying Charge Density","Parameter");
      //bulkFixQParamList.set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
      //						     Teuchos::rcp(new panzer::ParamLib));
      //bulkFixQParamList.set<Teuchos::RCP<panzer::ParamLib> >("ParamLib", this->getGlobalData()->pl);  

      user_value =
	panzer::createAndRegisterScalarParameter<EvalT>(
							std::string("Varying Charge Density"),
							*bulkFixQParamList.get<RCP<panzer::ParamLib> >("ParamLib"));
    }

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BulkFixCharge_Function<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
BulkFixCharge_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // bulk fixed charge density at IPs
    for (int ip = 0; ip < num_ips; ++ip)
    {
      double x = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,0);
      double y = 0.0, z = 0.0;
      if (num_dims == 2)
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
      if (num_dims == 3)
      {
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
        z = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,2);
      }

      // std::cout << "IP: cell=" << cell << ", ip=" << ip << ", x=" << x << ", y=" << y << ", z=" << std::endl; 

      // evaluate the bulk fixed charge density
      double value = evaluateBulkFixCharge(x,y,z);  // in cm^{-3}
      fixed_charge(cell,ip) = value / C0;  // scale the charge density
    }

    // bulk fixed charge density at basis points
    for (int basis = 0; basis < num_basis; ++basis)
    {
      double x = (workset.bases[basis_index])->basis_coordinates(cell,basis,0);
      double y = 0.0, z = 0.0;
      if (num_dims == 2)
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
      if (num_dims == 3)
      {
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
        z = (workset.bases[basis_index])->basis_coordinates(cell,basis,2);
      }

      // std::cout << "Basis: cell=" << cell << ", basis=" << basis << ", x=" << x << ", y=" << y << ", z=" << std::endl; 
      
      double value = evaluateBulkFixCharge(x,y,z);
      fixed_charge_basis(cell,basis) = value / C0;
    }

  } // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateBulkFixCharge()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double BulkFixCharge_Function<EvalT, Traits>::evaluateBulkFixCharge
  (const double& x, const double& y, const double& z)
{
  using std::string;
  using Teuchos::ParameterList;

  double qDensity = 0.0;
  double tempVal = 0.0;

  for (std::size_t i = 0; i < umfp_vec.size(); ++i)
  {
    tempVal = evalUniformBulkFixQ(x,y,z,umfp_vec[i]);
    qDensity += tempVal;
  }
  
  // to be done for other spatially dependent fixed charges 
  
  return qDensity;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalUniformBulkFixQ()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double BulkFixCharge_Function<EvalT, Traits>::evalUniformBulkFixQ
  (const double& x, const double& y, const double& z, const uniformBulkFixQParams& umfp)
{
  using std::string;
  using Teuchos::ParameterList;

  double qDensity = 0.0;
  double varyingCharge;
  const double val = umfp.value;
  if(varyingFixCharge)
    {
      varyingCharge = user_value->getRealValue();
    }
  else
    varyingCharge = val;
  const double xmin = umfp.xmin;
  const double ymin = umfp.ymin;
  const double zmin = umfp.zmin;
  const double xmax = umfp.xmax;
  const double ymax = umfp.ymax;
  const double zmax = umfp.zmax;

  if ( (x >= xmin) && (x <= xmax) && (y >= ymin) && (y <= ymax) && (z >= zmin) && (z <= zmax) )
    qDensity = varyingCharge;

  // return 0 if (x,y,z) is outside the box region

  return qDensity;
}

} // namespace charon

#endif
