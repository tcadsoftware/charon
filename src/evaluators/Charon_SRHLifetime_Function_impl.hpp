
#ifndef CHARON_SRHLIFETIME_FUNCTION_IMPL_HPP
#define CHARON_SRHLIFETIME_FUNCTION_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Names.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
SRHLifetime_Function<EvalT, Traits>::
SRHLifetime_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using Teuchos::ParameterList;

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Input from closure model
  const string& carrType = p.get<string>("Carrier Type");
  const string& matName = p.get<string>("Material Name");

  // Lifetime parameterlist
  const ParameterList& ltParamList = p.sublist("Lifetime ParameterList");

  // Obtain Material_Properties instance
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Obtain Tau0
  if (ltParamList.isParameter("Tau0"))
    tau0 = ltParamList.get<double>("Tau0");  // user specified
  else  // otherwise, retrieve from Material_Properties
  {
    if (carrType == "Electron")
      tau0 = matProperty.getPropertyValue(matName, "Electron Lifetime");
    else if (carrType == "Hole")
      tau0 = matProperty.getPropertyValue(matName, "Hole Lifetime");
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
        << "Invalid Carrier Type ! Must be either Electron or Hole !");
  }

  isConcDep = false;
  isTempDep = false;
  isExpTempDep = false;

  // Loop over lifetime functions
  for (ParameterList::ConstIterator model_it = ltParamList.begin();
       model_it !=ltParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");

      // concentration dependence
      if (funcType == "ConcDep")
      {
        isConcDep = true;
        if (funcParamList.isParameter("Nsrh"))
          nsrh = funcParamList.get<double>("Nsrh");
        else
        {
          if (carrType == "Electron")
            nsrh = matProperty.getPropertyValue(matName, "Electron SRH Conc");
          else if (carrType == "Hole")
            nsrh = matProperty.getPropertyValue(matName, "Hole SRH Conc");
        }
      }

      // power law temperature dependence
      else if (funcType == "TempDep")
      {
        isTempDep = true;
        if (funcParamList.isParameter("TPowerLaw"))
          Tpowlaw = funcParamList.get<double>("TPowerLaw");
        else
        {
          if (carrType == "Electron")
            Tpowlaw = matProperty.getPropertyValue(matName, "Electron SRH TPowerLaw");
          else if (carrType == "Hole")
            Tpowlaw = matProperty.getPropertyValue(matName, "Hole SRH TPowerLaw");
        }
      }

      // exponential temperature dependence
      else if (funcType == "ExpTempDep")
      {
        isExpTempDep = true;
        if (funcParamList.isParameter("TExponential"))
          Texp = funcParamList.get<double>("TExponential");
        else
        {
          if (carrType == "Electron")
            Texp = matProperty.getPropertyValue(matName, "Electron SRH TExponential");
          else if (carrType == "Hole")
            Texp = matProperty.getPropertyValue(matName, "Hole SRH TExponential");
        }
      }

      // field dependence
      else if (funcType == "FieldDep")
      {
        // TODO
      }

      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Invalid SRH lifetime Function Type! Must be either ConcDep, "
          << "TempDep, ExpTempDep, or FieldDep !" << std::endl );

    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop

  // TempDep and ExpTempDep cannot be specified at the same time
  if (isTempDep && isExpTempDep)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "TempDep and ExpTempDep CANNOT be specified together !");

  // Carrier-dependent field
  if (carrType == "Electron")
    lifetime = MDField<ScalarT,Cell,Point>(n.field.elec_lifetime,scalar);
  else if (carrType == "Hole")
    lifetime = MDField<ScalarT,Cell,Point>(n.field.hole_lifetime,scalar);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");

  // Evaluated field
  this->addEvaluatedField(lifetime);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;

  // Add dependent doping fields when isConcDep = true
  if (isConcDep)
  {
    acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor_raw,scalar);
    donor = MDField<const ScalarT,Cell,Point>(n.field.donor_raw,scalar);
    C0 = scaleParams->scale_params.C0;
    this->addDependentField(acceptor);
    this->addDependentField(donor);
  }

  // Add dependent lattice temp. when temp. dep. is specified
  if (isTempDep || isExpTempDep)
  {
    latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
    T0 = scaleParams->scale_params.T0;
    this->addDependentField(latt_temp);
  }

  std::string name = "SRHLifetime_Function";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
SRHLifetime_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      ScalarT ltValue = tau0;

      if (isConcDep)
      {
        const ScalarT& Na = acceptor(cell,point);
        const ScalarT& Nd = donor(cell,point);
        ScalarT ntot = (Na + Nd) * C0;
        ltValue *= 1.0 / (1.0 + ntot/nsrh );
      }

      if (isTempDep)
      {
        ScalarT lattT = latt_temp(cell,point)*T0;
        ltValue *= pow(lattT/300.0, Tpowlaw);
      }

      if (isExpTempDep)
      {
        ScalarT lattT = latt_temp(cell,point)*T0;
        ltValue *= exp( Texp * (lattT/300.0 - 1.0) );
      }

      lifetime(cell,point) = ltValue / t0;

    }
  }
}

}

#endif
