
#ifndef CHARON_NEUMANNBC_THERMALCONTACT_IMPL_HPP
#define CHARON_NEUMANNBC_THERMALCONTACT_IMPL_HPP

#include <cmath>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_String_Utilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"

/*
The Neumann Thermal Contact sets kappa \grad_T \cdot \norm either to a constant
specified by Power in [W/cm^2], or to (Text-T)/Rth specified by Surface Resistance
with Rth being in the units of [K.cm^2/W], or to Gth(Text-T) specified by Surface
Conductance with Gth in [W/(K.cm^2)]. Only one of the three is allowed to be used.
When either Surface Resistance or Conductance is given, Temperature must also be given.
The specification in the input xml takes the following form:

<ParameterList>
    <Parameter name="Type" type="string" value="Neumann"/>
    <Parameter name="Sideset ID" type="string" value="anode"/>
    <Parameter name="Element Block ID" type="string" value="silicon"/>
    <Parameter name="Equation Set Name" type="string" value="Lattice Temperature"/>
    <Parameter name="Strategy" type="string" value="Thermal Contact"/>
    <ParameterList name="Data">
        <Parameter name="Power" type="double" value="1e5"/>
        // <Parameter name="Surface Resistance" type="double" value="0.1"/>
        // <Parameter name="Surface Conductance" type="double" value="10"/>
        // <Parameter name="Temperature" type="double" value="300"/>
    </ParameterList>
</ParameterList>

*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
NeumannBC_ThermalContact<EvalT, Traits>::
NeumannBC_ThermalContact(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // get data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // get names and values
  string flux_name = p.get<string>("Flux Name");
  string dof_name = p.get<string>("DOF Name");
  paramName = p.get<string>("Parameter Name");
  value = p.get<double>("Value");
  temp = p.get<double>("Temperature");

  // evaluated field
  heat_flux = MDField<ScalarT,Cell,Point>(flux_name, scalar);
  this->addEvaluatedField(heat_flux);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  H0 = scaleParams->scale_params.H0;
  X0 = scaleParams->scale_params.X0;

  // dependent fields
  if ( (paramName == "Surface Resistance") || (paramName == "Surface Conductance") )
  {
    latt_temp = MDField<const ScalarT,Cell,Point>(dof_name, scalar);
    this->addDependentField(latt_temp);
    T0 = scaleParams->scale_params.T0;
  }

  std::string n = "NeumannBC Thermal Contact";
  this->setName(n);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
NeumannBC_ThermalContact<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  ScalarT scaling = H0*X0;

  if (paramName == "Power")
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (int point = 0; point < num_points; ++point)
        heat_flux(cell,point) = value / scaling;
  }

  else if (paramName == "Surface Resistance")
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        // obtain temperature in [K]
        ScalarT latT = latt_temp(cell,point) * T0;

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latT) <= 0.0)  latT = 300.0;

        // compute (Text-T)/Rth in [W/cm^2]
        ScalarT power = (temp - latT) / value;

        heat_flux(cell,point) = power / scaling;
      }
    }
  }

  else if (paramName == "Surface Conductance")
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int point = 0; point < num_points; ++point)
      {
        // obtain temperature in [K]
        ScalarT latT = latt_temp(cell,point) * T0;

        // lattT should be always > 0, but it could become <= 0 due to numerical errors
        // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
        if (Sacado::ScalarValue<ScalarT>::eval(latT) <= 0.0)  latT = 300.0;

        // compute (Text-T)*Gth in [W/cm^2]
        ScalarT power = (temp - latT) * value;

        heat_flux(cell,point) = power / scaling;
      }
    }
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
       << "Error: Wrong parameter name is specified ! Must be either Power, or Surface Resistance "
       << " or Surface Conductance. But you specified " << paramName << " !");
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
NeumannBC_ThermalContact<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<std::string>("Flux Name", "???");
  p->set<std::string>("DOF Name", "???");
  p->set<std::string>("Parameter Name", "???");

  p->set<double>("Value", 0.0);
  p->set<double>("Temperature", 0.0);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif

