
#ifndef CHARON_GATETUNNELINGCURRENTDENSITY_IMPL_HPP
#define CHARON_GATETUNNELINGCURRENTDENSITY_IMPL_HPP

#include <cmath>
#include <map>

#include "Teuchos_TestForException.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_FieldLibrary.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"

#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include <iostream>


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
GateTunnelingCurrentDensity<EvalT, Traits>::
GateTunnelingCurrentDensity(const Teuchos::ParameterList& p)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using panzer::IntegrationRule;
  using std::string;
  using std::vector; 
  using Teuchos::rcp;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // get output data layout
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  num_ips = ip_scalar->dimension(1);

  // obtain basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis_scalar->dimension(1);
  basis_name = basis->name();
  
  // get charon::Names
  const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  //J0 = scaleParams->scale_params.J0;
  //T0 = scaleParams->scale_params.J0;

  // sideset ID
  sidesetID  = p.get<string>("Sideset ID");
  
  // gate sideset ID
  gate_sidesetID = p.get<string>("Sideset ID");

  // gate distance
  gate_dist = p.get<double>("Gate Distance");

  // geometry
  blockID = p.get<string>("Block ID");

  // get flux names
  tunnel_current_name = p.get<string>("Tunneling Current Density");

  // data layouts
  RCP<DataLayout> output_scalar = ip_scalar;
  RCP<DataLayout> input_scalar = basis_scalar; 
  
  // evaluated field
  tunnel_current = MDField<ScalarT,Cell,Point>(tunnel_current_name, output_scalar);
  this->addEvaluatedField(tunnel_current);
  
  // dependent fields
  latt_temp = MDField<ScalarT,Cell,Point>(names.field.latt_temp, input_scalar);
  pot = MDField<ScalarT,Cell,Point>(names.dof.phi, input_scalar);
  if (tunnel_current_name == "eGateTunnelingCurrentDensity") 
    carr_dens = MDField<ScalarT,Cell,Point>(names.dof.edensity, input_scalar);
  else 
    carr_dens = MDField<ScalarT,Cell,Point>(names.dof.hdensity, input_scalar);
  this->addDependentField(carr_dens);

  string n = "Gate Tunneling Current Density";
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void GateTunnelingCurrentDensity<EvalT, Traits>::postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */) {
  using std::vector;
  using Teuchos::RCP;

  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void GateTunnelingCurrentDensity<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset) {
  using panzer::index_t;

  //const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  //double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  tunnel_current.deep_copy(ScalarT(0.0));

  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (index_t point = 0; point < num_ips; ++point) {
        tunnel_current(cell,point) = 0.0;
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
GateTunnelingCurrentDensity<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector; 

  RCP<ParameterList> p = rcp(new ParameterList);

  RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  
  RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  RCP<const charon::Names> n;
  p->set("Names", n);

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set<string>("Tunneling Current Density", "?", "Tunneling Current Density");
  
  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));

  p->set<std::string>("Sideset ID","?");

  p->set<std::string>("Gate Sideset ID","?");

  p->set<double>("Gate Distance", 0.0, "Distance to tunneling gate in cm");

  p->set<std::string>("Block ID","?");

  return p;
}

}

#endif

