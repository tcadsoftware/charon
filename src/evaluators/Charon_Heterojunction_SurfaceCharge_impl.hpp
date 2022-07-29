
#ifndef CHARON_HETEROJUNCTION_SURFACECHARGE_IMPL_HPP
#define CHARON_HETEROJUNCTION_SURFACECHARGE_IMPL_HPP

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
Heterojunction_SurfaceCharge<EvalT, Traits>::
Heterojunction_SurfaceCharge(const Teuchos::ParameterList& p)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using std::string;
  using std::vector; 
  using Teuchos::rcp;

  auto valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  // get output data layout
  auto output_dl = p.get< RCP<DataLayout> >("Output Data Layout");
  num_ips = output_dl->dimension(1);

  // obtain basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> input_dl = basis->functional;
  basis_name = basis->name();
  num_nodes = input_dl->dimension(1);
  
  // get charon::Names
  //const charon::Names& names = *(p.get< RCP<const charon::Names> >("Names"));

  // scaling parameters
  auto scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;
  X0 = scaleParams->scale_params.X0;
  
  // get flux names
  fluxSurfCharge = p.get<string>("Flux Surface Charge");

  // read in user-specified charge
  fixedCharge = rcp(new panzer::ScalarParameterEntry<EvalT>);
  fixedCharge->setRealValue(0.0);
  fixedCharge->setRealValue(p.get<double>("Fixed Charge"));  // in unit of cm^{-2}
  
  surface_charge = MDField<ScalarT,Cell,Point>(fluxSurfCharge, output_dl);
  this->addEvaluatedField(surface_charge);
 
  string n = "Heterojunction Surface Charge";
  this->setName(n);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Heterojunction_SurfaceCharge<EvalT, Traits>::postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Heterojunction_SurfaceCharge<EvalT, Traits>::evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  double fixedChargeLocal = fixedCharge->getRealValue();
  double scaling4Charge = C0*X0;
  double scaledFixCharge = fixedChargeLocal / scaling4Charge; 

  //const charon::PhysicalConstants & phyConst = charon::PhysicalConstants::Instance();
  //double kb = phyConst.kb;  // Boltzmann constant in [eV/K]
  
  // zero out the arrays
  surface_charge.deep_copy(ScalarT(0.0));

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
    for (int point = 0; point < num_ips; ++point) {
      surface_charge(cell,point) += scaledFixCharge; 
    }
  
}




///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Heterojunction_SurfaceCharge<EvalT, Traits>::getValidParameters() const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::string;
  using std::vector; 

  RCP<ParameterList> p = rcp(new ParameterList);

  RCP<PHX::DataLayout> dl;
  p->set("Output Data Layout", dl);
  
  RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  p->set<string>("Flux Surface Charge", "?", "Name for the computed scaled surface charge");

  p->set<Teuchos::RCP<panzer::ParamLib> >("ParamLib",
       Teuchos::rcp(new panzer::ParamLib));

  p->set<double>("Fixed Charge", 0.0, "Fixed Charge in unit of cm^(-2)"); 

  return p;
}

}

#endif

