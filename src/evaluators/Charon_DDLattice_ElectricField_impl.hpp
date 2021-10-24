
#ifndef CHARON_DDLATTICE_ELECTRICFIELD_IMPL_HPP
#define CHARON_DDLATTICE_ELECTRICFIELD_IMPL_HPP

#include "Kokkos_ViewFactory.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"


/*
Evaluate the electric field at IPs for the DD+Ion+Lattice FEM formulations. The
expression for the field depends on the model of choice. Currently,
Fn = -\grad_phi - 0.5 * \grad(bgn/q) / V0,
Fp = -\grad_phi + 0.5 * \grad(bgn/q) / V0,
Fi = -\grad_phi. (i for ion)
Fn and Fp include the BGN contribution. They can be expanded to include other
more complex models. All quantities are scaled here.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLattice_ElectricField<EvalT, Traits>::
DDLattice_ElectricField(
  const Teuchos::ParameterList& p)
{
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<panzer::IntegrationRule> ir =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_cells = vector->dimension(0);
  num_ips = vector->dimension(1);
  num_dims = vector->dimension(2);

  // BASIS
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();
  num_basis = basis_scalar->dimension(1);

  carrType = p.get<string>("Carrier Type");
  modelName = p.get<string>("Electric Field Model");
  haveBGN = p.get<bool>("Band Gap Narrowing");

  // Evaluated field
  if (carrType == "Electron")
  {
    electric_field = MDField<ScalarT,Cell,IP,Dim>(n.field.elec_efield,vector);
    sign = -1.0;
  }
  else if (carrType == "Hole")
  {
    electric_field = MDField<ScalarT,Cell,IP,Dim>(n.field.hole_efield,vector);
    sign = 1.0;
  }
  else if (carrType == "Ion")
  {
    electric_field = MDField<ScalarT,Cell,IP,Dim>(n.field.ion_efield,vector);
    sign = 1.0;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: Invalid Carrier Type !"
        <<" Must be Electron, Hole, or Ion. You entered " << carrType << " ! \n");

  this->addEvaluatedField(electric_field);

  // Dependent fields
  if (modelName == "Potential Gradient")
  {
    grad_potential = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.phi,vector);
    this->addDependentField(grad_potential);
  }

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");

  if (haveBGN)
  {
    bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap,basis_scalar);
    eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,basis_scalar);
    this->addDependentField(bandgap);
    this->addDependentField(eff_bandgap);

    V0 = scaleParams->scale_params.V0;
  }

  std::string name = "DDLattice_ElectricField";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLattice_ElectricField<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);
  dEg = Kokkos::createDynRankView(electric_field.get_static_view(),"dEg",num_cells, num_basis);
  grad_dEg = Kokkos::createDynRankView(electric_field.get_static_view(),"grad_dEg",num_cells, num_ips, num_dims);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLattice_ElectricField<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

 if (modelName == "Potential Gradient")
 {
  if (haveBGN)  // BGN = ON
  {
   if (carrType == "Ion")  // NO BGN for Ion
   {
     for (index_t cell = 0; cell < workset.num_cells; ++cell)
       for (std::size_t ip = 0; ip < num_ips; ++ip)
         for (std::size_t dim = 0; dim < num_dims; ++dim)
           electric_field(cell,ip,dim) = -grad_potential(cell,ip,dim);
   }

   else  // for Electron and Hole
   {
     // zero out arrays
     for (std::size_t cell = 0; cell < grad_dEg.extent(0); ++cell)
       for (std::size_t ip = 0; ip < grad_dEg.extent(1); ++ip)
         for (std::size_t dim = 0; dim < grad_dEg.extent(2); ++dim)
           grad_dEg(cell,ip,dim) = ScalarT(0.0);

     // compute dEg at BASIS
     for (index_t cell = 0; cell < workset.num_cells; ++cell)
       for (std::size_t basis = 0; basis < num_basis; ++basis)
         dEg(cell,basis) = (bandgap(cell,basis) - eff_bandgap(cell,basis)) / V0;

     // compute grad_dEg at IPs
     if(workset.num_cells > 0)
       Intrepid2::FunctionSpaceTools<PHX::exec_space>::evaluate(grad_dEg,dEg,(workset.bases[basis_index])->grad_basis.get_view());

     // compute electric field at IPs
     for (index_t cell = 0; cell < workset.num_cells; ++cell)
       for (std::size_t ip = 0; ip < num_ips; ++ip)
         for (std::size_t dim = 0; dim < num_dims; ++dim)
           electric_field(cell,ip,dim) = -grad_potential(cell,ip,dim) + sign*0.5*grad_dEg(cell,ip,dim);
   }

  }  // end of if (haveBGN)

  else  // BGN = OFF
  {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t ip = 0; ip < num_ips; ++ip)
        for (std::size_t dim = 0; dim < num_dims; ++dim)
          electric_field(cell,ip,dim) = -grad_potential(cell,ip,dim);
  }
 }  // end of "Potential Gradient"

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DDLattice_ElectricField<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  p->set("Carrier Type", "?");
  p->set("Electric Field Model", "Potential Gradient");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set("Band Gap Narrowing", false);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;

}

}

#endif

