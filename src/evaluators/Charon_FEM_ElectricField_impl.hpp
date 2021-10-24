
#ifndef CHARON_FEM_ELECTRICFIELD_IMPL_HPP
#define CHARON_FEM_ELECTRICFIELD_IMPL_HPP

#include "Kokkos_ViewFactory.hpp"

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"
#include "Charon_Physical_Constants.hpp"


/*
The effective electric field Fn\p,eff = grad(negEffPot_n\p) (scaled), where
negEffPot_n\p = [Ei -\+ 0.5*dEg -\+ kbT*log(nie/nie0)] /V0, (- for n, + for p)
              = [Ei -\+ 0.5*dEg -\+ 0.5*kbT*log(gamma_n*gamma_p)] /V0.
nie = nie0*sqrt(gamma_n*gamma_p), nie0 = equilibrium effective intrinsic conc.,
and gamma_n\p = electron\hole degeneracy factor.
The calculation is valid for Boltzmann and Fermi-Dirac statistics, and
for BGN = On and Off cases. dEg is equal to 0 if BGN = Off.

Reference: D. Schroeder, T. Ostermann and O. Kalz, "Comparison of
transport models for the simulation of degenerate semiconductors,"
Semicond. Sci. Technol.9 (1994) 364-369.

This evaluator basically implements Eqn.(30)-(31) in the paper except that we
evaluate the scaled version. Note that Eqn.(30)-(31) are applicable to general
devices (homo- and hetero-geneous, the latter has constant mole fraction)
for Boltzmann and Fermi-Dirac statistics.

As demonstrated in the paper, the Boltzmann statistics works well in most
situations and even in certain degenerate regions, provided they are charge
neutral and the BGN effect is included through certain BGN models.
*/


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
FEM_ElectricField<EvalT, Traits>::
FEM_ElectricField(
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

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<panzer::IntegrationRule> ir =
    p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;
  num_points = vector->dimension(1);
  num_dims = vector->dimension(2);

  // BASIS
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();
  num_basis = basis_scalar->dimension(1);

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    efield = MDField<ScalarT,Cell,Point,Dim>(n.field.elec_efield,vector);
    grad_qfp = MDField<ScalarT,Cell,Point,Dim>(n.field.elec_grad_qfp,vector);
    grad_density = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.edensity,vector);
    density = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
    sign = -1.0;
  }
  else if (carrType == "Hole")
  {
    efield = MDField<ScalarT,Cell,Point,Dim>(n.field.hole_efield,vector);
    grad_qfp = MDField<ScalarT,Cell,Point,Dim>(n.field.hole_grad_qfp,vector);
    grad_density = MDField<const ScalarT,Cell,Point,Dim>(n.grad_dof.hdensity,vector);
    density = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
    sign = 1.0;
  }

  // Carrier-independent fields
  potential = MDField<const ScalarT,Cell,Point>(n.dof.phi,basis_scalar);
  intrinfermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi,basis_scalar);
  bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap,basis_scalar);
  effbandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap,basis_scalar);

  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, basis_scalar);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  V0 = scaleParams->scale_params.V0;
  T0 = scaleParams->scale_params.T0;

  // Degeneracy factors
  elec_degfactor = MDField<const ScalarT,Cell,Point>(n.field.elec_deg_factor,basis_scalar);
  hole_degfactor = MDField<const ScalarT,Cell,Point>(n.field.hole_deg_factor,basis_scalar);

  // Evaluated field
  this->addEvaluatedField(efield);
  this->addEvaluatedField(grad_qfp);

  // Dependent fields
  this->addDependentField(grad_density);
  this->addDependentField(density);
  this->addDependentField(potential);
  this->addDependentField(intrinfermi);
  this->addDependentField(bandgap);
  this->addDependentField(effbandgap);

  this->addDependentField(latt_temp);

  this->addDependentField(elec_degfactor);
  this->addDependentField(hole_degfactor);

  std::string name = "FEM_ElectricField";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_ElectricField<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  negEffPot = Kokkos::createDynRankView(potential.get_static_view(),"negEffPot",intrinfermi.dimension(0), num_basis);

  negPot = Kokkos::createDynRankView(potential.get_static_view(),"negPot",intrinfermi.dimension(0), num_basis);
  gradNegPot = Kokkos::createDynRankView(potential.get_static_view(),"gradNegPot",efield.dimension(0), num_points, num_dims);

}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
FEM_ElectricField<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain kb
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kbBoltz = cpc.kb;      // Boltzmann constant in [eV/K]

  // zero out arrays (important! otherwise, the Intrepid2::FunctionSpaceTools<PHX::exec_space>::evaluate
  // appears to sum in the destination field, leading to simulation divergence !)
  efield.deep_copy(ScalarT(0.0));
  for (std::size_t cell = 0; cell < gradNegPot.extent(0); ++cell)
    for (std::size_t pt = 0; pt < gradNegPot.extent(1); ++pt)
      for (std::size_t dim = 0; dim < gradNegPot.extent(2); ++dim)
        gradNegPot(cell,pt,dim) = ScalarT(0.0);

  // compute the effective potential at BASIS
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int basis = 0; basis < num_basis; ++basis)
    {
      // Lattice temperature in [K]
      //ScalarT lattT = latt_temp(cell,basis)*T0(0,0);
      ScalarT lattT = latt_temp(cell,basis)*T0;
      ScalarT kbT = kbBoltz*lattT;  // [eV]

      // obtain the band gap narrowing
      const ScalarT& Eg = bandgap(cell,basis);       // [eV]
      const ScalarT& effEg = effbandgap(cell,basis); // [eV]
      ScalarT deltaEg = Eg - effEg;

      // obtain the intrinsic fermi level
      const ScalarT& Ei = intrinfermi(cell,basis);   // [eV]

      // obtain the degeneracy factor
      const ScalarT& gamma_n = elec_degfactor(cell,basis); // unitless
      const ScalarT& gamma_p = hole_degfactor(cell,basis);

      // compute the effective potential
      negEffPot(cell,basis) = (Ei+sign*0.5*deltaEg+sign*0.5*kbT*log(gamma_n*gamma_p)) /V0;  // scaled

      const ScalarT& phi = potential(cell,basis);
      negPot(cell,basis) = -phi;

      //if (abs(negPot(cell,basis)-negEffPot(cell,basis)) > 1e-10)
      //  std::cout << "cell=" << cell <<", basis=" << basis << ", -phi= " <<
      //    negPot(cell,basis) << ", negEffPot=" << negEffPot(cell,basis) << std::endl;

    }
  }

  // compute the effective electric field at IPs (gradient of negEffPot)
  if(workset.num_cells>0)
  {
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::evaluate(efield.get_view(),negEffPot,(workset.bases[basis_index])->grad_basis.get_view());
    Intrepid2::FunctionSpaceTools<PHX::exec_space>::evaluate(gradNegPot,negPot,(workset.bases[basis_index])->grad_basis.get_view());
  }

  // compute the gradient of quasifermi potential at IPs
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t point = 0; point < num_points; ++point)
    {
      const ScalarT& dens = density(cell,point);

      for (std::size_t dim = 0; dim < num_dims; ++ dim)
      {
        const ScalarT& grad_dens = grad_density(cell,point,dim);
        grad_qfp(cell,point,dim) = sign * grad_dens/dens - efield(cell,point,dim);
      }
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
FEM_ElectricField<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;

}

}

#endif

