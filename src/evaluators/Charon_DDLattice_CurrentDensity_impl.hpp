
#ifndef CHARON_DDLATTICE_CURRENTDENSITY_IMPL_HPP
#define CHARON_DDLATTICE_CURRENTDENSITY_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Charon_Names.hpp"


/*
Evaluate the current density at IPs for the DDLattice, DDIon, and DDIonLattice
FEM formulations, where
Jn = n*mun*Fn + Dn*\grad_n + n*mun*\grad_T,  (n for electrons)
Jp = p*mup*Fp - Dp*\grad_p - p*mup*\grad_T,  (p for holes)
Ji = Ni*vi - Dichem*\grad_Ni - mui*T*Si*Ni*\grad_T
   = Ni*vi - Dichem*\grad_Ni - DTi*Ni*\grad_T, (i for ions or vacancies)
with all quantities in scaled units.

Note Dichem = Di/(1-N/Nmax), evaluated by charon::DiffCoeff_IonDep,
while Di = mui*T, evaluated according to the RigidPointIon model.

When bTempGrad = false, turn off the \grad_T term. By default, bTempGrad = true.
*/


namespace charon {


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DDLattice_CurrentDensity<EvalT, Traits>::
DDLattice_CurrentDensity(
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
  num_ips = vector->dimension(1);
  num_dims = vector->dimension(2);

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Current density name
  string currentName = p.get<string>("Current Name");

  // Determine if want to include the \grad_T contribution
  bTempGrad = true;  // default
  if (p.isParameter("Temperature Gradient"))
    bTempGrad = p.get<bool>("Temperature Gradient");

  // Evaluated field
  current_density = MDField<ScalarT,Cell,IP,Dim>(currentName,vector);
  this->addEvaluatedField(current_density);

  // Carrier-dependent fields
  if (carrType == "Electron")
  {
    electric_field = MDField<const ScalarT,Cell,IP,Dim>(n.field.elec_efield,vector);
    grad_carr_dens = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.edensity,vector);
    carr_dens = MDField<const ScalarT,Cell,IP>(n.dof.edensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,IP>(n.field.elec_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.elec_mobility,scalar);
    sign = 1.0;

    this->addDependentField(electric_field);
    this->addDependentField(grad_carr_dens);
    this->addDependentField(carr_dens);
    this->addDependentField(diff_coeff);
    this->addDependentField(mobility);
  }

  else if (carrType == "Hole")
  {
    electric_field = MDField<const ScalarT,Cell,IP,Dim>(n.field.hole_efield,vector);
    grad_carr_dens = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.hdensity,vector);
    carr_dens = MDField<const ScalarT,Cell,IP>(n.dof.hdensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,IP>(n.field.hole_diff_coeff,scalar);
    mobility = MDField<const ScalarT,Cell,IP>(n.field.hole_mobility,scalar);
    sign = -1.0;

    this->addDependentField(electric_field);
    this->addDependentField(grad_carr_dens);
    this->addDependentField(carr_dens);
    this->addDependentField(diff_coeff);
    this->addDependentField(mobility);
  }

  else if (carrType == "Ion")
  {
    velocity = MDField<const ScalarT,Cell,IP,Dim>(n.field.ion_velocity,vector);
    grad_carr_dens = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.iondensity,vector);
    carr_dens = MDField<const ScalarT,Cell,IP>(n.dof.iondensity,scalar);
    diff_coeff = MDField<const ScalarT,Cell,IP>(n.field.ion_diff_coeff,scalar);

    this->addDependentField(velocity);
    this->addDependentField(grad_carr_dens);
    this->addDependentField(carr_dens);
    this->addDependentField(diff_coeff);

    if (bTempGrad)
    {
      thermodiff_coeff = MDField<const ScalarT,Cell,IP>(n.field.ion_thermodiff_coeff,scalar);
      this->addDependentField(thermodiff_coeff);
    }

    // output to Exodus
    ion_curr_dens_x = MDField<ScalarT,Cell,IP>(n.field.ion_curr_dens_x,scalar);
    ion_curr_dens_y = MDField<ScalarT,Cell,IP>(n.field.ion_curr_dens_y,scalar);
    this->addEvaluatedField(ion_curr_dens_x);
    this->addEvaluatedField(ion_curr_dens_y);
  }

  // Carrier-independent fields
  if (bTempGrad)
  {
    grad_latt_temp = MDField<const ScalarT,Cell,IP,Dim>(n.grad_dof.latt_temp,vector);
    this->addDependentField(grad_latt_temp);
  }

  std::string name = "DDLattice_CurrentDensity";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DDLattice_CurrentDensity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // ***************************************************************************
  // include the \grad_T contribution
  // ***************************************************************************

  if (bTempGrad)
  {
   if ( (carrType == "Electron") || (carrType == "Hole") )
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (std::size_t ip = 0; ip < num_ips; ++ip)
      {
        const ScalarT& carr = carr_dens(cell,ip);
        const ScalarT& diff = diff_coeff(cell,ip);
        const ScalarT& mob = mobility(cell,ip);

        for (std::size_t dim = 0; dim < num_dims; ++dim)
        {
          const ScalarT& field = electric_field(cell,ip,dim);
          const ScalarT& gradcarr = grad_carr_dens(cell,ip,dim);
          const ScalarT& gradtemp = grad_latt_temp(cell,ip,dim);
          current_density(cell,ip,dim) = carr*mob*field + sign*diff*gradcarr + sign*carr*mob*gradtemp;
        }
      }
    }
   }  // carrType = Electron or Hole

   else if (carrType == "Ion")
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (std::size_t ip = 0; ip < num_ips; ++ip)
      {
        const ScalarT& carr = carr_dens(cell,ip);
        const ScalarT& diff = diff_coeff(cell,ip);
        const ScalarT& thermodiff = thermodiff_coeff(cell,ip);

        for (std::size_t dim = 0; dim < num_dims; ++dim)
        {
          const ScalarT& vel = velocity(cell,ip,dim);
          const ScalarT& gradcarr = grad_carr_dens(cell,ip,dim);
          const ScalarT& gradtemp = grad_latt_temp(cell,ip,dim);
          current_density(cell,ip,dim) = carr*vel - diff*gradcarr - thermodiff*carr*gradtemp;
        }

        ion_curr_dens_x(cell,ip) = current_density(cell,ip,0); // for outputting
        ion_curr_dens_y(cell,ip) = current_density(cell,ip,1);
      }
    }
   }  // end of carrType = Ion

  }  // end of the bTempGrad = true block


  // ***************************************************************************
  // exclude the \grad_T contribution
  // ***************************************************************************

  else
  {
   if ( (carrType == "Electron") || (carrType == "Hole") )
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (std::size_t ip = 0; ip < num_ips; ++ip)
      {
        const ScalarT& carr = carr_dens(cell,ip);
        const ScalarT& diff = diff_coeff(cell,ip);
        const ScalarT& mob = mobility(cell,ip);

        for (std::size_t dim = 0; dim < num_dims; ++dim)
        {
          const ScalarT& field = electric_field(cell,ip,dim);
          const ScalarT& gradcarr = grad_carr_dens(cell,ip,dim);
          current_density(cell,ip,dim) = carr*mob*field + sign*diff*gradcarr;
        }
      }
    }
   }  // end of carrType = Electron or Hole

   else if (carrType == "Ion")
   {
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (std::size_t ip = 0; ip < num_ips; ++ip)
      {
        const ScalarT& carr = carr_dens(cell,ip);
        const ScalarT& diff = diff_coeff(cell,ip);

        for (std::size_t dim = 0; dim < num_dims; ++dim)
        {
          const ScalarT& vel = velocity(cell,ip,dim);
          const ScalarT& gradcarr = grad_carr_dens(cell,ip,dim);
          current_density(cell,ip,dim) = carr*vel - diff*gradcarr;
        }

        ion_curr_dens_x(cell,ip) = current_density(cell,ip,0); // for outputting
        ion_curr_dens_y(cell,ip) = current_density(cell,ip,1);
      }
    }
   }  // end of carrType = Ion

  }  // end of the bTempGrad = false block

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
DDLattice_CurrentDensity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Carrier Type", "?");
  p->set("Current Name", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  p->set<bool>("Temperature Gradient", true, "Turn on the temperature gradient contribution by default !");

  return p;

}

}

#endif

