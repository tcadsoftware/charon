
#ifndef CHARON_DEGENERACY_FACTOR_IMPL_HPP
#define CHARON_DEGENERACY_FACTOR_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

/*
Compute the Fermi-Dirac degeneracy factor, see the reference by D. Schroeder
et. al, "Comparison of transport models for the simulation of degenerate
semiconductors," Semicond. Sci. Technol. 9, 364 (1994).
degeneracy factor = 1 for Boltzmann statistics,
degeneracy factor < 1 and > 0 for Fermi-Dirac statistics.

Inverse of the Fermi-Dirac integral of 1/2 order is computed using the
approximate expression by N. G. Nilsson, "An Accurate Approximation of the
Generalized Einstein Relation for Degenerate Semiconductors,"
Phys. Stat. Sol. (a) 19, K75 (1973). This approximation has an error of less
than 0.6% over the entire argument range.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Degeneracy_Factor<EvalT, Traits>::
Degeneracy_Factor(
  const Teuchos::ParameterList& p)
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

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Determine if Fermi Dirac is turned on
  bUseFD = p.get<bool>("Fermi Dirac");

  // Get the FD Formula
  fdFormula = p.get<string>("FD Formula");

  // Evaluated fields
  elec_degfactor = MDField<ScalarT,Cell,Point>(n.field.elec_deg_factor,scalar);
  hole_degfactor = MDField<ScalarT,Cell,Point>(n.field.hole_deg_factor,scalar);

  this->addEvaluatedField(elec_degfactor);
  this->addEvaluatedField(hole_degfactor);

  if (bUseFD)  // turn on the Fermi Dirac statistics
  {
    // Dependent fields
    elec_density = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
    hole_density = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);
    elec_effdos = MDField<const ScalarT,Cell,Point>(n.field.elec_eff_dos,scalar);
    hole_effdos = MDField<const ScalarT,Cell,Point>(n.field.hole_eff_dos,scalar);

    this->addDependentField(elec_density);
    this->addDependentField(hole_density);
    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
  }

  std::string name = "Degeneracy_Factor";
  this->setName(name);

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral =
    Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Degeneracy_Factor<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

 // Compute the FD degeneracy factors when bUseFD=true and fdFormula="Schroeder"
 if ( (bUseFD) && (fdFormula == "Schroeder") )
 {
  // loop over the cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      // obtain the e/h density (scale)
      const ScalarT& n = elec_density(cell,point);
      const ScalarT& p = hole_density(cell,point);

      // obtain the effective density of states (scaled)
      const ScalarT& Nc = elec_effdos(cell,point);
      const ScalarT& Nv = hole_effdos(cell,point);

      // compute the ratio
      ScalarT n_over_Nc = n/Nc;  // note: both n and Nc are scaled, n/Nc could be < 0,
      ScalarT p_over_Nv = p/Nv;  // since n could be < 0 during simulation.

      // compute the degeneracy factor
      if (n_over_Nc <= 1e-4)
        elec_degfactor(cell,point) = 1.0; // use Boltzmann statistics
      else
      {
        // compute inverse of the Fermi-Dirac integral of 1/2 order
        // ScalarT eta = fIntOneHalf->Inverse_Nilsson(n_over_Nc);
        ScalarT eta = (*inverseFermiIntegral)(n_over_Nc);
        elec_degfactor(cell,point) = n_over_Nc*exp(-eta);
      }

      if (p_over_Nv <= 1e-4)
        hole_degfactor(cell,point) = 1.0;
      else
      {
        // ScalarT nv = pow(sqrt(M_PI)*p_over_Nv*3.0/4.0, 2.0/3.0);
        // ScalarT nvterm = nv/(1.0+pow(0.24+1.08*nv,-2.0));
        // ScalarT eta = log(p_over_Nv)/(1.0-pow(p_over_Nv,2.0)) + nvterm;

        // ScalarT eta = fIntOneHalf->Inverse_Nilsson(p_over_Nv);
        ScalarT eta = (*inverseFermiIntegral)(p_over_Nv);
        hole_degfactor(cell,point) = p_over_Nv*exp(-eta);
      }

      // std::cout << "gamma_n=" << elec_degfactor(cell,point) <<", gamma_p=" << hole_degfactor(cell,point) << std::endl;

    } // end of point loop
  }  // end of cell loop
 }  // end of the if (bUseFD) block


 // otherwise, use the Boltzmann statistics and degeneracy factor = 1
 else
 {
  // loop over the cells
  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point = 0; point < num_points; ++point)
    {
      elec_degfactor(cell,point) = 1.0;
      hole_degfactor(cell,point) = 1.0;

    } // end of point loop
  }  // end of cell loop
 }  // end of the else block

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Degeneracy_Factor<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<bool>("Fermi Dirac", "false", "Use the Fermi-Dirac statistics if true");

  p->set<std::string>("FD Formula", "Schroeder", "Can be either Schroeder or Diffusion");

  return p;
}

}

#endif
