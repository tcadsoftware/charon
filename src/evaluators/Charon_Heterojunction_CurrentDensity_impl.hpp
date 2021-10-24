
#ifndef CHARON_HETEROJUNCTION_CURRENTDENSITY_IMPL_HPP
#define CHARON_HETEROJUNCTION_CURRENTDENSITY_IMPL_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellTopologyInfo.hpp"
#include "Shards_CellTopology.hpp"

#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Charon_Names.hpp"

namespace charon {

/**
 * @brief This evaluator computes the thermionic emission (TE) current density
 * normal to a heterojunction due to contribution from one side. The net normal
 * TE current density is obtained by calling this evaluator twice (due to two
 * sides of a HJ) and subtracting the two values using the charon::Subtract
 * evaluator.
 * The TE model is implemented for both FEM-SUPG and CVFEM-SG discretization
 * schemes. For FEM-SUGP, the TE current density is computed at the FEM integration
 * points (IPs) of a HJ; while for CVFEM-SG, the TE current density is computed
 * at the CVFEM IPs of a HJ.
*/

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Heterojunction_CurrentDensity<EvalT, Traits>::
Heterojunction_CurrentDensity(
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

  const charon::Names& n = *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // Retrieve the data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Integration rule
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  int_rule_degree = ir->cubature_degree;
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  num_ips = ip_scalar->dimension(1);

  // Obtain the BASIS layout
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  basis_name = basis->name();
  num_nodes = basis_scalar->dimension(1);

  // Retrieve parameters
  dof_name = p.get<string>("DOF Name");
  other_dof_name = p.get<string>("Other DOF Name");
  flux_name = p.get<string>("Flux Name");
  richConst = p.get<double>("Richardson Constant");
  bandOffset = p.get<double>("Band Offset");
  fdDensity = p.get<double>("Fermi Dirac Density");
  detailIndex = p.get<int>("Details Index");
  discMethod = p.get<string>("Discretization Method");

  if (p.get<string>("Primary Side") == "Left")
    isPrimarySideLeft = true;
  else
    isPrimarySideLeft = false;

  // Determine the effective density of states
  std::size_t foundedens = dof_name.find("ELECTRON_DENSITY");
  std::size_t foundhdens = dof_name.find("HOLE_DENSITY");
  string eff_dos_name = "";

  if (foundedens != string::npos)       // for ELECTRON_DENSITY
    eff_dos_name = n.field.elec_eff_dos;
  else if (foundhdens != string::npos)  // for HOLE_DENSITY
    eff_dos_name = n.field.hole_eff_dos;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error ! dof_name must contain either ELECTRON_DENSITY or HOLE_DENSITY !");

  // Evaluated field
  hj_flux = MDField<ScalarT,Cell,IP>(flux_name, ip_scalar);
  this->addEvaluatedField(hj_flux);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  C0 = scaleParams->scale_params.C0;
  J0 = scaleParams->scale_params.J0;

  // Dependent field
  carr_dens = MDField<const ScalarT,Cell,Point>(dof_name, scalar);
  other_carr_dens = MDField<const ScalarT,Cell,Point>(other_dof_name, scalar);
  eff_dos = MDField<const ScalarT,Cell,Point>(eff_dos_name, scalar);
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp, scalar);

  this->addDependentField(carr_dens);
  this->addDependentField(other_carr_dens);
  this->addDependentField(eff_dos);
  this->addDependentField(latt_temp);

  std::string name = "Heterojunction_CurrentDensity";
  this->setName(name);

  // instantiate the FermiDiracIntegral class
  inverseFermiIntegral =
    Teuchos::rcp(new charon::FermiDiracIntegral<EvalT>(charon::FermiDiracIntegral<EvalT>::inverse_PlusOneHalf));

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Heterojunction_CurrentDensity<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0], this->wda);
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0], this->wda);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Heterojunction_CurrentDensity<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Obtain physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb;   // [eV/K]

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    Kokkos::DynRankView<ScalarT,PHX::Device> hj_flux_point = Kokkos::createDynRankView(hj_flux.get_static_view(),"hj_flux_point",num_points);

    for (int point = 0; point < num_points; ++point)
    {
      // obtain lattice temperature in [K]
      ScalarT latT = latt_temp(cell,point) * T0;

      // lattT should be always > 0, but it could become <= 0 due to numerical errors
      // when the temperature eqn is solved, so reset it to 300 K to avoid unphysical parameters
      if (Sacado::ScalarValue<ScalarT>::eval(latT) <= 0.0)  latT = 300.0;

      ScalarT kbT = kb * latT;  // [eV]

      // obtain carrier density and effective dos in [cm^-3]
      const ScalarT& carrdens = carr_dens(cell,point) * C0;
      const ScalarT& other_carrdens = other_carr_dens(cell,point) * C0;
      const ScalarT& effdos = eff_dos(cell,point) * C0;

      hj_flux_point(point) = 0.0; // equal to 0.0 when carrdens <= 0.0 or other_carrdens <= 0.0

      // compute hj_flux [scaled] without the exponential factor
      if (Sacado::ScalarValue<ScalarT>::eval(carrdens) > std::numeric_limits<double>::epsilon() &&
          Sacado::ScalarValue<ScalarT>::eval(other_carrdens) > std::numeric_limits<double>::epsilon() )
      {
        // this logic makes more sense than he3
        if (Sacado::ScalarValue<ScalarT>::eval(carrdens) < fdDensity)  // use Boltzmann statistics
          hj_flux_point(point) = richConst * latT * latT * carrdens / effdos / J0;

        // note that std::exp(eta) could become infinity when carrdens becomes unphysically large during iteration
        else if (Sacado::ScalarValue<ScalarT>::eval(carrdens) <= (effdos*100.0))  // use Fermi-Dirac statistis
        {
          ScalarT eta = (*inverseFermiIntegral)(carrdens/effdos);
          hj_flux_point(point) = richConst * latT * latT * std::exp(eta) / J0;
        }

        else  // use the Boltzmann expression when carrdens becomes unphysically large due to numerical instability during Newton iteration
          hj_flux_point(point) = richConst * latT * latT * carrdens / effdos / J0;
      }

      ScalarT expFactor = 1.0;

      // determine when to add the exponential factor
      if (bandOffset >= 0.0)
      {
        if (isPrimarySideLeft)  // Primary Side is on the left, di=1 corresponds to side 2 (right side)
        {
          if (detailIndex == 1)  // add to side 2 when bandOffset >= 0.0
            expFactor = std::exp(-bandOffset / kbT);  // unitless
        }
        else  // Primary Side is on the right, di=0 corresponds to side 2 (right side)
        {
          if (detailIndex == 0)  // add to side 2 when bandOffset >= 0.0
            expFactor = std::exp(-bandOffset / kbT);  // unitless
        }
      }

      else  // negative bandOffset
      {
        if (isPrimarySideLeft)  // Primary Side is on the left, di=0 corresponds to side 1 (left side)
        {
          if (detailIndex == 0)  // add to side 1 when bandOffset < 0.0
            expFactor = std::exp(bandOffset / kbT);  // unitless
        }
        else  // Primary Side is on the right, di=1 corresponds to side 1 (left side)
        {
          if (detailIndex == 1)  // add to side 1 when bandOffset < 0.0
            expFactor = std::exp(bandOffset / kbT);  // unitless
        }
      }

      // multiply the exponential factor
      hj_flux_point(point) = hj_flux_point(point) * expFactor;

    }  // end of loop over points

    for (int ip = 0; ip < num_ips; ++ip)
    {

       if (discMethod == femsupg)
          hj_flux(cell,ip) = hj_flux_point(ip);
       else if (discMethod == cvfemsg)
       {
          ScalarT temp_flux = 0;
          for (int basis = 0; basis < num_nodes; ++basis)
          {
             temp_flux += hj_flux_point(basis) * (workset.details(detailIndex).bases[basis_index])->basis_scalar(cell, basis, ip);

             // double xn = workset.details(detailIndex).cell_vertex_coordinates(cell,basis,0);
             // double yn = workset.details(detailIndex).cell_vertex_coordinates(cell,basis,1);
             // std::cout << "c=" << cell << ",ip=" << ip << ", node=" << basis <<", xn= " << xn << ", yn= "<< yn <<", "<< isPrimarySideLeft << ",di=" << detailIndex << ", " << std::setprecision(20) << ", basis = " <<(workset.details(detailIndex).bases[basis_index])->basis_scalar(cell,basis,ip) <<", hj_flux=" << hj_flux_point(basis) << std::endl;
          }
          hj_flux(cell,ip) = temp_flux;
       }

      // get ip coordinates for debugging
      // double x = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,0);
      // double y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
      // const panzer::IntegrationValues2<double> &iv = *this->wda(workset).int_rules[int_rule_index];
      // double x = iv.ip_coordinates(cell,ip,0);
      // double y = iv.ip_coordinates(cell,ip,1);

      // std::cout << std::setprecision(20) <<"c=" << cell << ", ip=" << ip << ", x=" << x << ", y=" << y << ", " << isPrimarySideLeft << ", di=" << detailIndex << ", " << dof_name << ", " << flux_name << ", hj_flux=" <<  hj_flux(cell,ip) << std::endl;

    } // end loop over ips

  }  // end of loop over cells

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Heterojunction_CurrentDensity<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->set<std::string>("DOF Name", "");
  p->set<std::string>("Other DOF Name", "");
  p->set<std::string>("Flux Name", "", "Determine if HJ_CurrentDensity_Side1 or HJ_CurrentDensity_Side2");
  p->set<std::string>("Discretization Method", "");

  p->set<double>("Richardson Constant", 0., "Richardson constant is in unit of [A/(cm^2.K^2)]");
  p->set<double>("Band Offset", 0., "Either conduction or valence band offset in unit of [eV], can be a positive or negative value");
  p->set<double>("Fermi Dirac Density", 1e12, "Carrier density above which the Fermi-Dirac thermionic emission expression is used");
  p->set<int>("Details Index", 0, "Determine if on side 1 or on side 2");

  p->set<std::string>("Primary Side", "", "Determine if Primary Side is Left or Right");

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
