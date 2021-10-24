
#ifndef CHARON_RECOMBRATE_TOTAL_IMPL_HPP
#define CHARON_RECOMBRATE_TOTAL_IMPL_HPP

#include <cmath>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_Total<EvalT, Traits>::
RecombRate_Total(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using panzer::BasisIRLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // Get the equation set type
  const string eqset = p.get<string>("Equation Set Type");
  isSGCVFEM = false; 
  if (eqset == "SGCVFEM Drift Diffusion") isSGCVFEM = true; 

  // Get IR - It is FEM IP for SUPG-FEM and EFFPG-FEM, but is the SubCV centroid for SGCVFEM 
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> ip_scalar = ir->dl_scalar;
  num_points = ip_scalar->dimension(1);

  // Get Basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> basis_scalar = basis->functional;
  num_nodes = basis->functional->dimension(1);
  basis_name = basis->name();

  // Set the boolean values according to input from the closure model factory
  bSRH = false;
  bTrapSRH = false;
  bRadiative = false;
  bAuger = false;
  bAvalanche = false;
  bDefect = false;
  bEmpiricalDefect = false;
  bParticleStrike = false;
  bOptGen = false;
  if (p.get<string>("SRH") == "On") bSRH = true;
  if (p.get<string>("Trap SRH") == "On") bTrapSRH = true;
  if (p.get<string>("Radiative") == "On") bRadiative = true;
  if (p.get<string>("Auger") == "On") bAuger = true;
  if (p.get<string>("Avalanche") == "On") bAvalanche = true;
  if (p.get<string>("Defect Cluster") == "On") bDefect = true;
  if (p.get<string>("Empirical Defect") == "On") bEmpiricalDefect = true;
  if (p.get<string>("Particle Strike") == "On") bParticleStrike = true;
  if (p.get<string>("Optical Generation") == "On") bOptGen = true;

  // Evaluated fields (computed at IPs)
  total_rate = MDField<ScalarT,Cell,Point>(n.field.total_recomb, ip_scalar);
  total_deriv_e = MDField<ScalarT,Cell,Point>(n.field.recomb_deriv_e, ip_scalar);
  total_deriv_h = MDField<ScalarT,Cell,Point>(n.field.recomb_deriv_h, ip_scalar);
  this->addEvaluatedField(total_rate);
  this->addEvaluatedField(total_deriv_e);
  this->addEvaluatedField(total_deriv_h);

  // Determine the data layout for input scalar fields
  RCP<DataLayout> input_scalar = ip_scalar;    // for SUPG-FEM and EFFPG-FEM
  if (isSGCVFEM) input_scalar = basis_scalar;  // for SGCVFEM

  // Dependent fields
  if (bSRH)
  {
    srh_rate = MDField<const ScalarT,Cell,Point>(n.field.srh_recomb, input_scalar);
    srh_deriv_e = MDField<const ScalarT,Cell,Point>(n.field.srh_deriv_e, input_scalar);
    srh_deriv_h = MDField<const ScalarT,Cell,Point>(n.field.srh_deriv_h, input_scalar);
    this->addDependentField(srh_rate);
    this->addDependentField(srh_deriv_e);
    this->addDependentField(srh_deriv_h);
  }
  if (bRadiative)
  {
    rad_rate = MDField<const ScalarT,Cell,Point>(n.field.rad_recomb, input_scalar);
    rad_deriv_e = MDField<const ScalarT,Cell,Point>(n.field.rad_deriv_e, input_scalar);
    rad_deriv_h = MDField<const ScalarT,Cell,Point>(n.field.rad_deriv_h, input_scalar);
    this->addDependentField(rad_rate);
    this->addDependentField(rad_deriv_e);
    this->addDependentField(rad_deriv_h);
  }
  if (bAuger)
  {
    auger_rate = MDField<const ScalarT,Cell,Point>(n.field.auger_recomb, input_scalar);
    auger_deriv_e = MDField<const ScalarT,Cell,Point>(n.field.auger_deriv_e, input_scalar);
    auger_deriv_h = MDField<const ScalarT,Cell,Point>(n.field.auger_deriv_h, input_scalar);
    this->addDependentField(auger_rate);
    this->addDependentField(auger_deriv_e);
    this->addDependentField(auger_deriv_h);
  }
  if (bDefect)
  {
    defect_cluster_rate = MDField<const ScalarT,Cell,Point>(n.field.defect_cluster_recomb, input_scalar);
    this->addDependentField(defect_cluster_rate);
  }
  if (bEmpiricalDefect)
  {
    empirical_defect_rate = MDField<const ScalarT,Cell,Point>(n.field.empirical_defect_recomb, input_scalar);
    this->addDependentField(empirical_defect_rate);
  }
  if (bParticleStrike)
  {
    ionization_particle_strike_rate = MDField<const ScalarT,Cell,Point>(n.field.ionization_particle_strike_rate, input_scalar);
    this->addDependentField(ionization_particle_strike_rate);
  }
  if (bOptGen)
  {
    opt_gen_rate = MDField<const ScalarT,Cell,Point>(n.field.opt_gen, input_scalar);
    this->addDependentField(opt_gen_rate);
  }
  if (bAvalanche)
  {
    ava_rate = MDField<const ScalarT,Cell,Point>(n.field.avalanche_rate, ip_scalar);
    this->addDependentField(ava_rate);
  }
  if (bTrapSRH)
  {
    trap_srh_rate = MDField<const ScalarT,Cell,Point>(n.field.trap_srh_recomb, ip_scalar);
    trap_srh_deriv_e = MDField<const ScalarT,Cell,Point>(n.field.trap_srh_deriv_e, ip_scalar);
    trap_srh_deriv_h = MDField<const ScalarT,Cell,Point>(n.field.trap_srh_deriv_h, ip_scalar);
    this->addDependentField(trap_srh_rate);
    this->addDependentField(trap_srh_deriv_e);
    this->addDependentField(trap_srh_deriv_h);
  }

  std::string name = "Total_Recombination_Rate";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_Total<EvalT, Traits>::
postRegistrationSetup(
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
void
RecombRate_Total<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    Kokkos::DynRankView<ScalarT,PHX::Device> srh_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> srh_deriv_e_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> srh_deriv_h_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> rad_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> rad_deriv_e_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> rad_deriv_h_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> auger_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> auger_deriv_e_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> auger_deriv_h_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> defect_cluster_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> empirical_defect_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> ionization_particle_strike_rate_cpts;
    Kokkos::DynRankView<ScalarT,PHX::Device> opt_gen_rate_cpts;

    if (isSGCVFEM) 
    {
      // for SGCVFEM, interpolate all rates but avalanche and trap srh to centroids(IPs) from
      // the corresponding basis values; this is necessary because field-independent rates are 
      // not computed at IPs for SGCVFEM
      const int num_ips = num_points;

      // initialization
      if (bSRH) 
      {
        srh_rate_cpts = 
          createDynRankView(srh_rate.get_static_view(),"srh_rate_cpts",num_ips);
        Kokkos::deep_copy(srh_rate_cpts,ScalarT(0.0));
        srh_deriv_e_cpts =
          createDynRankView(srh_deriv_e.get_static_view(),"srh_deriv_e_cpts",num_ips);
        Kokkos::deep_copy(srh_deriv_e_cpts,ScalarT(0.0));
        srh_deriv_h_cpts =
          createDynRankView(srh_deriv_h.get_static_view(),"srh_deriv_h_cpts",num_ips);
        Kokkos::deep_copy(srh_deriv_h_cpts ,ScalarT(0.0));
      }
      if (bRadiative) 
      {
        rad_rate_cpts =
          createDynRankView(rad_rate.get_static_view(),"rad_rate_cpts",num_ips);
        Kokkos::deep_copy(rad_rate_cpts,ScalarT(0.0));
        rad_deriv_e_cpts =
          createDynRankView(rad_deriv_e.get_static_view(),"rad_deriv_e_cpts",num_ips);
        Kokkos::deep_copy(rad_deriv_e_cpts,ScalarT(0.0));
        rad_deriv_h_cpts =
          createDynRankView(rad_deriv_h.get_static_view(),"rad_deriv_h_cpts",num_ips);
        Kokkos::deep_copy(rad_deriv_h_cpts,ScalarT(0.0));
      }
      if (bAuger) 
      {
        auger_rate_cpts =
          createDynRankView(auger_rate.get_static_view(),"auger_rate_cpts",num_ips);
        Kokkos::deep_copy(auger_rate_cpts,ScalarT(0.0));
        auger_deriv_e_cpts =
          createDynRankView(auger_deriv_e.get_static_view(),"auger_deriv_e_cpts",num_ips);
        Kokkos::deep_copy(auger_deriv_e_cpts,ScalarT(0.0));
        auger_deriv_h_cpts =
          createDynRankView(auger_deriv_h.get_static_view(),"auger_deriv_h_cpts",num_ips);
        Kokkos::deep_copy(auger_deriv_h_cpts,ScalarT(0.0));
      }
      if (bDefect) 
      {
        defect_cluster_rate_cpts =
          createDynRankView(defect_cluster_rate.get_static_view(), "defect_cluster_rate_cpts",num_ips);
        Kokkos::deep_copy(defect_cluster_rate_cpts,ScalarT(0.0));
      }
      if (bEmpiricalDefect) 
      {
        empirical_defect_rate_cpts =
          createDynRankView(empirical_defect_rate.get_static_view(), "empirical_defect_rate_cpts",num_ips);
        Kokkos::deep_copy(empirical_defect_rate_cpts,ScalarT(0.0));
      }
      if (bParticleStrike) 
      {
        ionization_particle_strike_rate_cpts =
          createDynRankView(ionization_particle_strike_rate.get_static_view(),
                            "ionization_particle_strike_rate_cpts",num_ips);
        Kokkos::deep_copy(ionization_particle_strike_rate_cpts,ScalarT(0.0));
      }
      if (bOptGen) 
      {
        opt_gen_rate_cpts =
          createDynRankView(opt_gen_rate.get_static_view(),"opt_gen_rate_cpts",num_ips);
        Kokkos::deep_copy(opt_gen_rate_cpts,ScalarT(0.0));
      }

      // interpolate scalar fields at nodes to the centroids of subcv for SGCVFEM
      for (int inode = 0; inode < num_nodes; ++inode)
      {
        for (int ip = 0; ip < num_ips; ++ip) 
        {
          ScalarT tmp = (workset.bases[basis_index])->basis_scalar(cell,inode,ip);
          if (bSRH) 
          {
            srh_rate_cpts(ip) += tmp * srh_rate(cell,inode);
            srh_deriv_e_cpts(ip) += tmp * srh_deriv_e(cell,inode);
            srh_deriv_h_cpts(ip) += tmp * srh_deriv_h(cell,inode);
          }
          if (bRadiative) 
          {
            rad_rate_cpts(ip) += tmp * rad_rate(cell,inode);
            rad_deriv_e_cpts(ip) += tmp * rad_deriv_e(cell,inode);
            rad_deriv_h_cpts(ip) += tmp * rad_deriv_h(cell,inode);
          }
          if (bAuger) 
          {
            auger_rate_cpts(ip) += tmp * auger_rate(cell,inode);
            auger_deriv_e_cpts(ip) += tmp * auger_deriv_e(cell,inode);
            auger_deriv_h_cpts(ip) += tmp * auger_deriv_h(cell,inode);
          }
          if (bDefect) 
            defect_cluster_rate_cpts(ip) += tmp * defect_cluster_rate(cell,inode);
          if (bEmpiricalDefect)
            empirical_defect_rate_cpts(ip) += tmp * empirical_defect_rate(cell,inode);
          if (bParticleStrike)
            ionization_particle_strike_rate_cpts(ip) += tmp * ionization_particle_strike_rate(cell,inode);
          if (bOptGen) 
            opt_gen_rate_cpts(ip) += tmp * opt_gen_rate(cell,inode);
        }
      }
    }  // end of the node-to-centroid interpolation for SGCVFEM 

    for (int point = 0; point < num_points; ++point)
    {
      ScalarT total = 0.0;
      ScalarT e_deriv = 0.0;
      ScalarT h_deriv = 0.0;

      if (bSRH)  // add srh recomb.
      {
        if (isSGCVFEM) 
        {
          total += srh_rate_cpts(point);   
          e_deriv += srh_deriv_e_cpts(point);
          h_deriv += srh_deriv_h_cpts(point);
        } 
        else 
        {
          total += srh_rate(cell,point);   
          e_deriv += srh_deriv_e(cell,point);
          h_deriv += srh_deriv_h(cell,point);
        }
      }
      if (bRadiative)  // add radiative recomb.
      {
        if (isSGCVFEM) 
        {
          total += rad_rate_cpts(point);   
          e_deriv += rad_deriv_e_cpts(point);
          h_deriv += rad_deriv_h_cpts(point);
        } 
        else
        {
          total += rad_rate(cell,point);   
          e_deriv += rad_deriv_e(cell,point);
          h_deriv += rad_deriv_h(cell,point);
        }
      }
      if (bAuger)  // add auger recomb.
      {
        if (isSGCVFEM) 
        {
          total += auger_rate_cpts(point);  
          e_deriv += auger_deriv_e_cpts(point);
          h_deriv += auger_deriv_h_cpts(point);
        } 
        else 
        {
          total += auger_rate(cell,point);  
          e_deriv += auger_deriv_e(cell,point);
          h_deriv += auger_deriv_h(cell,point);
        }
      }
      if (bDefect)  // add defect cluster recomb.
      {
        if (isSGCVFEM)
          total += defect_cluster_rate_cpts(point);  
        else
          total += defect_cluster_rate(cell,point);  
      }
      if (bEmpiricalDefect)  // add empirical defect recomb.
      {
        if (isSGCVFEM)   
          total += empirical_defect_rate_cpts(point); 
        else
          total += empirical_defect_rate(cell,point); 
      }
      if (bParticleStrike)  // add empirical defect recomb.
      {
        if (isSGCVFEM)
          total += ionization_particle_strike_rate_cpts(point);  
        else
          total += ionization_particle_strike_rate(cell,point);  
      }
      if (bOptGen)  // subtract optical generation
      {
        if (isSGCVFEM)
          total -= opt_gen_rate_cpts(point); 
        else
          total -= opt_gen_rate(cell,point);  
      }

      if (bAvalanche)  // subtract avalanche generation
        total -= ava_rate(cell,point);   

      if (bTrapSRH)  // add trap srh recomb.
      {
        total += trap_srh_rate(cell,point); 
        e_deriv += trap_srh_deriv_e(cell,point); 
        h_deriv += trap_srh_deriv_h(cell,point); 
      }

      total_rate(cell,point) = total;
      total_deriv_e(cell,point) = e_deriv;
      total_deriv_h(cell,point) = h_deriv;
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
RecombRate_Total<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("SRH", "Off");
  p->set<std::string>("Trap SRH", "Off");
  p->set<std::string>("Radiative", "Off");
  p->set<std::string>("Auger", "Off");
  p->set<std::string>("Avalanche", "Off");
  p->set<std::string>("Defect Cluster", "Off");
  p->set<std::string>("Empirical Defect", "Off");
  p->set<std::string>("Particle Strike", "Off");
  p->set<std::string>("Optical Generation", "Off");
  p->set<std::string>("Equation Set Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  return p;
}

}

#endif
