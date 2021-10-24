
#ifndef CHARON_KIMPTONTID_DECL_HPP
#define CHARON_KIMPTONTID_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_interp.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Panzer_STK_GatherFields.hpp"

using panzer::Cell;
using panzer::Point;
using panzer::Dim;

namespace charon {

//! compute TID charge in insulators
template<typename EvalT, typename Traits>
class KimptonTID
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    KimptonTID(const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // output
  // insulator hole trapped charge density
  PHX::MDField<ScalarT,Cell,Point> ins_trappedcharge;

  // insulator generated pairs density
  PHX::MDField<ScalarT,Cell,Point> ins_genpair_dens;

  // input
  PHX::MDField<const ScalarT,Cell,Point,Dim> grad_pot;
  PHX::MDField<const ScalarT,Cell,Point> potential;
  double dose;
  double DEF;
  double Eform;
  double pow_dep;
  bool isSGCVFEM;

  bool WithInterfaceTraps;
  std::string sidesetID;
  std::string blockID;
  double Nti;
  double fill_facti;
  double sigma_i;

  bool WithVolumeTraps;
  double Ntv;
  double fill_factv;
  double sigma_v;
  ScalarT trap_vol; 
  double sigma_v_crit;
  ScalarT trap_crit_vol; 
 
  bool withVaryingV;
  double V_freeze;

  // material properties
  double mass_dens;
  
  // geometry/topology
  Teuchos::RCP<const panzer_stk::STK_Interface> dev_mesh;
  std::map<stk::mesh::Entity,double> area_map;
  std::map<stk::mesh::Entity,std::vector<double> > normal_map;
  
  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ips;
  int num_dims;
  int num_basis;
  int num_edges;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  // for basis points
  std::string hcurl_basis_name;
  std::size_t hcurl_basis_index;

  // primary cell topology
  Teuchos::RCP<const shards::CellTopology> cellType;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0;  // concentration scaling
  double E0;  // electric field scaling

  double gen_const;

  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
 
  Teuchos::RCP<Teuchos::Comm<int> const> comm;
  double ins_vol; // insulator volume in um^3

  void computeCentroidField(typename Traits::EvalData& workset, 
                            panzer::index_t cell, std::vector<ScalarT>& E);

  void comp_geo_info(const Teuchos::RCP<const panzer_stk::STK_Interface> mesh);
  
  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

}; // end of class KimptonTID


}


#endif
