
#ifndef CHARON_IC_EQUILIBRIUM_POTENTIAL_IMPL_HPP
#define CHARON_IC_EQUILIBRIUM_POTENTIAL_IMPL_HPP

#include "Panzer_Workset.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Charon_Names.hpp"
#include "Charon_Material_Properties.hpp"
#include "Charon_Physical_Constants.hpp"
#include "Panzer_STK_SetupUtilities.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
IC_Equilibrium_Potential<EvalT, Traits>::
IC_Equilibrium_Potential(
  const Teuchos::ParameterList& p) {
  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::BasisIRLayout;
  using Teuchos::ParameterList;

  m_names = p.get< Teuchos::RCP< const charon::Names> >("Names");

  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> scalar = basis->functional;
  num_basis = scalar->dimension(1);

  // Get device mesh
  dev_mesh = p.get<RCP<const panzer_stk::STK_Interface>>("Mesh");
  // Get Schottky contacts info
  sc_cnts = p.get<RCP<vector<string>>>("SchottkyContacts");
  sc_blks = p.get<RCP<vector<string>>>("SchottkyBlocks");
  sc_wf = p.get<RCP<vector<double>>>("SchottkyWF");
  sc_vapp = p.get<RCP<vector<double>>>("SchottkyVapp");
  
  RCP<stk::mesh::BulkData> bulkData = dev_mesh->getBulkData();
  for (std::size_t i=0;i<sc_cnts->size(); ++i) {
    const string cnt_name = (*sc_cnts)[i];
    const string blk_name = (*sc_blks)[i];
    sch_cnt_nodes.insert(std::make_pair(cnt_name, vector<stk::mesh::Entity>()));
    vector<stk::mesh::Entity> side_entities;
    dev_mesh->getAllSides(cnt_name,blk_name,side_entities);
    for(std::size_t et = 0; et < side_entities.size(); et++) { // 2D or 3D
      if( (dev_mesh->getDimension() == 2 and bulkData->entity_rank(side_entities[et]) == 
           stk::topology::EDGE_RANK) or (dev_mesh->getDimension() == 3 and
	   bulkData->entity_rank(side_entities[et]) == stk::topology::FACE_RANK) ) {
	stk::mesh::Entity const *et_nodes_i = bulkData->begin_nodes(side_entities[et]);
	stk::mesh::Entity const *et_nodes_e = bulkData->end_nodes(side_entities[et]);
	for(; et_nodes_i != et_nodes_e; ++et_nodes_i)
	  sch_cnt_nodes[cnt_name].push_back(*et_nodes_i);
      }
    }
  }

  // Get equation set type
  eqnSetType = p.get<string>("Equation Set Type");

  dof_name = p.get<string>("DOF Name");
  const string pot_name = "ELECTRIC_POTENTIAL";

  if(dof_name.length() == pot_name.length())
    haveSuffix = false;
  else
    haveSuffix = true;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters(haveSuffix);
  p.validateParameters(*valid_params);

  // Evaluated field
  potential = MDField<ScalarT,Cell,BASIS>(dof_name, scalar);
  this->addEvaluatedField(potential);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  T0 = scaleParams->scale_params.T0;
  V0 = scaleParams->scale_params.V0;


  // Dependent fields
  elec_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.elec_eff_dos, scalar);
  hole_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.hole_eff_dos, scalar);
  acceptor = MDField<const ScalarT,Cell,BASIS>(m_names->field.acceptor_raw, scalar);
  donor = MDField<const ScalarT,Cell,BASIS>(m_names->field.donor_raw, scalar);
  eff_affinity = MDField<const ScalarT,Cell,BASIS>(m_names->field.eff_affinity, scalar);
  eff_bandgap = MDField<const ScalarT,Cell,BASIS>(m_names->field.eff_band_gap, scalar);
  latt_temp = MDField<const ScalarT,Cell,BASIS>(m_names->field.latt_temp, scalar);
  if(eqnSetType == "Laplace" ||
     eqnSetType == "SGCVFEM Laplace" ||
     eqnSetType == "NLPoisson" ||
     eqnSetType == "SGCVFEM NLPoisson" ||
     eqnSetType == "Drift Diffusion" ||
     eqnSetType == "EFFPG Drift Diffusion" ||
     eqnSetType == "SGCVFEM Drift Diffusion" ||
     eqnSetType == "SGCharon1 Drift Diffusion" ||
     eqnSetType == "DDIon") {
    ref_energy = MDField<const ScalarT,Cell,BASIS>(m_names->field.ref_energy, scalar);
    this->addDependentField(ref_energy);
  }

  this->addDependentField(elec_effdos);
  this->addDependentField(hole_effdos);
  this->addDependentField(acceptor);
  this->addDependentField(donor);
  this->addDependentField(eff_affinity);
  this->addDependentField(eff_bandgap);
  this->addDependentField(latt_temp);


  std::string name = "IC_Equilibrium_Potential";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
IC_Equilibrium_Potential<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset) {
  const std::vector<stk::mesh::Entity>& localElements = *dev_mesh->getElementsOrderedByLID();
  const std::vector<std::size_t>& localCellIds = workset.cell_local_ids;

  // obtain physical constants
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double kb = cpc.kb; // Boltzmann constant in [eV/K]
  bool isothermal = false;
  if(eqnSetType == "Laplace" ||
     eqnSetType == "SGCVFEM Laplace" ||
     eqnSetType == "NLPoisson" ||
     eqnSetType == "SGCVFEM NLPoisson" ||
     eqnSetType == "Drift Diffusion" ||
     eqnSetType == "EFFPG Drift Diffusion" ||
     eqnSetType == "SGCVFEM Drift Diffusion" ||
     eqnSetType == "SGCharon1 Drift Diffusion" ||
     eqnSetType == "DDIon" ||
     eqnSetType == "DDLattice" )
    isothermal = true;
  else if(eqnSetType == "DDIonLattice" ||
          eqnSetType == "EFFPG DDIonLattice")
    isothermal = false;

  for(int cell = 0; cell < workset.num_cells; ++cell) {
    std::size_t cellLocalId = localCellIds[cell];
    stk::mesh::Entity const* relations = 
      dev_mesh->getBulkData()->begin_nodes(localElements[cellLocalId]);

    for(int basis = 0; basis < num_basis; ++basis) {
      ScalarT sch_WF = 0.0;
      ScalarT sch_Vapp = 0.0;
      stk::mesh::Entity node = relations[basis];
      // check to see if node is a Schottky contact
      for (std::size_t i=0;i<sc_cnts->size(); ++i) {
	const string cnt_name = (*sc_cnts)[i];
	std::vector<stk::mesh::Entity>::iterator it = 
	  std::find(sch_cnt_nodes[cnt_name].begin(),sch_cnt_nodes[cnt_name].end(),node);  
	if(it != sch_cnt_nodes[cnt_name].end()) {
	  // retrieve info about Schottky contact
	  sch_WF = (*sc_wf)[i];
	  sch_Vapp = (*sc_vapp)[i];
	  break;
	}
      }

      ScalarT kbT = latt_temp(cell,basis) * T0 * kb; // [eV]
      const ScalarT& Nc = elec_effdos(cell,basis); // scaled
      const ScalarT& Nv = hole_effdos(cell,basis); // scaled
      const ScalarT& Chieff = eff_affinity(cell,basis); // [eV]
      const ScalarT& Egeff = eff_bandgap(cell,basis); // [eV]
      const ScalarT dop = donor(cell,basis) - acceptor(cell,basis); // scaled

      if(isothermal) { // isothermal
        const ScalarT Eref = ref_energy(0,0);
        if(dop >= 0) { // n-type
	  if (sch_WF > 0.0) {
	    potential(cell,basis) = Eref - sch_WF + sch_Vapp; // [V]
	  } else {
	    ScalarT y = 0.5*dop/Nc +
	      std::sqrt( (0.5*dop/Nc)*(0.5*dop/Nc) + Nv/Nc*std::exp(-Egeff/kbT) );
	    potential(cell,basis) = Eref - Chieff + kbT*std::log(y); // [V]  
	  }
        } else { // p-type
	  if (sch_WF > 0.0) {
	    potential(cell,basis) = Eref - sch_WF + sch_Vapp; // [V]
	  } else {
	    ScalarT y = -0.5*dop/Nv +
            std::sqrt( (0.5*dop/Nv)*(0.5*dop/Nv) + Nc/Nv*std::exp(-Egeff/kbT) );
	    potential(cell,basis) = Eref - Chieff - Egeff - kbT*std::log(y);
	  }
        }
        potential(cell,basis) /= V0; // scaled
      } else { // non-isothermal
        ScalarT qtheta = Chieff + 0.5*Egeff + 0.5*kbT*std::log(Nc/Nv);  // [eV]
        if(dop >= 0) { // n-type
	  if (sch_WF > 0.0) {
	    potential(cell,basis) = qtheta - sch_WF + sch_Vapp; // [V]
	  } else { 
	    ScalarT y = 0.5*dop/Nc +
	      std::sqrt( (0.5*dop/Nc)*(0.5*dop/Nc) + Nv/Nc*std::exp(-Egeff/kbT) );
	    potential(cell,basis) = qtheta - Chieff + kbT*std::log(y); // [V]
	  }
        } else { // p-type
	  if (sch_WF > 0.0) {
	    potential(cell,basis) = qtheta - sch_WF + sch_Vapp; // [V]
	  } else {
	    ScalarT y = -0.5*dop/Nv +
	      std::sqrt( (0.5*dop/Nv)*(0.5*dop/Nv) + Nc/Nv*std::exp(-Egeff/kbT));
	    potential(cell,basis) = qtheta - Chieff - Egeff - kbT*std::log(y); // [V]
	  }
        }
        potential(cell,basis) /= V0; // scaled
      }
    } // loop over basis
  } // loop over cell
}

///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
IC_Equilibrium_Potential<EvalT, Traits>::getValidParameters(bool /* haveSuffix */) const {
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("DOF Name", "?");

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  p->set("Equation Set Type", "?");

  p->set<bool>("Fermi Dirac", "false", "Use the Fermi-Dirac statistics if true");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  Teuchos::RCP<const panzer_stk::STK_Interface> mesh;
  p->set("Mesh", mesh);

  Teuchos::RCP<std::vector<std::string>> sc_cnts;
  p->set("SchottkyContacts", sc_cnts);

  Teuchos::RCP<std::vector<std::string>> sc_blks;
  p->set("SchottkyBlocks", sc_blks);

  Teuchos::RCP<std::vector<double>> sc_wf;
  p->set("SchottkyWF", sc_wf);

  Teuchos::RCP<std::vector<double>> sc_vapp;
  p->set("SchottkyVapp", sc_vapp);
  
  return p;
}

}

#endif

