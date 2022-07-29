
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
  Material_Properties& matProperty = Material_Properties::getInstance();

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

  // Get gate contacts info
  g_cnts = p.get<RCP<vector<string>>>("GateContacts");
  g_blks = p.get<RCP<vector<string>>>("GateBlocks");
  g_wf = p.get<RCP<vector<double>>>("GateWF");
  g_vapp = p.get<RCP<vector<double>>>("GateVapp");

  for (std::size_t i=0;i<g_cnts->size(); ++i) {
    const string cnt_name = (*g_cnts)[i];
    const string blk_name = (*g_blks)[i];
    g_cnt_nodes.insert(std::make_pair(cnt_name, vector<stk::mesh::Entity>()));
    vector<stk::mesh::Entity> side_entities;
    dev_mesh->getAllSides(cnt_name,blk_name,side_entities);
    for(std::size_t et = 0; et < side_entities.size(); et++) { // 2D or 3D
      if( (dev_mesh->getDimension() == 2 and bulkData->entity_rank(side_entities[et]) == 
           stk::topology::EDGE_RANK) or (dev_mesh->getDimension() == 3 and
	   bulkData->entity_rank(side_entities[et]) == stk::topology::FACE_RANK) ) {
	stk::mesh::Entity const *et_nodes_i = bulkData->begin_nodes(side_entities[et]);
	stk::mesh::Entity const *et_nodes_e = bulkData->end_nodes(side_entities[et]);
	for(; et_nodes_i != et_nodes_e; ++et_nodes_i)
	  g_cnt_nodes[cnt_name].push_back(*et_nodes_i);
      }
    }
  }

  // Get equation set type
  eqnSetType = p.get<string>("Equation Set Type");

  blockId = p.get<string>("Block ID");
  matName = p.get<string>("Material");
  matType = matProperty.getMaterialType(matName);
  semBlockIds = p.get<RCP<vector<string>>>("Semiconductor Blocks");
  semMaterials = p.get<RCP<vector<string>>>("Semiconductor Materials");

  if (matType == "Insulator" and semBlockIds->size() > 0) {
    for(std::size_t bi = 0; bi < semBlockIds->size(); bi++) {
      string semBlockId = (*semBlockIds)[bi];
      string semMaterial = (*semMaterials)[bi];
      std::set<stk::mesh::Entity> inter_nodes_set;
      std::vector<stk::mesh::Entity> ins_elems;
      std::vector<stk::mesh::Entity> sem_elems;
      dev_mesh->getMyElements(blockId,ins_elems);
      dev_mesh->getMyElements(semBlockId,sem_elems);
      for(std::size_t el = 0; el < ins_elems.size(); el++) {
	if(dev_mesh->getDimension() == 2) { // 2D
	  stk::mesh::Entity const *edge_i = bulkData->begin_edges(ins_elems[el]);
	  stk::mesh::Entity const *edge_e = bulkData->end_edges(ins_elems[el]);
	  for(; edge_i != edge_e; ++edge_i) {
	    for(std::size_t el1 = 0; el1 <sem_elems.size(); el1++) {
	      stk::mesh::Entity const *edge1_i = bulkData->begin_edges(sem_elems[el1]);
	      stk::mesh::Entity const *edge1_e = bulkData->end_edges(sem_elems[el1]);
	      bool found = false;
	      for(; edge1_i != edge1_e; ++edge1_i) {
		if (* edge1_i == *edge_i) {
		  stk::mesh::Entity const *node_i = bulkData->begin_nodes(*edge_i);
		  stk::mesh::Entity const *node_e = bulkData->end_nodes(*edge_i);
		  for(; node_i != node_e; ++node_i) {
		    inter_nodes_set.insert(*node_i);
		  }
		  found = true;
		  break;
		}
	      }
	      if (found) break;
	    }
	  }
	} else if(dev_mesh->getDimension() == 3) { // 3D
	  stk::mesh::Entity const *face_i = bulkData->begin_faces(ins_elems[el]);
	  stk::mesh::Entity const *face_e = bulkData->end_faces(ins_elems[el]);
	  for(; face_i != face_e; ++face_i) {
	    for(std::size_t el1 = 0; el1 <sem_elems.size(); el1++) {
	      stk::mesh::Entity const *face1_i = bulkData->begin_faces(sem_elems[el1]);
	      stk::mesh::Entity const *face1_e = bulkData->end_faces(sem_elems[el1]);
	      bool found = false;
	      for(; face1_i != face1_e; ++face1_i) {
		if (* face1_i == *face_i) {
		  stk::mesh::Entity const *node_i = bulkData->begin_nodes(*face_i);
		  stk::mesh::Entity const *node_e = bulkData->end_nodes(*face_i);
		  for(; node_i != node_e; ++node_i) {
		    inter_nodes_set.insert(*node_i);
		  }
		  found = true;
		  break;
		}
	      }
	      if (found) break;
	    }
	  }
	}
      }
      inter_nodes.push_back(inter_nodes_set);

      Nc300.push_back(matProperty.getPropertyValue(semMaterial, "Electron Effective DOS at 300 K"));
      Nv300.push_back(matProperty.getPropertyValue(semMaterial, "Hole Effective DOS at 300 K"));
      Nc_F.push_back(matProperty.getPropertyValue(semMaterial, "Electron Effective DOS Exponent"));
      Nv_F.push_back(matProperty.getPropertyValue(semMaterial, "Hole Effective DOS Exponent"));
    } // semBlockId
  }

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
  C0 = scaleParams->scale_params.C0;

  // Dependent fields  
  latt_temp = MDField<const ScalarT,Cell,BASIS>(m_names->field.latt_temp, scalar);
  this->addDependentField(latt_temp);

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

  if (matType == "Semiconductor") {
    elec_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.elec_eff_dos, scalar);
    hole_effdos = MDField<const ScalarT,Cell,BASIS>(m_names->field.hole_eff_dos, scalar);
    acceptor = MDField<const ScalarT,Cell,BASIS>(m_names->field.acceptor_raw, scalar);
    donor = MDField<const ScalarT,Cell,BASIS>(m_names->field.donor_raw, scalar);
    eff_affinity = MDField<const ScalarT,Cell,BASIS>(m_names->field.eff_affinity, scalar);
    eff_bandgap = MDField<const ScalarT,Cell,BASIS>(m_names->field.eff_band_gap, scalar);
    this->addDependentField(elec_effdos);
    this->addDependentField(hole_effdos);
    this->addDependentField(eff_affinity);
    this->addDependentField(acceptor);
    this->addDependentField(donor);
    this->addDependentField(eff_bandgap);
  } else if (matType == "Insulator") {
    eff_affinity = MDField<const ScalarT,Cell,BASIS>(m_names->field.affinity, scalar);
    eff_bandgap = MDField<const ScalarT,Cell,BASIS>(m_names->field.band_gap, scalar);
    this->addDependentField(eff_affinity);
    this->addDependentField(eff_bandgap);
    acceptor = MDField<const ScalarT,Cell,BASIS>(m_names->field.acceptor_raw, scalar);
    donor = MDField<const ScalarT,Cell,BASIS>(m_names->field.donor_raw, scalar);
    this->addDependentField(acceptor);
    this->addDependentField(donor);
  }
  
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

    // semiconductor
    if (matType == "Semiconductor" ) {
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
    } else {
      // insulator
      for(int basis = 0; basis < num_basis; ++basis) {
	stk::mesh::Entity node = relations[basis];

	bool isNodeOnContact = false;
	ScalarT bcValue = 0.0;
	for (std::size_t i=0; i<g_cnts->size(); ++i) {
	  const string cnt_name = (*g_cnts)[i];
	  std::vector<stk::mesh::Entity>::iterator it = 
	    std::find(g_cnt_nodes[cnt_name].begin(),g_cnt_nodes[cnt_name].end(),node);
	  if (it != g_cnt_nodes[cnt_name].end()) {
	    // node on a gate contact
	    isNodeOnContact = true;
	    // retrieve info about gate contact
	    ScalarT g_WF = (*g_wf)[i];
	    ScalarT g_Vapp = (*g_vapp)[i];  
	    ScalarT offsetDueToWF = (g_WF - ref_energy(0,0))/1.0;  // 1.0 converts from [eV] to [V]
	    bcValue = (g_Vapp - offsetDueToWF)/V0;
	    break;
	  }
	}

	bool isNodeOnInsSemInterface = false;
	ScalarT interValue = 0.0;
	if (semBlockIds->size() > 0) { 
	  for(std::size_t bi = 0; bi < semBlockIds->size(); bi++) {
	    // semiconductor/insulator interface
	    if ( inter_nodes[bi].find(node) != inter_nodes[bi].end() ) {
	      isNodeOnInsSemInterface = true;
	      ScalarT lattT = latt_temp(cell,basis) * T0;
	      ScalarT kbT = lattT * kb; // [eV]
	      const ScalarT& Chi = eff_affinity(cell,basis); // [eV]
	      const ScalarT& Eg = eff_bandgap(cell,basis); // [eV]
	      const ScalarT Eref = ref_energy(0,0);
	      const ScalarT dop = donor(cell,basis) - acceptor(cell,basis); // scaled
	      const ScalarT Nc = Nc300[bi] * pow(lattT/300.0, Nc_F[bi]) / C0; // scaled
	      const ScalarT Nv = Nv300[bi] * pow(lattT/300.0, Nv_F[bi]) / C0; // scaled
	      if(dop >= 0) { // n-type
		ScalarT y = 0.5*dop/Nc +
		  std::sqrt( (0.5*dop/Nc)*(0.5*dop/Nc) + Nv/Nc*std::exp(-Eg/kbT) );
		interValue = Eref - Chi + kbT*std::log(y); // [V] 
	      } else { // p-type
		ScalarT y = -0.5*dop/Nv +
		  std::sqrt( (0.5*dop/Nv)*(0.5*dop/Nv) + Nc/Nv*std::exp(-Eg/kbT) );
		interValue = Eref - Chi - Eg - kbT*std::log(y);
	      }
	    } // semiconductor/insulator interface
	  } // semiconductor/insulator interfaces
	}

	ScalarT val = 0.0;

	if (isNodeOnContact)
	  val = bcValue;
	if (isNodeOnInsSemInterface) 
	  val = interValue;
	potential(cell,basis) = val;
      } // insulator nodes
    } // insulator
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

  Teuchos::RCP<std::vector<std::string>> g_cnts;
  p->set("GateContacts", g_cnts);

  Teuchos::RCP<std::vector<std::string>> g_blks;
  p->set("GateBlocks", g_blks);

  Teuchos::RCP<std::vector<double>> g_wf;
  p->set("GateWF", g_wf);

  Teuchos::RCP<std::vector<double>> g_vapp;
  p->set("GateVapp", g_vapp);

  p->set("Material", "?");

  p->set<std::string>("Block ID", "?");

  Teuchos::RCP<std::vector<std::string>> semBlocks;
  p->set("Semiconductor Blocks", semBlocks);

  Teuchos::RCP<std::vector<std::string>> semMaterials;
  p->set("Semiconductor Materials", semMaterials);
  
  return p;
}

}

#endif

