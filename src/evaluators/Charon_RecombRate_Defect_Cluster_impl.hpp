
#ifndef CHARON_RECOMBRATE_DEFECT_CLUSTER_IMPL_HPP
#define CHARON_RECOMBRATE_DEFECT_CLUSTER_IMPL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <cmath>
#include <iostream>
#include <typeinfo> // Part of the cluster KLUDGE - LCM

// Charon
#include "Charon_Names.hpp"
#include <Charon_Vector.hpp>

// Panzer
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
RecombRate_Defect_Cluster<EvalT, Traits>::
RecombRate_Defect_Cluster(
  const Teuchos::ParameterList& p)
{
  using charon::Vector;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;
  using PHX::DataLayout;
  using PHX::MDField;
  using std::string;
  using Teuchos::RCP;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  isIPset = p.get<bool>("Is IP Set");

  if(isIPset)
  {
    //IP
    RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
    int_rule_degree = ir->cubature_degree;
  }
  else
  {
    // basis
    RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
    RCP<DataLayout> data_layout = basis->functional;
    basis_name = basis->name();
  }

  interpolator = p.get<Teuchos::RCP<charon::clusterInterpolator> >("cluster interpolator");


  // Input from closure model factory
  methodName = p.get<std::string>("Interpolant Method");
  //Set default Shepard Power, should typically be 1 < p < 3

  shepardPneg = 2.0;
  shepardPneg = p.get<double>("Shepard Power");
  bool set_Method = interpolator->setMethod(methodName,shepardPneg);
  if (!set_Method)
  {
    std::stringstream msg;
    msg << "Failed to set interpolation method: " << methodName << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!set_Method, std::logic_error, msg.str());
  }
  influenceRadius = p.get<double>("Influence Radius");

  inputType = p.get<std::string>("Input File Type");
  bool fileTypeSet = interpolator->setFileType(inputType);
  if (!fileTypeSet)
  {
    std::stringstream msg;
    msg << "Unable to set file type: " << inputType << " Check type."<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!fileTypeSet, std::logic_error, msg.str());
  }
  if(p.isParameter("Number of Input Files"))
    numberOfFiles = p.get<int>("Number of Input Files");

  cascadeDensity = p.get<double>("Cascade Density");
  if (!cascadeDensity)
  {
    std::stringstream msg;
    msg << "No cascade density set"<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!cascadeDensity, std::logic_error, msg.str());
  }

  //Set cascade density to 1.0 if in situ clusters are present
  if(inputType == "IN SITU")
    cascadeDensity = 1.0;

  // Read Cluster data from file
  // vector of filenames that contain data
  std::vector<std::string> fileNames;
  std::stringstream filenoSS;
  std::string fileno;

  for (int i=1 ; i<numberOfFiles+1 ; ++i)
  {
    filenoSS<<"ClusterData/Cluster1D."<<i<<".dat";
    fileno = filenoSS.str();
    fileNames.push_back(fileno.c_str());
    filenoSS.str("\0");
  }

  if(inputType != "IN SITU")
    interpolator->readFiles(fileNames);

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);


  // fields
  defect_cluster_rate = MDField<ScalarT,Cell,Point>(n.field.defect_cluster_recomb,scalar);

  intrin_conc = MDField<const ScalarT,Cell,Point>(n.field.intrin_conc,scalar);
  edensity = MDField<const ScalarT,Cell,Point>(n.dof.edensity,scalar);
  hdensity = MDField<const ScalarT,Cell,Point>(n.dof.hdensity,scalar);

  //THis is part of the cluster KLUDGE  LCM
  acceptor = MDField<const ScalarT,Cell,Point>(n.field.acceptor,scalar);
  donor = MDField<const ScalarT,Cell,Point>(n.field.donor,scalar);

  // Evaluated fields
  this->addEvaluatedField(defect_cluster_rate);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  C0 = scaleParams->scale_params.C0;

  // Dependent fields
  this->addDependentField(intrin_conc);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);
  //THis is part of the cluster KLUDGE  LCM
  this->addDependentField(acceptor);
  this->addDependentField(donor);


  std::string name = "Defect_Cluster_Recombination_Rate";
  this->setName(name);

}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_Defect_Cluster<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if (isIPset)
  {
    int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  }
  else
  {
    basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
  }


}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
RecombRate_Defect_Cluster<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // retrieve time scaling factor in [s]
  double tscale = Sacado::ScalarValue<ScalarT>::eval(t0);

  // get the current (present) time in [s]
  double curr_time = workset.time * tscale;

  double sinkLifetime;
  ScalarT scale_const = t0/C0;
  int dimension = workset.subcell_dim;

  double xn = 0.0;
  double yn = 0.0;
  double zn = 0.0;

  TEUCHOS_TEST_FOR_EXCEPTION(dimension == 3, std::logic_error, "The defect cluster recombination code "\
                             "is not currently valid for three-dimensional geometries");


  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {

    for (int point = 0; point < num_points; ++point)
    {
      if (isIPset)
      {
        xn = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
        yn = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
      }
      else
      {
        xn = (workset.bases[basis_index])->basis_coordinates(cell,point,0);
        yn = (workset.bases[basis_index])->basis_coordinates(cell,point,1);
      }

      const ScalarT& n = edensity(cell,point)*C0;
      const ScalarT& p = hdensity(cell,point)*C0;
      const ScalarT& ni = intrin_conc(cell,point)*C0;

      interpolator->interpolateToPoint(xn, yn, zn, curr_time, sinkLifetime);

      defect_cluster_rate(cell,point) = sinkLifetime * cascadeDensity *
        (p * n - ni*ni)/(p + n + 2.0*ni) * scale_const;

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
RecombRate_Defect_Cluster<EvalT, Traits>::getValidParameters() const
{
  using charon::Vector;

  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<std::string>("Interpolant Method", "ONED_LINEAR_INTERPOLATION");
  p->set<std::string>("Input File Type", "CLUSTER1DX");
  p->set<double>("Cascade Density", 0.0);
  p->set<int>("Number of Input Files", 0);
  p->set<bool>("Is IP Set", false);
  p->set<double>("Shepard Power", 2.0);
  p->set<double>("Influence Radius",-1.0);


  Teuchos::RCP<charon::clusterInterpolator> cInterp;
  p->set<Teuchos::RCP<charon::clusterInterpolator> >("cluster interpolator", cInterp);
  p->set<std::string>("Cluster Netlist Template","");

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);




  return p;
}

  //For in situ clustering

}

#endif //CHARON_RECOMBRATE_DEFECT_CLUSTER_IMPL_HPP
