
#ifndef CHARON_DOPING_STEPJUNCTION_IMPL_HPP
#define CHARON_DOPING_STEPJUNCTION_IMPL_HPP

#include <cmath>
#include "Teuchos_Assert.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Charon_Names.hpp"


namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Doping_StepJunction<EvalT, Traits>::
Doping_StepJunction(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  // IP
  RCP<IntegrationRule> ir = p.get< RCP<IntegrationRule> >("IR");
  RCP<DataLayout> scalar = ir->dl_scalar;
  RCP<DataLayout> vector = ir->dl_vector;

  int_rule_degree = ir->cubature_degree;
  num_ip = vector->dimension(1);
  num_dim = vector->dimension(2);

  // basis
  RCP<BasisIRLayout> basis = p.get<RCP<BasisIRLayout> >("Basis");
  RCP<DataLayout> data_layout = basis->functional;
  basis_name = basis->name();

  // input from user
  acceptorValue = p.get<double>("Acceptor Value");
  donorValue = p.get<double>("Donor Value");
  junctionLoc = p.get<double>("Junction Location");
  config = p.get<std::string>("Configuration");
  direction = p.get<std::string>("Direction");

  // fields
  doping_raw = MDField<ScalarT,Cell,IP>(n.field.doping_raw,scalar);
  acceptor_raw = MDField<ScalarT,Cell,IP>(n.field.acceptor_raw,scalar);
  donor_raw = MDField<ScalarT,Cell,IP>(n.field.donor_raw,scalar);

  doping_raw_basis = MDField<ScalarT,Cell,BASIS>(n.field.doping_raw,data_layout);
  acceptor_raw_basis = MDField<ScalarT,Cell,BASIS>(n.field.acceptor_raw,data_layout);
  donor_raw_basis = MDField<ScalarT,Cell,BASIS>(n.field.donor_raw,data_layout);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  C0 = scaleParams->scale_params.C0;

  this->addEvaluatedField(doping_raw);
  this->addEvaluatedField(acceptor_raw);
  this->addEvaluatedField(donor_raw);

  this->addEvaluatedField(doping_raw_basis);
  this->addEvaluatedField(acceptor_raw_basis);
  this->addEvaluatedField(donor_raw_basis);

  std::string name = "Doping_StepJunction";
  this->setName(name);
}

///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Doping_StepJunction<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Doping_StepJunction<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;
  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  size_type num_basis = doping_raw_basis.dimension(1);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    // doping at IPs
    for (int ip = 0; ip < num_ip; ++ip)
    {
      double x = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,0);
      double y = 0.0, z = 0.0;
      if (num_dim == 2)
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
      if (num_dim == 3)
      {
        y = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
        z = (workset.int_rules[int_rule_index])->ip_coordinates(cell,ip,2);
      }

      std::vector<double> dopValue = evaluateDoping(x,y,z);

      acceptor_raw(cell,ip) = dopValue[0]/C0;
      donor_raw(cell,ip) = dopValue[1]/C0;
      doping_raw(cell,ip) = (dopValue[1] - dopValue[0])/C0;
    }

    // doping at basis points
    for (size_type basis = 0; basis < num_basis; ++basis)
    {
      double x = (workset.bases[basis_index])->basis_coordinates(cell,basis,0);
      double y = 0.0, z = 0.0;
      if (num_dim == 2)
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
      if (num_dim == 3)
      {
        y = (workset.bases[basis_index])->basis_coordinates(cell,basis,1);
        z = (workset.bases[basis_index])->basis_coordinates(cell,basis,2);
      }

      std::vector<double> dopValue = evaluateDoping(x,y,z);

      acceptor_raw_basis(cell,basis) = dopValue[0]/C0;
      donor_raw_basis(cell,basis) = dopValue[1]/C0;
      doping_raw_basis(cell,basis) = (dopValue[1] - dopValue[0])/C0;
    }
  } // end of loop over cells
}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Doping_StepJunction<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  p->set<double>("Acceptor Value", 0.0);
  p->set<double>("Donor Value", 0.0);
  p->set<double>("Junction Location", 0.0);
  p->set<std::string>("Configuration", "??");
  p->set<std::string>("Direction", "??");
  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> Doping_StepJunction<EvalT, Traits>::
evaluateDoping(const double& x, const double& y, const double& z)
{
  std::vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  // PN junction
  if (config == "PN")
  {
    // PN junction along x direction
    if (direction == "X")
    {
      if (x < junctionLoc)
      {
        dopValue[0] = acceptorValue;  // hold acceptor value
        dopValue[1] = 0.0;            // hold donor value
      }
      else if (x > junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    // PN junction along y direction
    else if (direction == "Y")
    {
      if (y < junctionLoc)
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = 0.0;
      }
      else if (y > junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    // PN junction along z direction
    else if (direction == "Z")
    {
      if (z < junctionLoc)
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = 0.0;
      }
      else if (z > junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Invalid step junction direction and it has to be X, or Y, or Z!");

  }  // end of if (config == "PN")


  // NP junction
  else if (config == "NP")
  {
    // NP junction along x direction
    if (direction == "X")
    {
      if (x < junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else if (x > junctionLoc)
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = 0.0;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    // NP junction along y direction
    else if (direction == "Y")
    {
      if (y < junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else if ( y > junctionLoc)
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = 0.0;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    // NP junction along z direction
    else if (direction == "Z")
    {
      if (z < junctionLoc)
      {
        dopValue[0] = 0.0;
        dopValue[1] = donorValue;
      }
      else if (z > junctionLoc)
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = 0.0;
      }
      else  // at the junction
      {
        dopValue[0] = acceptorValue;
        dopValue[1] = donorValue;
      }
    }

    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Invalid step junction direction and it has to be X, or Y, or Z!");

  }  // end of else if (config == "NP")

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Invalid step junction configuration and it has to be either PN or NP !");

  return dopValue;
}

}

#endif
