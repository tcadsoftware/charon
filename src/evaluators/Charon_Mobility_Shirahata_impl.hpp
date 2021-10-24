
#ifndef CHARON_MOBILITY_SHIRAHATA_IMPL_HPP
#define CHARON_MOBILITY_SHIRAHATA_IMPL_HPP

#include <cmath>
#include <algorithm>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_ViewFactory.hpp"

#include <iostream>

#include "Charon_Material_Properties.hpp"
#include "Charon_Names.hpp"

/*
Shirahata mobility model -- reference by Shirahata, Hauser, and Roulston,
"Electron and Hole Mobilities in Silicon as a Function of
Concentration and Temperature," IEEE Transactions on Electron
Devices, Vol. ED-29, pp.292-295, 1982.

This mobility model is a low-field model, contains effects from phonon scattering,
ionized impurity, and carrier-carrier scatterings. It does not differ majority from
minority mobility and does not include screening. The electron mobility expression
is good for Phosphorous-doped Silicon, 250-500 K and 1e13-1e20 cm^-3 doping.
The hole expresson is good for Boron-doped Silicon, 200-400 K and 1e15-1e20 cm^-3.
Both expressions have about 13% error w.r.t. reported experimental values.
*/

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Mobility_Shirahata<EvalT, Traits>::
Mobility_Shirahata(
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


  criticalDistance = 1e12;
  if(p.sublist("Mobility ParameterList").isParameter("Mobility Boundary Layer"))
    criticalDistance = p.sublist("Mobility ParameterList").get<double>("Mobility Boundary Layer");

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  isEdgedl = p.get<bool>("Is Edge Data Layout");
  isNodalDL = p.get<bool>("Nodal Evaluator");

  RCP<DataLayout> output_dl;
  RCP<DataLayout> vector, scalar;
  RCP<BasisIRLayout> basis;
  RCP<IntegrationRule> ir;

  if(isNodalDL)
    {
      // BASIS
      basis = p.get<RCP<BasisIRLayout> >("Basis");
      scalar = basis->functional;
      num_points = scalar->dimension(1);
      basis_name = basis->name();
      output_dl = scalar;

      //just to get dimension
      ir = p.get< RCP<IntegrationRule> >("IR");
      num_dims = ir->dl_vector->dimension(2);

      if(isEdgedl)
	{
	  RCP<const panzer::CellTopologyInfo> cellTopoInfo = basis->getCellTopologyInfo();
	  output_dl = cellTopoInfo->edge_scalar;
	  num_edges = output_dl->dimension(1);
	  cellType = cellTopoInfo->getCellTopology();
	}
    }
  else
    {
      // IP
      ir = p.get< RCP<IntegrationRule> >("IR");
      scalar = ir->dl_scalar;
      vector = ir->dl_vector;
      output_dl = ir->dl_scalar;
      num_points = scalar->dimension(1);
      num_dims = vector->dimension(2);
      int_rule_degree = ir->cubature_degree;
    }

  // Obtain material name
  const string& matName = p.get<string>("Material Name");

  // Obtain carrier type
  carrType = p.get<string>("Carrier Type");

  // Mobility ParameterList
  const ParameterList& mobParamList = p.sublist("Mobility ParameterList");

  // Initialize Shirahata mobility parameters, must be called here since
  // it uses the above variables in the code
  initMobilityParams(matName, mobParamList);

  // Evaluated field
  if (carrType == "Electron")
    {
      mobility = MDField<ScalarT,Cell,Point>(n.field.elec_shirahata_mobility, output_dl);
      potentialSign =  1.0;
    }
  else if (carrType == "Hole")
    {
      mobility = MDField<ScalarT,Cell,Point>(n.field.hole_shirahata_mobility, output_dl);
      potentialSign = -1.0;
    }
  this->addEvaluatedField(mobility);

  // Scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  Mu0 = scaleParams->scale_params.Mu0;
  C0 = scaleParams->scale_params.C0;
  T0 = scaleParams->scale_params.T0;
  E0 = scaleParams->scale_params.E0;
  X0 = scaleParams->scale_params.X0;

  // Depedent fields
  latt_temp = MDField<const ScalarT,Cell,Point>(n.field.latt_temp,scalar);
  this->addDependentField(latt_temp);

  if(!isNodalDL)
    {
      // SUPG-FEM
      if (carrType == "Electron")
	eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.elec_efield, vector);
      else if (carrType == "Hole")
	eff_field = MDField<const ScalarT,Cell,Point,Dim>(n.field.hole_efield, vector);
      this->addDependentField(eff_field);
    }
  else
    {
      // Always need these fields when isEdgedl = true (for CVFEM-SG and EFFPG-FEM)
      intrin_fermi = MDField<const ScalarT,Cell,Point>(n.field.intrin_fermi, scalar);
      bandgap = MDField<const ScalarT,Cell,Point>(n.field.band_gap, scalar);
      eff_bandgap = MDField<const ScalarT,Cell,Point>(n.field.eff_band_gap, scalar);

      this->addDependentField(intrin_fermi);
      this->addDependentField(bandgap);
      this->addDependentField(eff_bandgap);
    }

  //Calculate a scale factor based on distance from the oxide interface also convert um to cm
  if(p.sublist("Mobility ParameterList").isParameter("Start X"))
    xOIStart = p.sublist("Mobility ParameterList").get<double>("Start X");
  if(p.sublist("Mobility ParameterList").isParameter("Start Y"))
    yOIStart = p.sublist("Mobility ParameterList").get<double>("Start Y");
  if(p.sublist("Mobility ParameterList").isParameter("Start Z"))
    zOIStart = p.sublist("Mobility ParameterList").get<double>("Start Z");
  if(p.sublist("Mobility ParameterList").isParameter("End X"))
    xOIEnd = p.sublist("Mobility ParameterList").get<double>("End X");
  if(p.sublist("Mobility ParameterList").isParameter("End Y"))
    yOIEnd = p.sublist("Mobility ParameterList").get<double>("End Y");
  if(p.sublist("Mobility ParameterList").isParameter("End Z"))
    zOIEnd = p.sublist("Mobility ParameterList").get<double>("End Z");

  if(num_dims == 2)
    {
      double normx = -(yOIEnd - yOIStart);
      double normy = xOIEnd - xOIStart;
      double normMag = sqrt(normx*normx + normy*normy);
      oxideNorm.push_back(normx/normMag);
      oxideNorm.push_back(normy/normMag);
      oxideNorm.push_back(0.0);
    }
  else if (num_dims == 3)
    {
      // P1, P2, P3 points in the plane
      double xp1=xOIStart, yp1=yOIStart, zp1=zOIStart;
      double xp2=xOIEnd, yp2=yOIEnd, zp2=zOIStart;
      double xp3=xOIStart, yp3=yOIStart, zp3=zOIEnd;
      // P1P2, P1P3 vectors in the plane
      double x12 = xp2-xp1, y12 = yp2-yp1, z12 = zp2-zp1;
      double x13 = xp3-xp1, y13 = yp3-yp1, z13 = zp3-zp1;
      // find the normal to the plane (P1P2, P1P3) defined by P1P2 x P1P3
      double normx = y12*z13-y13*z12;
      double normy = x13*z12-x12*z13;
      double normz = x12*y13-x13*y12;
      double normMag = sqrt(normx*normx + normy*normy + normz*normz);
      oxideNorm.push_back(normx/normMag);
      oxideNorm.push_back(normy/normMag);
      oxideNorm.push_back(normz/normMag);
    }


  std::string name = "Shirahata_Mobility_Model";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Shirahata<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if (isNodalDL)
    basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);
  else
    int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
Mobility_Shirahata<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{
  using panzer::index_t;

  // Reference temperature [K]
  ScalarT refTemp = 300.0;

  // FieldContainer to temporarily hold mobility values
  Kokkos::DynRankView<ScalarT,PHX::Device> tmpMob = 
    Kokkos::createDynRankView(mobility.get_static_view(),"tmpMob",workset.num_cells, num_points);

  // Compute the mobility
  if(!isEdgedl) 
    {
      for (index_t cell = 0; cell < workset.num_cells; ++cell)
	{
	  for (int point = 0; point < num_points; ++point)
	    {
	      // Obtain temperature [K]
	      ScalarT lattTemp = latt_temp(cell,point) * T0;
	      
	      // Compute temperature-dependent parameters
	      ScalarT normTemp = lattTemp / refTemp;
	      
	      ScalarT Field = 0.0;
	      for (int dim = 0; dim < num_dims; ++dim)
		{
		  Field += oxideNorm[dim] * eff_field(cell,point,dim) * E0;
		}
	      Field = fabs(Field);
	      
	      ScalarT mob1 = pow(1.0+Field/E1,P1) + pow(Field/E2,P2);
	      tmpMob(cell,point) = muo * pow(normTemp,-theta) / mob1 / Mu0;
	    }
	}
    }

  double xn,yn,zn,xn0,xn1,yn0,yn1,zn0,zn1,distance;

  // Want mobility available at the center of primary edges for CVFEM-SG and EFFPG-FEM
  if (isEdgedl)
  {
    // compute mobility at the center of primary edge
    for (index_t cell = 0; cell < workset.num_cells; ++cell)
    {
      for (int edge = 0; edge < num_edges; ++edge)
      {
        // get local node ids: first index 1 for edge (0 for vertex, 2 for face, 3 for volume)
        int node0 = cellType->getNodeMap(1,edge,0);
        int node1 = cellType->getNodeMap(1,edge,1);

	xn0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,0);
	xn1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,0);
	yn0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,1);
	yn1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,1);
	zn0 = 0.0;
	zn1 = 0.0;


	if(num_dims == 3)
	  {
	    zn0 = (workset.bases[basis_index])->basis_coordinates(cell,node0,2);
	    zn1 = (workset.bases[basis_index])->basis_coordinates(cell,node1,2);
	  }

	xn = 0.5*( xn0 + xn1);
	yn = 0.5*( yn0 + yn1);
	zn = 0.0;
	if(num_dims == 3)
	  zn = 0.5*( zn0 + zn1);

	double edgeLength = sqrt((xn1-xn0)*(xn1-xn0) + (yn1-yn0)*(yn1-yn0) + (zn1-zn0)*(zn1-zn0));

	// compute the effective potential at local nodes
	ScalarT dEg0 = bandgap(cell,node0) - eff_bandgap(cell,node0) ;  // [eV]
	ScalarT dEg1 = bandgap(cell,node1) - eff_bandgap(cell,node1) ;
	ScalarT iEf0 = intrin_fermi(cell,node0);  // [eV]
	ScalarT iEf1 = intrin_fermi(cell,node1);

	ScalarT effPot0 = (potentialSign*0.5*dEg0 - iEf0); // 1.0 converts [eV] to [V]
	ScalarT effPot1 = (potentialSign*0.5*dEg1 - iEf1);

	// primary edge effective electric field in [V/cm], Eij = -(Potj-Poti)/hij
	ScalarT edgeField = -(effPot1-effPot0) / (edgeLength*X0);

	// This field is directed along the edge.  For perp field, dot edge vector 
        // with interface norm
	double xvec = (xn1-xn0)/edgeLength;
	double yvec = (yn1-yn0)/edgeLength;
	double zvec = 0.0;
	if(num_dims == 3)
	  zvec = (zn1-zn0)/edgeLength;

	double FieldScaleFactor = xvec*oxideNorm[0] + yvec*oxideNorm[1] + zvec*oxideNorm[2];
	ScalarT perpField = fabs(FieldScaleFactor * edgeField);

	// Obtain temperature [K]
	ScalarT lattTemp0 = latt_temp(cell,node0) * T0;
	ScalarT lattTemp1 = latt_temp(cell,node1) * T0;

	// Compute temperature-dependent parameters
	ScalarT normTemp0 = lattTemp0 / refTemp;
	ScalarT normTemp1 = lattTemp1 / refTemp;
	ScalarT normMidTemp = 0.5 * (normTemp0 + normTemp1);

	ScalarT mob1 = pow(1.0+perpField/E1,P1) + pow(perpField/E2,P2);
	ScalarT mob  = muo * pow(normMidTemp,-theta) / mob1;

        /*
	if(xn < xOIStart || xn > xOIEnd)
	  distance = 1;
	else
	  distance = yOIStart - yn;
	
	double distanceFactor = exp(distance/criticalDistance);
	if(distanceFactor < 1e-16)
	  distanceFactor = 1e-16;

	if(distance > 2*criticalDistance)
	  mob = 0.0;

        mobility(cell,edge) = mob/distanceFactor/Mu0;
        */
    
	distance = yOIStart - yn;
	double distanceFactor = std::exp(-distance/criticalDistance);
	
	double latdistanceFactor = 1.;
	double latcriticalDistance = 0.01; // 10 nm
	if(xn < xOIStart ||  xn > xOIEnd) {
	  double lat_distance = (xn < xOIStart) ? xOIStart - xn : xn - xOIEnd;
	  latdistanceFactor = std::exp(-lat_distance/latcriticalDistance);
	}

       if (distance > 2*criticalDistance)
         mob = 0.0;

	mobility(cell,edge) = mob*latdistanceFactor*distanceFactor/Mu0;
      }
    }
  }  // end of if block

  // Want mobility available at IP or BASIS
  else
  {
   for (index_t cell = 0; cell < workset.num_cells; ++cell)
     for (int point = 0; point < num_points; ++point)
       {
	 //Calculate a scale factor based on distance from the oxide interface
	 if (isNodalDL)
	   {
	     xn = (workset.bases[basis_index])->basis_coordinates(cell,point,0);
	     yn = (workset.bases[basis_index])->basis_coordinates(cell,point,1);
	     if(num_dims == 3)
	       zn = (workset.bases[basis_index])->basis_coordinates(cell,point,2);
	   }
	 else
	   {
	     xn = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,0);
	     yn = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,1);
	     zn = 0.0;
	     if(num_dims == 3)
	       zn = (workset.int_rules[int_rule_index])->ip_coordinates(cell,point,2);
	   }

	 zn = zn;  //This exists for the sole purpose of squelching a warning for the time being.

         /*
	 if(xn < xOIStart || xn > xOIEnd)
	   distance = 1;
	 else
	   distance = yOIStart - yn;

	 double distanceFactor = exp(distance/criticalDistance);
	 if(distanceFactor < 1e-16)
	   distanceFactor = 1e-16;

	 ScalarT mob = tmpMob(cell,point);

	 if(distance > criticalDistance)
	   mob = 0.0;

	 mobility(cell,point) = mob/distanceFactor;
	 */
       
         distance = yOIStart - yn;
	 double distanceFactor = std::exp(-distance/criticalDistance);

	 double latdistanceFactor = 1;
	 double latcriticalDistance = 0.01; // 10 nm
	 if(xn < xOIStart || xn > xOIEnd) {
	   double lat_distance = (xn < xOIStart) ? xOIStart - xn : xn - xOIEnd;
	   latdistanceFactor = std::exp(-lat_distance/latcriticalDistance);
	 }
	 
	 ScalarT mob = tmpMob(cell,point);
	 mobility(cell,point) = latdistanceFactor*distanceFactor*mob;
       }
  }

}


///////////////////////////////////////////////////////////////////////////////
//
//  initMobilityParams()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void Mobility_Shirahata<EvalT, Traits>::initMobilityParams
(const std::string& matName, const Teuchos::ParameterList& mobParamList)
{
  using std::string;

  // Obtain the instance of charon::Material_Properties
  charon::Material_Properties& matProperty = charon::Material_Properties::getInstance();

  // Set up parameters for Shirahata electron mobility model
  if (carrType == "Electron")
  {
    // Retrieve parameters from charon::Material_Properties by default
    muo  = matProperty.getPropertyValue(matName, "Shirahata Electron muo");
    theta = matProperty.getPropertyValue(matName, "Shirahata Electron theta");
    E1    = matProperty.getPropertyValue(matName, "Shirahata Electron E1");
    E2    = matProperty.getPropertyValue(matName, "Shirahata Electron E2");
    P1    = matProperty.getPropertyValue(matName, "Shirahata Electron P1");
    P2    = matProperty.getPropertyValue(matName, "Shirahata Electron P2");
  }

  // Set up parameters for Shirahata hole mobility model
  else if (carrType == "Hole")
  {
    // Retrieve parameters from charon::Material_Properties by default
    muo  = matProperty.getPropertyValue(matName, "Shirahata Hole muo");
    theta = matProperty.getPropertyValue(matName, "Shirahata Hole theta");
    E1    = matProperty.getPropertyValue(matName, "Shirahata Hole E1");
    E2    = matProperty.getPropertyValue(matName, "Shirahata Hole E2");
    P1    = matProperty.getPropertyValue(matName, "Shirahata Hole P1");
    P2    = matProperty.getPropertyValue(matName, "Shirahata Hole P2");
  }

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
      << "Invalid Carrier Type ! Must be either Electron or Hole !");


  // Overwrite parameters when specified by users
  if (mobParamList.isParameter("muo"))
    muo = mobParamList.get<double>("muo");
  if (mobParamList.isParameter("theta"))
    theta = mobParamList.get<double>("theta");

}


///////////////////////////////////////////////////////////////////////////////
//
//  getValidParameters()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList>
Mobility_Shirahata<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Material Name", "?");
  p->set("Carrier Type", "?");

  Teuchos::RCP<const charon::Names> n;
  p->set("Names", n);

  Teuchos::RCP<PHX::DataLayout> dl;
  p->set("Data Layout", dl);

  p->sublist("Mobility ParameterList").set("Start X",0.0);
  p->sublist("Mobility ParameterList").set("Start Y",0.0);
  p->sublist("Mobility ParameterList").set("Start Z",0.0);
  p->sublist("Mobility ParameterList").set("End X",0.0);
  p->sublist("Mobility ParameterList").set("End Y",0.0);
  p->sublist("Mobility ParameterList").set("End Z",0.0);
  p->sublist("Mobility ParameterList").set("Mobility Boundary Layer",1.0e12);

  p->sublist("Mobility ParameterList", false, "");
  p->sublist("Mobility ParameterList").set<std::string>("Value", "Shirahata", "Shirahata low-field mobility model");
  p->sublist("Mobility ParameterList").set<double>("muo", 0., "[cm^2/(V.s)]");
  p->sublist("Mobility ParameterList").set<double>("theta", 0., "");
  p->sublist("Mobility ParameterList").set<double>("E1", 0., "");
  p->sublist("Mobility ParameterList").set<double>("E2", 0., "");
  p->sublist("Mobility ParameterList").set<double>("P1", 0., "");
  p->sublist("Mobility ParameterList").set<double>("P2", 0., "");

  //p->sublist("Mobility ParameterList").set<std::string>("Driving Force", "ElectricField", "Electric Field");

  p->set("Is Edge Data Layout", false);
  p->set("Nodal Evaluator", false);

  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);

  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);

  Teuchos::RCP<charon::Scaling_Parameters> sp;
  p->set("Scaling Parameters", sp);

  return p;
}

}

#endif
