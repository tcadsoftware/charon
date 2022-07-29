
#ifndef CHARON_DOPINGRAW_FUNCTION_IMPL_HPP
#define CHARON_DOPINGRAW_FUNCTION_IMPL_HPP

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include "Charon_Names.hpp"
#include "lusolve.hpp"

namespace charon {


///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
DopingRaw_Function<EvalT, Traits>::
DopingRaw_Function(
  const Teuchos::ParameterList& p)
{
  using std::string;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  const charon::Names& n =
    *(p.get< Teuchos::RCP<const charon::Names> >("Names"));

  sweepingIsOn = false;
  if (p.isParameter("SweepingIsOn"))
    sweepingIsOn = p.get<bool>("SweepingIsOn");

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

  // doping parameterlist
  dopParamList = p.sublist("Doping ParameterList");

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

  //homotopy parameters

  user_value = rcp(new panzer::ScalarParameterEntry<EvalT>);
  user_value->setRealValue(1.0);
  if(dopParamList.isType<double>("Doping Scaling"))
    {
      if((dopParamList.isType<std::string>("Doping Homotopy")            )  &&
         (dopParamList.get<std::string>("Doping Homotopy") == "Parameter"))
        {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! Cannot specify both a scale factor and a homotopy in the doping block.\n\n");
        }
      else
        {
          user_value->setRealValue(dopParamList.get<double>("Doping Scaling"));
        }
    }

  if ((dopParamList.isType<std::string>("Doping Homotopy")            )  and
      (dopParamList.get<std::string>("Doping Homotopy") == "Parameter"))
    {
      user_value =
        panzer::createAndRegisterScalarParameter<EvalT>(
                                                        std::string("Doping Homotopy"),
                                                        *dopParamList.get<RCP<panzer::ParamLib> >("ParamLib"));
    }


  this->addEvaluatedField(doping_raw);
  this->addEvaluatedField(acceptor_raw);
  this->addEvaluatedField(donor_raw);

  this->addEvaluatedField(doping_raw_basis);
  this->addEvaluatedField(acceptor_raw_basis);
  this->addEvaluatedField(donor_raw_basis);

  store_wkst_doping = false;
  std::string name = "DopingRaw_Function";
  this->setName(name);

  // define the outermost vector size for gauss decay params vectors
  decayPos.resize(2);  // 0 for File1D, 1 for Uniform
  decayDir.resize(2);
  decayWidth.resize(2);
  int udCounter = 0; 

  for (ParameterList::ConstIterator model_it = dopParamList.begin();
       model_it !=dopParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");
       
      if (funcType == "File2D")
      {
        DV.resize(DV.size()+1);
        readDopingFile2D(funcParamList);
        store_wkst_doping = true;
      }
      else if (funcType == "File3D")
      {
        DV.resize(DV.size()+1);
        readDopingFile3D(funcParamList);
        store_wkst_doping = true;
      }
      else if (funcType == "File1D") 
      {
        DV1.resize(DV1.size()+1);
        readDopingFile1D(funcParamList);
        store_wkst_doping = true;
      }
      else if (funcType == "Uniform")
      {
        uniformDopingParams udp_;
        udp_.parseUniform(funcParamList);
        udp_vec.push_back(udp_);
 
	if (udp_.sweepMe)
	  {
	    user_value =
	      panzer::createAndRegisterScalarParameter<EvalT>(
							      std::string("Doping Value"),
							      *dopParamList.get<RCP<panzer::ParamLib> >("ParamLib"));
	  }
	// resize gauss decay vectors for uniform doping
        udCounter += 1; 
        int udIndex = udCounter-1;
        decayDir[1].resize(udCounter);
        decayPos[1].resize(udCounter);
        decayWidth[1].resize(udCounter);
  
        // read in Gauss decay parameters when specified
        if (funcParamList.isSublist("Gauss Decay"))
        {
          const ParameterList& gaussPList = funcParamList.sublist("Gauss Decay");
          readGaussDecayParams(1, udIndex, gaussPList);
        }
        else  // initialize the vectors
        {
          decayDir[1][udIndex].push_back("");
          decayPos[1][udIndex].push_back(0.0);
          decayWidth[1][udIndex].push_back(0.0);
        }
      }
      else if ( (funcType == "Gauss") || (funcType == "Gaussian") )
      {
        gaussianDopingParams gdp_;
        gdp_.parseGaussian(funcParamList,num_dim);
        gdp_vec.push_back(gdp_);
      }
      else if (funcType == "Linear")
      {
        linearDopingParams ldp_;
        ldp_.parseLinear(funcParamList,num_dim);
        ldp_vec.push_back(ldp_);
      }
      else if (funcType == "Erfc")
      {
        erfcDopingParams edp_;
        edp_.parseErfc(funcParamList,num_dim);
        edp_vec.push_back(edp_);
      }
      else if (funcType == "MGauss")
      {
        mgaussDopingParams mgdp_;
        mgdp_.parseMGauss(funcParamList,num_dim);
        mgdp_vec.push_back(mgdp_);
      }
      else if (funcType == "Halo")
	{
	  haloDopingParams hdp_;
	  hdp_.parseHalo(funcParamList,num_dim);
	  hdp_vec.push_back(hdp_);
	}
      else if (funcType == "MMS_RDH_1")
      {
        mms1 = Teuchos::rcp(new charon::MMS_DD_RDH_1_AnalyticFunction(p));
      }
      else if (funcType == "MMS_RDH_2")
      {
        mms2 = Teuchos::rcp(new charon::MMS_DD_RDH_2_AnalyticFunction(p));
      }
      else if (funcType == "mms_nlp_glh_1")
      {
        mms_nlp_func = Teuchos::rcp(new charon::MMS_NLP_GLH_1_AnalyticFunction());
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Invalid doping Function Type!"
          << "Must be Uniform, Erfc, Gauss/Gaussian, MGauss, Halo, "
          << "Linear, File1D, File2D or File3D." << std::endl );

    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop


  if (store_wkst_doping)
  {

     std::size_t num_wksts = p.get<std::size_t>("Max Worksets");
     acceptor_raw_wkst.resize(num_wksts);
     donor_raw_wkst.resize(num_wksts);
     acceptor_raw_basis_wkst.resize(num_wksts);
     donor_raw_basis_wkst.resize(num_wksts);
     for (std::size_t i = 0;i < num_wksts; ++i)
     {
        std::string acc_name = "acceptor_"+std::to_string(i);
        std::string don_name = "donor_"+std::to_string(i);
        acceptor_raw_wkst[i] = MDField<ScalarT,Cell,IP>(acc_name,scalar);
        donor_raw_wkst[i] = MDField<ScalarT,Cell,IP>(don_name,scalar);

        std::string acc_name_basis = "acceptor_b_"+std::to_string(i);
        std::string don_name_basis = "donor_b_"+std::to_string(i);
        acceptor_raw_basis_wkst[i] = MDField<ScalarT,Cell,BASIS>(acc_name_basis,data_layout);
        donor_raw_basis_wkst[i] = MDField<ScalarT,Cell,BASIS>(don_name_basis,data_layout);

        this->addEvaluatedField(acceptor_raw_wkst[i]);
        this->addEvaluatedField(donor_raw_wkst[i]);

        this->addEvaluatedField(acceptor_raw_basis_wkst[i]);
        this->addEvaluatedField(donor_raw_basis_wkst[i]);
     }
   }

  prof_eval = rcp(new ProfileEvals(num_dim));

}

template<typename EvalT, typename Traits>
void DopingRaw_Function<EvalT, Traits>::readDopingFile2D(const Teuchos::ParameterList& plist)
{
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::string;

   //read 2D doping profile from an external file and assign it to FileDoping array

   const string FileName = plist.get<string>("File Name");
   doping_struct r;

   ifstream in(FileName.c_str());

   if (!in)
     TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
     << "Error ! Cannot read doping file '" << FileName << "'" << std::endl);

   int dvIndex = DV.size()-1;
   double xmin = 0.0, xmax = 0.0;
   double ymin = 0.0, ymax = 0.0;

   // Size the extents vectors appropriately
   MinX.resize(DV.size());
   MinY.resize(DV.size());
   MaxX.resize(DV.size());
   MaxY.resize(DV.size());

   while( in >> r.x >> r.y >> r.d )
   {
     if (r.d < 0.0) 
       TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
         << "Error ! The doping value in a doping file cannot be negative: " << r.d << std::endl);
      
     //Set z coordinate to zero
     //This should only impact doping interpolation.  The mesh needn't be set on z=0 plane
     r.z = 0;

     if (DV[dvIndex].size() == 0) //initialize the boundaries
     {
       xmin = r.x;
       xmax = xmin;
       ymin = r.y;
       ymax = ymin;
     }
     else
     {
       if (xmin > r.x) xmin = r.x;
       if (xmax < r.x) xmax = r.x;
       if (ymin > r.y) ymin = r.y;
       if (ymax < r.y) ymax = r.y;
     }

     DV[dvIndex].push_back(r);
   }
   
   // save the min/max coordinates
   MinX[dvIndex] = xmin;
   MaxX[dvIndex] = xmax;
   MinY[dvIndex] = ymin;
   MaxY[dvIndex] = ymax;

   // sort the doping points by x and y axis
   std::sort(DV[dvIndex].begin(),DV[dvIndex].end());
   
   // eliminate duplicate (x,y) entries
   DV[dvIndex].resize( unique (DV[dvIndex].begin(), DV[dvIndex].end()) - DV[dvIndex].begin() );
}



template<typename EvalT, typename Traits>
void DopingRaw_Function<EvalT, Traits>::readDopingFile3D(const Teuchos::ParameterList& plist)
{
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::string;

  int me;
  MPI_Comm_rank(MPI_COMM_WORLD,&me);


   //read 3D doping profile from an external file and assign it to FileDoping array
   const string FileName = plist.get<string>("File Name");
   doping_struct r;

   ifstream in(FileName.c_str());

   if (!in)
     TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
     << "Error ! Cannot read doping file '" << FileName << "'" << std::endl);

   int dvIndex = DV.size()-1;
   double xmin = 0.0, xmax = 0.0;
   double ymin = 0.0, ymax = 0.0; 
   double zmin = 0.0, zmax = 0.0; 

   // Size the extents vectors appropriately
   MinX.resize(DV.size());
   MinY.resize(DV.size());
   MinZ.resize(DV.size());
   MaxX.resize(DV.size());
   MaxY.resize(DV.size());
   MaxZ.resize(DV.size());

   while( in >> r.x >> r.y >> r.z >> r.d )
   {
     if (r.d < 0.0) 
       TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
         << "Error ! The doping value in a doping file cannot be negative: " << r.d << std::endl);

     if (DV[dvIndex].size() == 0) //initialize the boundaries
     {
       xmin = r.x;
       xmax = xmin;
       ymin = r.y;
       ymax = ymin;
       zmin = r.z;
       zmax = zmin; 
     }
     else
     {
       if (xmin > r.x) xmin = r.x;
       if (xmax < r.x) xmax = r.x;
       if (ymin > r.y) ymin = r.y;
       if (ymax < r.y) ymax = r.y;
       if (zmin > r.y) zmin = r.z;
       if (zmax < r.y) zmax = r.z;
     }
 
     DV[dvIndex].push_back(r);
   }

   // save the min/max coordinates
   MinX[dvIndex] = xmin;
   MaxX[dvIndex] = xmax;
   MinY[dvIndex] = ymin;
   MaxY[dvIndex] = ymax;
   MinZ[dvIndex] = zmin;
   MaxZ[dvIndex] = zmax;

   // sort the doping points by x and y axis
   sort(DV[dvIndex].begin(),DV[dvIndex].end());

  // eliminate duplicate (x,y) entries
   DV[dvIndex].resize( unique (DV[dvIndex].begin(), DV[dvIndex].end()) - DV[dvIndex].begin() );

}


template<typename EvalT, typename Traits>
void DopingRaw_Function<EvalT, Traits>::readDopingFile1D(const Teuchos::ParameterList& plist)
{
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::string;
  using Teuchos::ParameterList; 

  //read 1D doping profile from an external file and assign it to FileDoping array
  const string FileName = plist.get<string>("File Name");
  TEUCHOS_ASSERT(FileName.size() > 0);
   
  doping_struct_1D r;

  ifstream in(FileName.c_str());

  if (!in)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Cannot read doping file '" << FileName << "'" << std::endl);

  int dvIndex = DV1.size()-1;

  while( in >> r.x >> r.d ) 
  {
    if (r.d < 0.0) 
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Error ! The doping value in a doping file cannot be negative: " << r.d << std::endl);
    DV1[dvIndex].push_back(r);
  }
   
  // sort the doping points by x
  sort(DV1[dvIndex].begin(),DV1[dvIndex].end());

  // eliminate duplicate (x,d) entries
  DV1[dvIndex].resize( unique (DV1[dvIndex].begin(), DV1[dvIndex].end()) - DV1[dvIndex].begin() );
  
  // resize the vectors
  decayDir[0].resize(DV1.size());
  decayPos[0].resize(DV1.size());
  decayWidth[0].resize(DV1.size());
  
  // read in Gauss decay parameters when specified
  if (plist.isSublist("Gauss Decay"))
  {
    const ParameterList& gaussPList = plist.sublist("Gauss Decay"); 
    readGaussDecayParams(0, dvIndex, gaussPList);
  }
  else  // initialize the vectors
  {
    decayDir[0][dvIndex].push_back(""); 
    decayPos[0][dvIndex].push_back(0.0);
    decayWidth[0][dvIndex].push_back(0.0);
  }
}


///////////////////////////////////////////////////////////////////////////////
//
//  postRegistrationSetup()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DopingRaw_Function<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
  basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

  
  if (store_wkst_doping)
  {
      TEUCHOS_TEST_FOR_EXCEPTION((acceptor_raw_wkst.size() < (*sd.worksets_).size()),std::logic_error,
                             "DopingRaw: Workset fields for storage too small.\n");

     int_rule_index = panzer::getIntegrationRuleIndex(int_rule_degree,(*sd.worksets_)[0]);
     basis_index = panzer::getBasisIndex(basis_name,(*sd.worksets_)[0]);

     using panzer::index_t;
     typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
     size_type num_basis = doping_raw_basis.dimension(1);

     int worksetnum = 0;
     for (const auto& wkst : (*sd.worksets_))
     {
        for (index_t cell = 0; cell < wkst.num_cells; ++cell)
        {
          // doping at IPs
          for (int ip = 0; ip < num_ip; ++ip)
          {
             double x = (wkst.int_rules[int_rule_index])->ip_coordinates(cell,ip,0);
             double y = 0.0, z = 0.0;
             if (num_dim == 2)
                y = (wkst.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
             if (num_dim == 3)
             {
                y = (wkst.int_rules[int_rule_index])->ip_coordinates(cell,ip,1);
                z = (wkst.int_rules[int_rule_index])->ip_coordinates(cell,ip,2);
             }

             // evaluate the acceptor and donor doping
             std::vector<double> dopValue = evaluateDoping(x,y,z);

             acceptor_raw_wkst[worksetnum](cell,ip) = dopValue[0]/C0;
             donor_raw_wkst[worksetnum](cell,ip) = dopValue[1]/C0;
          }

          // doping at basis points
          for (size_type basis = 0; basis < num_basis; ++basis)
          {
             double x = (wkst.bases[basis_index])->basis_coordinates(cell,basis,0);
             double y = 0.0, z = 0.0;
             if (num_dim == 2)
               y = (wkst.bases[basis_index])->basis_coordinates(cell,basis,1);
             if (num_dim == 3)
             {
               y = (wkst.bases[basis_index])->basis_coordinates(cell,basis,1);
               z = (wkst.bases[basis_index])->basis_coordinates(cell,basis,2);
             }

             std::vector<double> dopValue = evaluateDoping(x,y,z);

             acceptor_raw_basis_wkst[worksetnum](cell,basis) = dopValue[0]/C0;
             donor_raw_basis_wkst[worksetnum](cell,basis) = dopValue[1]/C0;
          }

        } // end of loop over cells

        worksetnum++;
     } // end of loop over worksets
  } // end if

}

///////////////////////////////////////////////////////////////////////////////
//
//  preEvaluate()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT,typename Traits>
void DopingRaw_Function<EvalT,Traits>::preEvaluate(typename Traits::PreEvalData /* d */)
{
  // initialize worksetId
  worksetId = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void
DopingRaw_Function<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ using std::cout;
  using std::endl;
  using std::ofstream;
  using std::ios;
  using panzer::index_t;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef typename PHX::MDField<ScalarT,Cell,BASIS>::size_type size_type;
  size_type num_basis = doping_raw_basis.dimension(1);

  double dopeHomotopy = 1.0;
  if (not sweepingIsOn)
    user_value->getRealValue();

  if (!store_wkst_doping)
  {
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

          // evaluate the acceptor and donor doping
          std::vector<double> dopValue = evaluateDoping(x,y,z);

          acceptor_raw(cell,ip) = dopeHomotopy*dopValue[0]/C0;
          donor_raw(cell,ip) = dopeHomotopy*dopValue[1]/C0;
          doping_raw(cell,ip) = dopeHomotopy*(dopValue[1] - dopValue[0])/C0;
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
          acceptor_raw_basis(cell,basis) = dopeHomotopy*dopValue[0]/C0;
          donor_raw_basis(cell,basis) = dopeHomotopy*dopValue[1]/C0;
          doping_raw_basis(cell,basis) = dopeHomotopy*(dopValue[1] - dopValue[0])/C0;
        }
     } // end of loop over cells
   } // end if

   else
   {
     for (index_t cell = 0; cell < workset.num_cells; ++cell)
     {

       // doping at IPs
       for (int ip = 0; ip < num_ip; ++ip)
       {
          acceptor_raw(cell,ip) = dopeHomotopy*acceptor_raw_wkst[worksetId](cell,ip);
          donor_raw(cell,ip) = dopeHomotopy*donor_raw_wkst[worksetId](cell,ip);
          doping_raw(cell,ip) = dopeHomotopy*(donor_raw(cell,ip) - acceptor_raw(cell,ip));
       }

       // doping at basis
       for (size_type basis = 0; basis < num_basis; ++basis)
       {
          acceptor_raw_basis(cell,basis) = dopeHomotopy*acceptor_raw_basis_wkst[worksetId](cell,basis);
          donor_raw_basis(cell,basis) = dopeHomotopy*donor_raw_basis_wkst[worksetId](cell,basis);
          doping_raw_basis(cell,basis) = dopeHomotopy*(donor_raw_basis(cell,basis) - acceptor_raw_basis(cell,basis));
       }

     }
     ++worksetId;
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evaluateDoping(
  const double& x, const double& y, const double& z)
{
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  std::vector<double> tempVal(2, 0.0);
  /*
  for (std::size_t i=0; i < udp_vec.size(); ++i)
  {
    tempVal=evalUniformDoping(x,y,z,udp_vec[i]);
      dopValue[0] += tempVal[0];  // acceptor
      dopValue[1] += tempVal[1];  // donor
  }
  */
  for (std::size_t i=0; i < gdp_vec.size(); ++i)
  {
    tempVal=evalGaussianDoping(x,y,z,gdp_vec[i]);
    dopValue[0] += tempVal[0];  // acceptor
    dopValue[1] += tempVal[1];  // donor
  }
  for (std::size_t i=0; i < ldp_vec.size(); ++i)
  {
    tempVal=evalLinearDoping(x,y,z,ldp_vec[i]);
    dopValue[0] += tempVal[0];  // acceptor
    dopValue[1] += tempVal[1];  // donor
  }
  for (std::size_t i=0; i < edp_vec.size(); ++i)
  {
    tempVal=evalErfcDoping(x,y,z,edp_vec[i]);
    dopValue[0] += tempVal[0];  // acceptor
    dopValue[1] += tempVal[1];  // donor
  }
  for (std::size_t i=0; i < mgdp_vec.size(); ++i)
  {
    tempVal=evalMGaussDoping(x,y,z,mgdp_vec[i]);
    dopValue[0] += tempVal[0];  // acceptor
    dopValue[1] += tempVal[1];  // donor
  }

  for (std::size_t i=0; i < hdp_vec.size(); ++i)
  {
    tempVal=evalHaloDoping(x,y,z,hdp_vec[i]);
    dopValue[0] += tempVal[0];  // acceptor
    dopValue[1] += tempVal[1];  // donor
  }

  int fileCounter1D = 0;
  int fileCounter2D = 0; 
  int fileCounter3D = 0;
  int udCounter = 0;  

  for (ParameterList::ConstIterator model_it = dopParamList.begin();
       model_it != dopParamList.end(); ++model_it)
  {
    const string key = model_it->first;
    
    if (key.compare(0, 8, "Function") == 0)
    {
      const Teuchos::ParameterEntry& entry = model_it->second;
      const ParameterList& funcParamList = Teuchos::getValue<Teuchos::ParameterList>(entry);
      const string funcType = funcParamList.get<string>("Function Type");
      std::vector<double> tmpVal(2, 0.0);


      if (funcType == "MMS_RDH_1")
      {
      // This assumes the x coordinate is given in microns. The function
      // itself was originally formulated in CGS units thus the
      // multiplication by 1.0e-4.
        tmpVal = mms1->doping(x*1.0e-4);
      }

      if (funcType == "MMS_RDH_2")
      {
        tmpVal[0] = mms2->Na(x*1.0e-4);
        tmpVal[1] = mms2->Nd(x*1.0e-4);
      }
      
      if (funcType == "mms_nlp_glh_1")
      {
        // Convert microns to cm and invoke the function that calculates
        // the doping for the MMS NLP
        tmpVal[1] = mms_nlp_func->net_doping(x*1.0e-4);
      }

      if (funcType == "File2D")
      {
        tmpVal = evalFile2DDoping(fileCounter2D,x,y,z,funcParamList);
        fileCounter2D++;
      }

      if (funcType == "File3D")
      {
        tmpVal = evalFile3DDoping(fileCounter3D,x,y,z,funcParamList);
        fileCounter3D++;
      }
      
      if (funcType == "File1D")
      {
        tmpVal = evalFile1DDoping(fileCounter1D,x,y,z,funcParamList);
        fileCounter1D++;
      }  

      // compute uniform doping here due to addition of gauss decay
      if (funcType == "Uniform")
      { 
        tmpVal = evalUniformDoping(x,y,z,udp_vec[udCounter],udCounter,funcParamList);
        udCounter++; 
      }
 
      // add up doping
      dopValue[0] += tmpVal[0];  // acceptor
      dopValue[1] += tmpVal[1];  // donor

    }  // end of if (key.compare(0, 8, "Function") == 0)

  }  // end of for loop

  return dopValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalUniformDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalUniformDoping
  (const double& x, const double& y, const double& z, const uniformDopingParams& udp,
   int udCounter, const Teuchos::ParameterList& plist)
{
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  const string dopType = udp.dopType;
  double dopVal = udp.dopVal;
  double xmin = udp.xmin;
  double ymin = udp.ymin;
  double zmin = udp.zmin;
  double xmax = udp.xmax;
  double ymax = udp.ymax;
  double zmax = udp.zmax;
  double initialDopVal = udp.initialDopVal;
  double finalDopVal = udp.finalDopVal;

  if (udp.sweepMe)
    {
      double fraction = user_value->getRealValue();
      dopVal = initialDopVal*(1.0-fraction) + finalDopVal*fraction;
    }


  if ( (x >= xmin) && (x <= xmax) && (y >= ymin) && (y <= ymax) && (z >= zmin) && (z <= zmax) )
  {
    // add Gaussian decay to doping profile
    if (plist.isSublist("Gauss Decay"))
    {
      double decayFactor =  evalGaussDecayFactor(1, udCounter, x, y, z);
      dopVal *= decayFactor;
    }
     
    if (dopType == "Acceptor")
    {
      dopValue[0] = dopVal;  // acceptor doping
      dopValue[1] = 0;       // donor doping
    }
    else if (dopType == "Donor")
    {
      dopValue[0] = 0;
      dopValue[1] = dopVal;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
        << "Invalid Doping Type ! Must be Acceptor or Donor !");
  }

  // return 0 if (x,y,z) is outside the box region

  return dopValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalGaussianDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalGaussianDoping
  (const double& x, const double& y, const double& z, const gaussianDopingParams& gdp)
{
  return prof_eval->evalGaussianProfile(x,y,z,gdp);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleGaussian()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double DopingRaw_Function<EvalT, Traits>::evalSingleGaussian
  (const std::string axis, bool& found, const double& coord, const double& minDopVal,
   const double& maxDopVal, const double& min, const double& max, const double& loc,
   const double& width, const bool& checkAxis, const std::string& dir)
{
  return prof_eval->evalSingleGaussian(axis,found,coord,minDopVal,maxDopVal,
                                      min,max,loc,width,checkAxis,dir);
}



///////////////////////////////////////////////////////////////////////////////
//
//  evalLinearDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalLinearDoping
  (const double& x, const double& y, const double& z, const linearDopingParams& ldp)
{
  return prof_eval->evalLinearProfile(x,y,z,ldp);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleLinear()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double DopingRaw_Function<EvalT, Traits>::evalSingleLinear
  (const std::string axis, bool& found, const double& coord, const double& min,
   const double& max, const bool& checkAxis)
{
  return prof_eval->evalSingleLinear(axis,found,coord,min,max,checkAxis);
}



///////////////////////////////////////////////////////////////////////////////
//
//  evalErfcDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalErfcDoping
  (const double& x, const double& y, const double& z, const erfcDopingParams& edp)
{
  return prof_eval->evalErfcProfile(x,y,z,edp);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleErfc()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double DopingRaw_Function<EvalT, Traits>::evalSingleErfc
  (const std::string axis, bool& found, const double& coord, const double& minDopVal,
   const double& maxDopVal, const double& min, const double& max, const double& loc,
   const double& width, const bool& checkAxis, const std::string& dir)
{
  return prof_eval->evalSingleErfc(axis,found,coord,minDopVal,maxDopVal,
                                   min,max,loc,width,checkAxis,dir);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalMGaussDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalMGaussDoping
  (const double& x, const double& y, const double& z, const mgaussDopingParams& mgdp)
{
  return prof_eval->evalMGaussProfile(x,y,z,mgdp);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalSingleMGauss()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double DopingRaw_Function<EvalT, Traits>::evalSingleMGauss(
   const std::string axis, bool& found, const double& coord,
   const double& minDopVal, const double& maxDopVal,
   const double& min, const double& max, const bool& checkErfc,
   const double& width, const bool& checkAxis)
{
  return prof_eval->evalSingleMGauss(axis,found,coord,minDopVal,maxDopVal,
				     min,max,checkErfc,width,checkAxis);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evalHaloDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::
evalHaloDoping(const double& x, const double& y, const double& z, const haloDopingParams& hdp)
{
  return prof_eval->evalHaloProfile(x,y,z,hdp);
}


//**********************************************************************
// Exports 2D doping profile from an external file and maps doping values
// to each (x,y,z) point if (x,y) is within the doping profile range
// Notes:
// 1) In 3D case, doping along z-axis is assumed to be homogeneous
// 2) Doping type is determined by Doping Type
// 3) Three-column file format: x y dop_value(x,y) (dop_value must be positive)
// 4) A nearest-neighbor search is performed for every (x,y,z)
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  evalFile2DDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalFile2DDoping
(int fileCounter, const double& x, const double& y, const double& z, const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::vector;

  doping_struct r;
  vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  double inverse_power = 0.0; // default is the nearest neighbor only
  double buffer = 0.0;

  if (plist.isParameter("Inverse Power")) 
    inverse_power = plist.get<double>("Inverse Power");
  if (plist.isParameter("Buffer")) 
    buffer = plist.get<double>("Buffer");
  
  string dopingType = plist.get<std::string>("Doping Type");

  if ((x+buffer >= MinX[fileCounter]) && (x-buffer <= MaxX[fileCounter]) &&
      (y+buffer >= MinY[fileCounter]) && (y-buffer <= MaxY[fileCounter]) )
  {
    r.x = x; r.y = y;

    double closest_distance=1e100, dop_value = 0.0, weight = 1.0;

    for (size_t i=0; i<=DV[fileCounter].size()-1; i++)
    {
     double distance = DV[fileCounter][i].distance(r);

     if (distance == 0.0)
     {
      weight=1.0;
      dop_value = DV[fileCounter][i].d;
      break;
     }

     if (inverse_power>0.0) // inverse weighted average
     {
       distance = pow(distance, inverse_power);
       weight+=1.0/distance;
       dop_value+= DV[fileCounter][i].d/distance;
     }
     else if (distance < closest_distance) //nearest neighbor
     {
      closest_distance = distance;
      dop_value = DV[fileCounter][i].d;
     }
    }
    dop_value = dop_value/weight;

    // assign value according to Doping Type
    if (dopingType == "Acceptor")
    {
      dopValue[0] = dop_value;
      dopValue[1] = 0.;
    }
    else if (dopingType == "Donor")
    { 
      dopValue[0] = 0.;
      dopValue[1] = dop_value;
    }
    else
    {
      //error out
    }
  }

  return dopValue;
}


//**********************************************************************
// Imports 3D doping profile from an external file and maps doping values
// to each (x,y,z) point if (x,y,z) is within the doping profile range
// Notes:
// 1) Doping type is determined by Doping Type
// 2) Four-column file format: x y z dop_value(x,y) (dop_value must be positive)
// 3) A nearest-neighbor search is performed for every (x,y,z)
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  evalFile3DDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalFile3DDoping(
      int fileCounter, const double& x, const double& y, 
      const double& z, const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::vector;

  doping_struct r;
  vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  double inverse_power = 0.0; // default is the nearest neighbor only
  double buffer=0.0;

  if (plist.isParameter("Inverse Power")) 
    inverse_power = plist.get<double>("Inverse Power");
  if (plist.isParameter("Buffer")) 
    buffer = plist.get<double>("Buffer");
    
  string dopingType = plist.get<std::string>("Doping Type");

  if ((x+buffer >= MinX[fileCounter]) && (x-buffer <= MaxX[fileCounter]) &&
      (y+buffer >= MinY[fileCounter]) && (y-buffer <= MaxY[fileCounter]) &&
      (z+buffer >= MinZ[fileCounter]) && (z-buffer <= MaxZ[fileCounter]) )
  {
    r.x = x;
    r.y = y;
    r.z = z;

    double closest_distance=1e100, dop_value = 0.0, weight = 1.0;

    for (size_t i=0; i<=DV[fileCounter].size()-1; i++)
    {
     double distance = DV[fileCounter][i].distance(r);

     if (distance == 0.0)
     {
      weight=1.0;
      dop_value = DV[fileCounter][i].d;
      break;
     }

     if (inverse_power>0.0) // inverse weighted average
     {
       distance = pow(distance, inverse_power);
       weight+=1.0/distance;
       dop_value+= DV[fileCounter][i].d/distance;
     }
     else if (distance < closest_distance) //nearest neighbor
     {
      closest_distance = distance;
      dop_value = DV[fileCounter][i].d;
     }
    }
    dop_value = dop_value/weight;

    // assign value according to Doping Type
    if (dopingType == "Acceptor")
    {
      dopValue[0] = dop_value;
      dopValue[1] = 0.;
    }
    else if(dopingType == "Donor")
    { 
      dopValue[0] = 0.;
      dopValue[1] = dop_value;
    }
    else
    {
      //error out
    }
  }

  return dopValue;
}


//**********************************************************************
// Exports 1D doping profile from an external file and maps doping values
// to each (x,y,z) point if x is within the user-defined range
// Notes:
// 1) The direction of doping is specified by Doping Direction and uniform in other 2 directions
// 2) Doping type is determined by Doping Type
// 3) Two-column file format: x dop_value(x) (dop_value must be positive)
// 4) A nearest-neighbor search is performed for every (x,y,z)
//**********************************************************************

///////////////////////////////////////////////////////////////////////////////
//
//  evalFile1DDoping()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
std::vector<double> DopingRaw_Function<EvalT, Traits>::evalFile1DDoping
  (int fileCounter, const double& x, const double& y, const double& z, const Teuchos::ParameterList& plist)
{
  using Teuchos::ParameterList;
  using std::string;
  using std::cout;
  using std::endl;
  using std::ifstream;
  using std::vector;

  doping_struct_1D r;

  vector<double> dopValue(2, 0.0);
  TEUCHOS_ASSERT(!(dopValue.size() > 2));

  string dopingType = plist.get<string>("Doping Type");
  TEUCHOS_ASSERT(dopingType.size() > 0);
  
  string direction = plist.get<string>("Doping Direction");
  TEUCHOS_ASSERT(direction.size() > 0);

  double buffer = 0.0;
  if (plist.isParameter("Buffer")) 
    buffer = plist.get<double>("Buffer");
  
  double xmin = -1e100, ymin = -1e100, zmin = -1e100;
  double xmax =  1e100, ymax =  1e100, zmax =  1e100; 

  if (plist.isParameter("Xmin"))  xmin = plist.get<double>("Xmin");
  if (plist.isParameter("Xmax"))  xmax = plist.get<double>("Xmax");
  if (plist.isParameter("Ymin"))  ymin = plist.get<double>("Ymin");
  if (plist.isParameter("Ymax"))  ymax = plist.get<double>("Ymax");
  if (plist.isParameter("Zmin"))  zmin = plist.get<double>("Zmin");
  if (plist.isParameter("Zmax"))  zmax = plist.get<double>("Zmax");

  size_t end = 0; 
  
  // determine the doping direction and modify min/max values
  if (direction == "X")
  {
    double minFile, maxFile; 
    r.x = x; 
    minFile = DV1[fileCounter][0].x - buffer;  // ascendingly sorted in readDopingFile1D
    end = DV1[fileCounter].size()-1;
    maxFile = DV1[fileCounter][end].x + buffer;
    if (minFile > xmin) xmin = minFile;   // use a smaller region between user-specified and file-contained
    if (maxFile < xmax) xmax = maxFile;
  }  
  else if (direction == "Y")
  {
    double minFile, maxFile; 
    r.x = y;
    minFile = DV1[fileCounter][0].x - buffer;  // ascendingly sorted in readDopingFile1D
    end = DV1[fileCounter].size()-1;
    maxFile = DV1[fileCounter][end].x + buffer;
    if (minFile > ymin) ymin = minFile;   // use a smaller region between user-specified and file-contained
    if (maxFile < ymax) ymax = maxFile;
  }
  else if (direction == "Z")
  {
    double minFile, maxFile; 
    r.x = z; 
    minFile = DV1[fileCounter][0].x - buffer;  // ascendingly sorted in readDopingFile1D
    end = DV1[fileCounter].size()-1;
    maxFile = DV1[fileCounter][end].x + buffer;
    if (minFile > zmin) zmin = minFile;   // use a smaller region between user-specified and file-contained
    if (maxFile < zmax) zmax = maxFile;
  }
  else
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! Doping Direction must be either X, Y, or Z!" << std::endl);

  if (x >= xmin && x <= xmax && 
      y >= ymin && y <= ymax &&
      z >= zmin && z <= zmax   )
  {
    size_t closest_point = 0;
    double closest_distance = DV1[fileCounter][0].distance(r);
    
    // find the closest point 
    for (size_t i = 1; i <= DV1[fileCounter].size()-1; i++)
    {
     double distance = DV1[fileCounter][i].distance(r);
     if (distance < closest_distance)
     {
      closest_distance = distance;
      closest_point = i;
     }
    }
    
    double signedDist = r.x - DV1[fileCounter][closest_point].x; 
    double dop_value = 0.0; 
    
    if (signedDist > 0.0) // r.x is located within closest_point and closest_point+1
    {
      double dop1 = DV1[fileCounter][closest_point].d;
      double x1 = DV1[fileCounter][closest_point].x;
      
      if (closest_point == DV1[fileCounter].size()-1)  // last point (also largest x point)
        dop_value = dop1;  // use the last doping value
        
      else  // perform linear interpolation
      {
        double dop2 = DV1[fileCounter][closest_point+1].d;
        double x2 = DV1[fileCounter][closest_point+1].x;
        double slope = (dop1-dop2)/(x1-x2);
        double intercept = (dop1+dop2 - slope*(x1+x2)) * 0.5; 
        dop_value = slope*r.x + intercept; 
      }
    }
    
    else if (signedDist < 0.0)  // r.x is located within closest_point and closest_point-1
    {
      double dop1 = DV1[fileCounter][closest_point].d;
      double x1 = DV1[fileCounter][closest_point].x;
      
      if (closest_point == 0)  // first point (also smallest x point)
        dop_value = dop1;  // use the first doping value
        
      else  // perform linear interpolation
      {
        double dop2 = DV1[fileCounter][closest_point-1].d;
        double x2 = DV1[fileCounter][closest_point-1].x;
        double slope = (dop1-dop2)/(x1-x2);
        double intercept = (dop1+dop2 - slope*(x1+x2)) * 0.5; 
        dop_value = slope*r.x + intercept; 
      }
    }
    else // signedDist == 0
      dop_value = DV1[fileCounter][closest_point].d;

    // add Gaussian decay to doping profile
    if (plist.isSublist("Gauss Decay"))
    {
      double decayFactor =  evalGaussDecayFactor(0, fileCounter, x, y, z);
      dop_value *= decayFactor; 
    }
    
    // assign value according to Doping Type
    if (dopingType == "Acceptor")
    {
      dopValue[0] = dop_value;
      dopValue[1] = 0.;
    }
    else if(dopingType == "Donor")
    { 
      dopValue[0] = 0.;
      dopValue[1] = dop_value;
    }
    else
    {
      //error out
    }
    
  }  // within a specified box

  return dopValue;
}


///////////////////////////////////////////////////////////////////////////////
//
//  readGaussDecayParams(), type=0 for File1D, 1 for Uniform
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
void DopingRaw_Function<EvalT, Traits>::readGaussDecayParams
  (int type, int counter, const Teuchos::ParameterList& plist)   
{
  using Teuchos::ParameterList;
  using std::string;
  using std::vector;

  // support 6 gauss decays at most (+/- x, +/- y, and/or +/- z)
  for (ParameterList::ConstIterator it = plist.begin(); it != plist.end(); ++it)
  {
    const string key = it->first;
    if (key == "Decay Direction 1")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 1"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 1"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 1"));
    }      
    else if (key == "Decay Direction 2")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 2"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 2"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 2"));
    }
    else if (key == "Decay Direction 3")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 3"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 3"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 3"));
    }
    else if (key == "Decay Direction 4")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 4"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 4"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 4"));
    }
    else if (key == "Decay Direction 5")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 5"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 5"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 5"));
    }
    else if (key == "Decay Direction 6")
    {
      decayDir[type][counter].push_back(plist.get<string>("Decay Direction 6"));
      decayPos[type][counter].push_back(plist.get<double>("Decay Position 6"));
      decayWidth[type][counter].push_back(plist.get<double>("Decay Width 6"));
    }
  }
}  


///////////////////////////////////////////////////////////////////////////////
//
//  evalGaussDecayFactor()
//
///////////////////////////////////////////////////////////////////////////////
template<typename EvalT, typename Traits>
double DopingRaw_Function<EvalT, Traits>::evalGaussDecayFactor(int type,
  int counter, const double& x, const double& y, const double& z)
{
  using std::string;

  // loop over the number of Gauss decays for a given doping file
  double decayFactor = 1.0; 
  for (unsigned int id = 0; id < decayDir[type][counter].size(); ++id)
  {
    string dir = decayDir[type][counter][id]; 
    double pos = decayPos[type][counter][id];
    double width = decayWidth[type][counter][id];

    if ((dir == "X Positive") && (x > pos))
      decayFactor *= std::exp(-(x-pos)*(x-pos)/(2.0 *width*width));
    
    else if ((dir == "X Negative") && (x < pos))
      decayFactor *= std::exp(-(x-pos)*(x-pos)/(2.0*width*width));
    
    else if ((dir == "Y Positive") && (y > pos))
      decayFactor *= std::exp(-(y-pos)*(y-pos)/(2.0*width*width));

    else if ((dir == "Y Negative") && (y < pos))
      decayFactor *= std::exp(-(y-pos)*(y-pos)/(2.0*width*width));
  
    else if ((dir == "Z Positive") && (z > pos))
      decayFactor *= std::exp(-(z-pos)*(z-pos)/(2.0*width*width));

    else if ((dir == "Z Negative") && (z < pos))
      decayFactor *= std::exp(-(z-pos)*(z-pos)/(2.0*width*width));
  }
 
  return decayFactor;  
}




////////////////////////////////////////////////////////////////////

// Generic Profile Calculator
ProfileEvals::ProfileEvals(int dim) {
  num_dim = dim;
}


// evaluate a single-axis linear (x, or y, or z)
double ProfileEvals::evalSingleLinear(
                          const std::string axis, bool& found, 
			  const double& coord, const double& min, 
                          const double& max, const bool& checkAxis) {
  // If Linear is NOT specified for a certain axis (say X), returns 1.0
  double LinearVal = 1.0;

  // Linear is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if Linear is set along an axis, then found = true

    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
     LinearVal = (coord-min)/(max-min);
    else
     LinearVal = -1.0; // -1 is a flag indicating that the point is outside the [min,max] range
  }

  return LinearVal;
}


// evaluate linear profile at given (x,y,z)
std::vector<double> ProfileEvals::evalLinearProfile(
	           const double& x, const double& y,
		   const double& z, const linearDopingParams& ldp) {
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> profValue(2, 0.0);
  TEUCHOS_ASSERT(!(profValue.size() > 2));

  const string profType = ldp.dopType;
  const double startVal = ldp.startVal;
  const double endVal = ldp.endVal;

  // x direction
  const double x_min = ldp.x_min;
  const double x_max = ldp.x_max;
  const bool x_checkAxis = ldp.x_checkAxis;

  // y direction
  const double y_min = ldp.y_min;
  const double y_max = ldp.y_max;
  const bool y_checkAxis = ldp.y_checkAxis;

  // z direction
  const double z_min = ldp.z_min;
  const double z_max = ldp.z_max;
  const bool z_checkAxis = ldp.z_checkAxis;

  bool found = false;
  double xLinearVal = 1.0, yLinearVal = 1.0, zLinearVal = 1.0;

  xLinearVal = evalSingleLinear("X", found, x, x_min, x_max, x_checkAxis);
  if (num_dim == 2)
    yLinearVal = evalSingleLinear("Y", found, y, y_min, y_max, y_checkAxis);
  if (num_dim == 3)
  {
    yLinearVal = evalSingleLinear("Y", found, y, y_min, y_max, y_checkAxis);
    zLinearVal = evalSingleLinear("Z", found, z, z_min, z_max, z_checkAxis);
  }

  // throw exception if NO Linear profile is specified
  if (!found)
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Linear function is specified "
     << "for doping Function Type of Linear! At least one Linear function along "
     << "x, y, or z must be specified! ");

  if (xLinearVal>=0.0 && yLinearVal>=0.0 && zLinearVal>=0.0)
  {
    // assign value according to Doping Type
    if (profType == "Acceptor")
    {
      profValue[0] = (endVal-startVal)*xLinearVal*yLinearVal*zLinearVal+startVal;
      profValue[1] = 0.;
    }
    else if (profType == "Donor")
    {
      profValue[0] = 0.;
      profValue[1] = (endVal-startVal)*xLinearVal*yLinearVal*zLinearVal+startVal;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Invalid Profile Type ! Must be Acceptor or Donor !");
  }
  else
  {
    // The point is outside the Linear definition [min..max] range: doping is not modified
    profValue[0] = 0.;
    profValue[1] = 0.;
  }

  return profValue;

}


// evaluate a single-axis gaussian (x, or y, or z)
double ProfileEvals::evalSingleGaussian(
			  const std::string axis, bool& found, 
			  const double& coord, const double& minProfVal,
			  const double& maxProfVal, const double& min, 
			  const double& max, const double& loc,
			  const double& width, const bool& checkAxis, 
			  const std::string& dir) {
  using std::string;
  using Teuchos::ParameterList;

  // If Gaussian is NOT specified for a certain axis (say X), returns 1.0
  double gaussVal = 1.0;

  // set gaussVal to 0 when coord is outside [min max] range
  if ( (coord < min) || (coord > max) )  gaussVal = 0.0;

  //Gaussian is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if a Gaussian is set along an axis, then found = true

    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
    {
      if (dir == "Both")
        gaussVal = exp(-log(maxProfVal/minProfVal) * pow((coord-loc)/width, 2.0) );
      else if (dir == "Positive")
      {
        if (coord >= loc)
          gaussVal = exp(-log(maxProfVal/minProfVal) * pow((coord-loc)/width, 2.0));
        else
          gaussVal = 1.0;
      }
      else if (dir == "Negative")
      {
        if (coord <= loc)
          gaussVal = exp(-log(maxProfVal/minProfVal) * pow((coord-loc)/width, 2.0));
        else
          gaussVal = 1.0;
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Direction must be either Both, Positive, or Negative !");
    }

    // If Gaussian is specified for a certain axis but coord is outside [min, max], returns 0.0
    else
      gaussVal = 0.0; // 0 outside the [min, max) range

  }  // end of if (plist.isParameter(axis+"Peak Location") ...)

  return gaussVal;
}


// evaluate gaussian profile at given (x,y,z)
std::vector<double> ProfileEvals::evalGaussianProfile(
	const double& x, const double& y, const double& z, 
	const gaussianDopingParams& gdp) {
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> profValue(2, 0.0);
  TEUCHOS_ASSERT(!(profValue.size() > 2));

  const string profType = gdp.dopType;
  const double maxVal = gdp.maxVal;
  const double minVal = gdp.minVal;

  // x direction
  const string x_dir = gdp.x_dir;
  const double x_loc = gdp.x_loc;
  const double x_width = gdp.x_width;
  const double x_min = gdp.x_min;
  const double x_max = gdp.x_max;
  const bool x_checkAxis = gdp.x_checkAxis;

  // y direction
  const string y_dir = gdp.y_dir;
  const double y_loc = gdp.y_loc;
  const double y_width = gdp.y_width;
  const double y_min = gdp.y_min;
  const double y_max = gdp.y_max;
  const bool y_checkAxis = gdp.y_checkAxis;

  // z direction
  const string z_dir = gdp.z_dir;
  const double z_loc = gdp.z_loc;
  const double z_width = gdp.z_width;
  const double z_min = gdp.z_min;
  const double z_max = gdp.z_max;
  const bool z_checkAxis = gdp.z_checkAxis;

  bool found = false;
  double xGaussVal = 1.0, yGaussVal = 1.0, zGaussVal = 1.0;

  xGaussVal = evalSingleGaussian("X", found, x, minVal, maxVal, x_min, x_max, x_loc, x_width, x_checkAxis, x_dir);
  if (num_dim == 2)
    yGaussVal = evalSingleGaussian("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
  if (num_dim == 3)
  {
    yGaussVal = evalSingleGaussian("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
    zGaussVal = evalSingleGaussian("Z", found, z, minVal, maxVal, z_min, z_max, z_loc, z_width, z_checkAxis, z_dir);
  }

  // throw exception if NO Gaussian profile is specified
  if (!found)
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Gaussian is specified "
     << "for profile Function Type of Gauss/Gaussian! At least one Gaussian along "
     << "x, y, or z must be specified! ");

  // assign value according to Doping Type
  if (profType == "Acceptor")
  {
    profValue[0] = maxVal*xGaussVal*yGaussVal*zGaussVal;
    profValue[1] = 0.;
  }
  else if (profType == "Donor")
  {
    profValue[0] = 0.;
    profValue[1] = maxVal*xGaussVal*yGaussVal*zGaussVal;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Invalid Profile Type ! Must be Acceptor or Donor !");
  return profValue;
}


// evaluate a single-axis erfc (x, or y, or z)
double ProfileEvals::evalSingleErfc(
	 const std::string axis, bool& found, 
	 const double& coord, const double& /* minProfVal */,
	 const double& /* maxProfVal */, const double& min, 
	 const double& max, const double& loc,
	 const double& width, const bool& checkAxis, 
	 const std::string& dir) {
  using std::string;
  using Teuchos::ParameterList;

  // If Erfc is NOT specified for a certain axis (say X), returns 1.0
  double ErfcVal = 1.0;

  // Erfc is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if a Erfc is set along an axis, then found = true

    // within [min, max] range
    if ( (coord >= min) && (coord <= max) )
    {
      if (dir == "Positive")
        ErfcVal = 0.5*erfc( (coord-loc)/width );
      else if (dir == "Negative")
        ErfcVal = 0.5*erfc(-(coord-loc)/width );
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
          << "Error ! " << axis << " Direction must be Positive or Negative !");
    }
    else
      ErfcVal = -1.0; // -1 is a flag indicating that the point is outside the [min,max] range

  }  // end of if (plist.isParameter(axis+"Bend Location") ...)

  return ErfcVal;
}


// evaluate erfc profile at given (x,y,z)
std::vector<double> ProfileEvals::evalErfcProfile(
        const double& x, const double& y, const double& z, 
	const erfcDopingParams& edp) {
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> profValue(2, 0.0);
  TEUCHOS_ASSERT(!(profValue.size() > 2));

  const string profType = edp.dopType;
  const double maxVal = edp.maxVal;
  const double minVal = edp.minVal;

  // x direction
  const string x_dir = edp.x_dir;
  const double x_loc = edp.x_loc;
  const double x_width = edp.x_width;
  const double x_min = edp.x_min;
  const double x_max = edp.x_max;
  const bool x_checkAxis = edp.x_checkAxis;

  // y direction
  const string y_dir = edp.y_dir;
  const double y_loc = edp.y_loc;
  const double y_width = edp.y_width;
  const double y_min = edp.y_min;
  const double y_max = edp.y_max;
  const bool y_checkAxis = edp.y_checkAxis;

  // z direction
  const string z_dir = edp.z_dir;
  const double z_loc = edp.z_loc;
  const double z_width = edp.z_width;
  const double z_min = edp.z_min;
  const double z_max = edp.z_max;
  const bool z_checkAxis = edp.z_checkAxis;

  if ((maxVal < 0.) || (minVal < 0.))
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Error ! Erfc profile Max and Min Values must be greater than 0.");

  bool found = false;
  double xErfcVal = 1.0, yErfcVal = 1.0, zErfcVal = 1.0;

  xErfcVal = evalSingleErfc("X", found, x, minVal, maxVal, x_min, x_max, x_loc, x_width, x_checkAxis, x_dir);
  if (num_dim == 2)
    yErfcVal = evalSingleErfc("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
  if (num_dim == 3)
  {
    yErfcVal = evalSingleErfc("Y", found, y, minVal, maxVal, y_min, y_max, y_loc, y_width, y_checkAxis, y_dir);
    zErfcVal = evalSingleErfc("Z", found, z, minVal, maxVal, z_min, z_max, z_loc, z_width, z_checkAxis, z_dir);
  }

  // throw exception if NO Erfc profile is specified
  if (!found)
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Erfc is specified "
     << "for profile Function Type of Erfc! At least one Erfc profile along "
     << "x, y, or z must be specified! ");

  if (xErfcVal>=0.0 && yErfcVal>=0.0 && zErfcVal>=0.0)
  {
    // assign value according to Doping Type
    if (profType == "Acceptor")
    {
      profValue[0] = minVal*pow(maxVal/minVal,xErfcVal*yErfcVal*zErfcVal);
      profValue[1] = 0.;
    }
    else if (profType == "Donor")
    {
      profValue[0] = 0.;
      profValue[1] = minVal*pow(maxVal/minVal,xErfcVal*yErfcVal*zErfcVal);;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
           << "Invalid Doping Type ! Must be Acceptor or Donor !");
  }
  else
  {
    // The point is outside the Erfc definition [min..max] range: doping is not modified
    profValue[0] = 0.;
    profValue[1] = 0.;
  }

  return profValue;
}


// evaluate a single-axis MGauss (x, or y, or z)
double ProfileEvals::evalSingleMGauss(
		 const std::string axis, bool& found, 
		 const double& coord, const double& minProfVal, 
		 const double& maxProfVal, const double& min, 
		 const double& max, const bool& checkErfc,
		 const double& width, const bool& checkAxis) {
  using std::string;
  using Teuchos::ParameterList;

  // If MGauss is NOT specified for a certain axis (say X), returns 1.0
  double MGaussVal = 1.0;

  // MGauss is specified along a given axis
  if (checkAxis)
  {
    found = true;  // if a MGauss is set along an axis, then found = true

    if (checkErfc)
    {
      // Erfc along this axis
      MGaussVal = 0.5*(erfc((coord-max)/width)-erfc((coord-min)/width));
    }
    else
    {
      // Gaussian profile along this axis
      if (coord < min)
      {
        if (minProfVal > 0.0)  // Profile Min Value is given
          MGaussVal = std::exp(-std::log(maxProfVal/minProfVal) * std::pow((coord-min)/width, 2.0));
        else
          MGaussVal = std::exp(-(coord-min)*(coord-min)/width/width);
      }
      else if (coord > max)
      {
        if (minProfVal > 0.0)  // Profile Min Value is given
          MGaussVal = std::exp(-std::log(maxProfVal/minProfVal) * std::pow((coord-max)/width, 2.0));
        else
          MGaussVal = std::exp(-(coord-max)*(coord-max)/width/width);
      }
      else
        MGaussVal = 1.0;
    }

  }  // end of if (plist.isParameter(axis+" Width"))

  return MGaussVal;

}


// evaluate gaussian (=MGauss) profile at given (x,y,z)
std::vector<double> ProfileEvals::evalMGaussProfile(
	const double& x, const double& y, const double& z, 
        const mgaussDopingParams& mgdp) {
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> profValue(2, 0.0);
  TEUCHOS_ASSERT(!(profValue.size() > 2));

  const string profType = mgdp.dopType;
  const double maxVal = mgdp.maxVal;
  const double minVal = mgdp.minVal;

  // x direction
  const double x_width = mgdp.x_width;
  const double x_min = mgdp.x_min;
  const double x_max = mgdp.x_max;
  const bool x_checkErfc = mgdp.x_checkErfc;
  const bool x_checkAxis = mgdp.x_checkAxis;

  // y direction
  const double y_width = mgdp.y_width;
  const double y_min = mgdp.y_min;
  const double y_max = mgdp.y_max;
  const bool y_checkErfc = mgdp.y_checkErfc;
  const bool y_checkAxis = mgdp.y_checkAxis;

  // z direction
  const double z_width = mgdp.z_width;
  const double z_min = mgdp.z_min;
  const double z_max = mgdp.z_max;
  const bool z_checkErfc = mgdp.z_checkErfc;
  const bool z_checkAxis = mgdp.z_checkAxis;

  bool found = false;
  double xMGaussVal = 1.0, yMGaussVal = 1.0, zMGaussVal = 1.0;

  xMGaussVal = evalSingleMGauss("X", found, x, minVal, maxVal, x_min, x_max, x_checkErfc, x_width, x_checkAxis);
  if (num_dim == 2)
    yMGaussVal = evalSingleMGauss("Y", found, y, minVal, maxVal, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
  if (num_dim == 3)
  {
    yMGaussVal = evalSingleMGauss("Y", found, y, minVal, maxVal, y_min, y_max, y_checkErfc, y_width, y_checkAxis);
    zMGaussVal = evalSingleMGauss("Z", found, z, minVal, maxVal, z_min, z_max, z_checkErfc, z_width, z_checkAxis);
  }

  // throw exception if NO MGauss profile is specified
  if (!found)
   TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Error! No Gaussian is specified "
     << "for doping Function Type of MGauss! At least one MGauss profile along "
     << "x, y, or z must be specified! ");

  // assign value according to Doping Type
  if (profType == "Acceptor")
  {
    profValue[0] = maxVal*xMGaussVal*yMGaussVal*zMGaussVal;
    profValue[1] = 0.;
  }
  else if (profType == "Donor")
  {
    profValue[0] = 0.;
    profValue[1] = maxVal*xMGaussVal*yMGaussVal*zMGaussVal;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, std::endl
       << "Invalid Doping Type ! Must be Acceptor or Donor !");

  return profValue;
}


// evaluate halo profile
std::vector<double> ProfileEvals::evalHaloProfile(
        const double& x, const double& y, const double& z, 
        const haloDopingParams& hdp) {
  using std::string;
  using Teuchos::ParameterList;

  std::vector<double> profValue(2, 0.0);
  TEUCHOS_ASSERT(!(profValue.size() > 2));

  const string profType = hdp.dopType;
  const double profileVal = hdp.dopingVal;
  const double minProfileVal = hdp.minDopingVal;
  const std::string distributionType = hdp.distributionType;

  // x direction
  const double x_center = hdp.x_center;
  //const bool x_checkAxis = hdp.x_checkAxis;

  // y direction
  const double y_center = hdp.y_center;
  //const bool y_checkAxis = hdp.y_checkAxis;

  // z direction
  const double z_center = hdp.z_center;
  //const bool z_checkAxis = hdp.z_checkAxis;

  //rotation and critical radii
  const double r1 = hdp.r1;
  const double r2 = hdp.r2;
  const double r3 = 1e10;
  const double rotation = hdp.rotation;

  //Gaussian parameters
  const double width = hdp.width;

  //It is presently assumend that rotation will be around the z axis.
  //Halo is usually created by ion implantation, so the doping is constant within and zero without

  //Translate for center
  double x_loc_c = x-x_center;
  double y_loc_c = y-y_center;
  double z_loc_c = z-z_center;
  
  //Rotate
  //Convert rotation angle from degrees to radians
  double pi = acos(-1.0);
  double angle = pi*rotation/180.0;

  //Map the coordinates into the reference through rotation.
  double x_loc =  x_loc_c*cos(-angle) + y_loc_c*sin(-angle);
  double y_loc = -x_loc_c*sin(-angle) + y_loc_c*cos(-angle);
  double z_loc = z_loc_c;

  double functionEvalX = (x_loc)*(x_loc)/(r1*r1);
  double functionEvalY = (y_loc)*(y_loc)/(r2*r2);
  double functionEvalZ = (z_loc)*(z_loc)/(r3*r3);

  double functionEval = functionEvalX + functionEvalY + functionEvalZ;

  double profMag = 0.0;
  if (functionEval <= 1.0)
    {
      profMag = profileVal;
      std::vector<std::vector<double> > a;
      std::vector<double> b;
      double minDist = 0.0;

      std::string xdir,ydir,zdir;
      bool found=true;

      if (distributionType == "Gaussian")
	{
	  std::vector<std::vector<double> > a;
	  std::vector<double> b;
	  std::vector<double> x;
	  x.resize(3);
	  if (x_loc >= 0)
	    x[0] = r1;
	  else
	    x[0] = -r1;
	  if (y_loc >=0)
	    x[1] = r2;
	  else
	    x[1] = -r2;
	  x[2] = 0.0;
	  b.resize(3,0.);
	  a.resize(3);
	  for (size_t i=0 ; i<a.size() ; ++i)
	    a[i].resize(a.size());
	  
	  //Iterate
	  int itnmax = 20;
	  int itn = 0;
	  double crit = 1e-8;
	  double res = 0.0;
	  
	  while (itn < itnmax)
	    {
	      //RHS
	      b[0] = -( 2.0*(x[0] - x_loc) + x[2] * 2.0 * x[0] /(r1*r1));
	      b[1] = -( 2.0*(x[1] - y_loc) + x[2] * 2.0 * x[1] /(r2*r2));
	      b[2] = -(  (x[0])*(x[0]) /(r1*r1) + (x[1])*(x[1]) /(r2*r2) - 1.0);
	      
	      res = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
	      minDist = sqrt((x[0] - x_loc) * (x[0] - x_loc)  + (x[1] - y_loc) * (x[1] - y_loc) );
	      if (res < crit) break;
	      
	      //Matrix 
	      a[0][0] = 2.0 + 2.0*x[2]/r1/r1;
	      a[0][1] = 0.0;
	      a[0][2] = 2.0 * x[0] / r1/r1;
	      a[1][0] = 0.0;
	      a[1][1] = 2.0 + 2.0*x[2]/r2/r2;
	      a[1][2] = 2.0 * x[1] / r2/r2;
	      a[2][0] = 2.0 * x[0] /r1/r1;
	      a[2][1] = 2.0 * x[1] /r2/r2;
	      a[2][2] = 0.0;
	      
	      int success = lusolve(a,3,b);
	      if (success == 1)
		for (size_t i=0 ; i<b.size() ; ++i)
		  x[i] += b[i];
	      ++itn;
	    }
	  
	  x.clear();
	  b.clear();
	  for(size_t i=0 ; i<a.size() ; ++i)
	    a[i].clear();
	  a.clear();
	  
	  std::string dir = "Positive";
	  double coord = 0.0;
	  if (minDist < width)
	    {
	      coord = 1.0 - minDist/width;
	      if (coord < 0) coord = 0.0;
	    }
	  double gauss = evalSingleGaussian("X",found,coord,minProfileVal,profileVal,0,1,0.0,1.0,found,dir);
	  profMag = profileVal*gauss + (1-gauss)*minProfileVal;

	}

      else if(distributionType == "Uniform")
	profMag = profileVal;
    }

  // assign value according to Profile Type
  if (profType == "Acceptor")
  {
    profValue[0] = profMag;
    profValue[1] = 0.;
  }
  else if (profType == "Donor")
  {
    profValue[0] = 0.;
    profValue[1] = profMag;
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, std::endl
      << "Invalid Profile Type ! Must be Acceptor or Donor !");

  return profValue;
}




} // namespace charon

#endif
