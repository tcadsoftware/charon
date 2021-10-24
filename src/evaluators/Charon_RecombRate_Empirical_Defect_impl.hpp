
#ifndef CHARON_RECOMBRATE_EMPIRICAL_DEFECT_IMPL_HPP
#define CHARON_RECOMBRATE_EMPIRICAL_DEFECT_IMPL_HPP

#include <cmath>
#include <iostream>
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"

#include "Panzer_GatherIntegrationCoordinates.hpp"
#include "Panzer_GatherBasisCoordinates.hpp"

#include "Charon_Names.hpp"
#include "Charon_empiricalConvolution.hpp"

namespace charon {

///////////////////////////////////////////////////////////////////////////////
//
//  getCoordField() (Point specialization)
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits, typename PointType>
PHX::MDField<const typename RecombRate_Empirical_Defect<EvalT,Traits,PointType>::ScalarT,panzer::Cell,PointType,panzer::Dim>
RecombRate_Empirical_Defect<EvalT,Traits,PointType>::
getCoordField(Type<panzer::Point>, panzer::PureBasis const& /* basis */, panzer::IntegrationRule const& ir)
{
  const std::string coordName = panzer::GatherIntegrationCoordinates<EvalT,Traits>::fieldName(ir.cubature_degree);
  return PHX::MDField<const ScalarT,panzer::Cell,PointType,panzer::Dim>(coordName, ir.dl_vector);
}

///////////////////////////////////////////////////////////////////////////////
//
//  getCoordField() (BASIS specialization)
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits, typename PointType>
PHX::MDField<const typename RecombRate_Empirical_Defect<EvalT,Traits,PointType>::ScalarT,panzer::Cell,PointType,panzer::Dim>
  RecombRate_Empirical_Defect<EvalT,Traits,PointType>::
getCoordField(Type<panzer::BASIS>, panzer::PureBasis const& basis, panzer::IntegrationRule const& /* ir */)
{
  const std::string coordName = panzer::GatherBasisCoordinates<EvalT,Traits>::fieldName(basis.name());
  return PHX::MDField<const ScalarT,panzer::Cell,PointType,panzer::Dim>(coordName, basis.coordinates);
}

///////////////////////////////////////////////////////////////////////////////
//
//  Constructor
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits, typename PointType>
RecombRate_Empirical_Defect<EvalT, Traits, PointType>::
RecombRate_Empirical_Defect(panzer::PureBasis const& basis,
                            panzer::IntegrationRule const& ir,
                            Teuchos::ParameterList const& p)
{
  using std::string;
  using Teuchos::RCP;
  using PHX::DataLayout;
  using PHX::MDField;
  using panzer::IntegrationRule;
  using panzer::BasisIRLayout;

  const charon::Names& n = *(p.get< RCP<const charon::Names> >("Names"));

  // Input from closure model factory

  //Bounding box for emitter base recombination region
  ebxlo = p.get<double>("eb x low");
  ebxhi = p.get<double>("eb x high");
  ebylo = p.get<double>("eb y low");
  ebyhi = p.get<double>("eb y high");
  ebzlo = p.get<double>("eb z low");
  ebzhi = p.get<double>("eb z high");

  //Get any user-overrid voltage values across the junctions
  ebOverrideBool = p.get<bool>("eb voltage override bool");
  cbOverrideBool = p.get<bool>("cb voltage override bool");

  ebVoltageOverride = p.get<double>("eb voltage override");
  cbVoltageOverride = p.get<double>("cb voltage override");

  //prefactors for the recombination rate in the emitter base region
  thermalVelocity = p.get<double>("thermal velocity");
  crossSection = p.get<double>("cross section");

  //intput data files for the emitter base pulse and mu data
  pulseType = p.get<std::string>("pulse type");
  muDataFile = p.get<std::string>("mu data file");

  //input number of delta pulses to resolve total pulse
  pulsePoints = p.get<int>("pulse resolution");

  //input number of delta pulses to resolve total pulse
  pulseIsRate = p.get<bool>("pulse is rate");

  // The pulse is specified via an input data file
  pulseDataFile = p.get<std::string>("pulse data file");

  // The specifies whether data in a pulse file is all points or a subset
  fileDataPulses = p.get<std::string>("file pulse sampling scheme");

  // Retrieve data layout
  RCP<DataLayout> scalar = p.get< RCP<DataLayout> >("Data Layout");
  num_points = scalar->dimension(1);

  // Get the damage spec out of the parameter list
  damage_spec = p.get<Teuchos::RCP<charon::PulseDamage_Spec> >("damage spec object");

  // Get the data for the empirical damage model
  damage_data = p.get<Teuchos::RCP<charon::EmpiricalDamage_Data> >("empirical damage data");

  //---------------------------------------------------------------------------------------
  //Perform some checks that we have everything we need for the empirical model to function
  //---------------------------------------------------------------------------------------

  if(!ebxlo && !ebxhi &&
     !ebylo && !ebyhi &&
     !ebzlo && !ebyhi)
  {
    std::stringstream msg;
    msg << "No bounding box parameters have been set for the emitter base empirical recombination.  Cannot Continute."<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  if (!crossSection)
  {
    std::stringstream msg;
    msg << "No cross-section has been set for the emitter base recombination.  Cannot continue."<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!crossSection, std::logic_error, msg.str());
  }

  if (!thermalVelocity)
  {
    std::stringstream msg;
    msg << "No thermal velocity has been set for the emitter base recombination.  Cannot continue."<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!crossSection, std::logic_error, msg.str());
  }

  if (muDataFile == "")
  {
    std::stringstream msg;
    msg << "No mu tables have been set for the emitter base recombination.  Cannot continue."<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  if (damage_spec->numberOfPoints() == 0)
  {
    std::stringstream msg;
    msg << "No pulse resolution specified for the emitter-base junction.  Need to specify the number of delta pulses to resolve total pulse waveform with the keyword \"pulse resolution.\""<< std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }


  //---------------------------------------------------------------------------------------
  // Create the object that computes the muNfp product
  //
  //---------------------------------------------------------------------------------------
  NfpMu = Teuchos::rcp(new charon::empiricalConvolution(*damage_spec, muDataFile, pulseIsRate));

  //---------------------------------------------------------------------------------------
  // fields
  //
  //---------------------------------------------------------------------------------------
  empirical_defect_rate = MDField<ScalarT,panzer::Cell,PointType>(n.field.empirical_defect_recomb,scalar);

  intrin_conc = MDField<const ScalarT,panzer::Cell,PointType>(n.field.intrin_conc,scalar);
  edensity = MDField<const ScalarT,panzer::Cell,PointType>(n.dof.edensity,scalar);
  hdensity = MDField<const ScalarT,panzer::Cell,PointType>(n.dof.hdensity,scalar);

  // scaling parameters
  scaleParams = p.get< RCP<charon::Scaling_Parameters> >("Scaling Parameters");
  t0 = scaleParams->scale_params.t0;
  C0 = scaleParams->scale_params.C0;
  X0 = scaleParams->scale_params.X0;

  coordinates = getCoordField(Type<PointType>(), basis, ir);

  // Evaluated fields
  this->addEvaluatedField(empirical_defect_rate);

  // Dependent fields
  this->addDependentField(intrin_conc);
  this->addDependentField(edensity);
  this->addDependentField(hdensity);
  this->addDependentField(coordinates);

  //Create an object to compute the SRH prefactor (pulse*mu)

  std::string name = "Empirical Defect Recombination";
  this->setName(name);
}


///////////////////////////////////////////////////////////////////////////////
//
//  evaluateFields()
//
///////////////////////////////////////////////////////////////////////////////
template <typename EvalT, typename Traits, typename PointType>
void RecombRate_Empirical_Defect<EvalT, Traits, PointType>::
evaluateFields(typename Traits::EvalData workset)
{
  using panzer::index_t;

  ScalarT scale_const = t0/C0;
  double tscale = Sacado::ScalarValue<ScalarT>::eval(t0);

  // get the current (present) time in [s]
  double currentTime = workset.time * tscale;

  //Get the timestep--N.B. Estimating the timestep here by 1/alpha.
  //This is potentially wrong.  But for some stupid reason, Rhythmos
  //won't just hand over the time step.

  //Also need to check if alpha=0.  This happens when currents are being calculated.
  double timestep = 0.0;
  if(workset.alpha != 0.0)
    timestep = tscale * 1/workset.alpha;

  double NfpMuProduct=0;

  // Calculate the applied bias across the emitter-base junction and
  // take the absolute value since the sign isn't relevant for the Mu
  // table. I think major breakage would occur if the emitter-base
  // weren't forward biased.
  // double voltage = (*appliedBiasAtContact)["base"] - (*appliedBiasAtContact)["emitter"];
  double v_emitter = damage_data->getVoltage("emitter");
  double v_base = damage_data->getVoltage("base");
  double voltage = v_base - v_emitter;
  if(ebOverrideBool)
    voltage = ebVoltageOverride;

  if (voltage < 0.0)
  {
    std::ostringstream err_msg;
    err_msg << "The base-emitter voltage is negative and "
            << "the empirical damage model is not currently "
            << "set up to handle that situation.";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, err_msg.str());
  }

  NfpMuProduct = NfpMu->computeNfpMu(currentTime, timestep, voltage);

  for (index_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (int point=0; point < num_points; ++point)
    {
      double xn = Sacado::ScalarValue<ScalarT>::eval(coordinates(cell,point,0));
      double yn = Sacado::ScalarValue<ScalarT>::eval(coordinates(cell,point,1));

      // If the point is within the bounding box evaluate the recombination
      if (xn > ebxlo && xn < ebxhi &&
          yn > ebylo && yn < ebyhi)
      {

        const ScalarT& n = edensity(cell,point)*C0;
        const ScalarT& p = hdensity(cell,point)*C0;
        const ScalarT& ni = intrin_conc(cell,point)*C0;

        empirical_defect_rate(cell,point) = thermalVelocity * crossSection *
          NfpMuProduct * (p * n - ni*ni)/(p + n + 2.0*ni) * scale_const;

      }
      else
      {
        empirical_defect_rate(cell,point) = 0.0;
      }
    }
  }

  return;
}

}

#endif
