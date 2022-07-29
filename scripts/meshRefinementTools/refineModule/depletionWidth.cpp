

#include "depletionWidth.hpp"
#include <cmath>
#include <iostream>


depletionWidth::depletionWidth() :
  relativePermittivity(11.9),  //default to Si
  intrinsicConcentration(1e10)  //default to Si--can be approximate
{}




void depletionWidth::calculateDepletionWidths(double rP_, double iC_,
					      double acceptorConcentration, double donorConcentration,
					      double& width, double& pWidth, double& nWidth)
{

  relativePermittivity = rP_;
  intrinsicConcentration = iC_;
  calculateDepletionWidths(acceptorConcentration, donorConcentration, width, pWidth, nWidth);

}


void depletionWidth::calculateDepletionWidths(double acceptorConcentration, double donorConcentration, 
					      double& width, double& pWidth, double& nWidth)
{

  //Calculate the built in potential

  double builtInPotential = kTq * log(acceptorConcentration*donorConcentration 
				      / (intrinsicConcentration*intrinsicConcentration));

  //Calculate depletion region thickness

  double concentrationRatio = (acceptorConcentration + donorConcentration)
    / (acceptorConcentration*donorConcentration);

  width = sqrt(2.0*relativePermittivity*freePermittivity/fundamentalCharge
	       *(concentrationRatio)*builtInPotential);

  double byConc = 1.0/(acceptorConcentration + donorConcentration);

  nWidth = sqrt(2.0*relativePermittivity*freePermittivity/fundamentalCharge
		*(acceptorConcentration/donorConcentration)*byConc*builtInPotential);

  pWidth = sqrt(2.0*relativePermittivity*freePermittivity/fundamentalCharge
		*(donorConcentration/acceptorConcentration)*byConc*builtInPotential);

}
