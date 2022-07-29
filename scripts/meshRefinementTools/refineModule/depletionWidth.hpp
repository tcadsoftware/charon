

#ifndef DEPLETIONWIDTH_HPP
#define DEPLETIONWIDTH_HPP


class depletionWidth
{

private:

  //Define some constants

  const double boltzmann = 1.38064852e-23;  //  J/K

  const double kTq = 0.026;  //  V

  const double freePermittivity = 8.8541878176e-14;  // F/cm = C/Vcm

  const double fundamentalCharge = 1.6021766208e-19 ;  // C

  //Other necessary things
  double relativePermittivity;

  double intrinsicConcentration;  //  1/cm^3

public:

  depletionWidth();

  inline void setIntrinsicConcentration(double iC_) {intrinsicConcentration = iC_;}
  inline void setRelativePermittivity(double rP_) {relativePermittivity = rP_;}

  void calculateDepletionWidths(double relativePermittivity, double intrinsicConcentration,
				double acceptorConcentration, double donorConcentration,
				double& width, double& pWidth, double& nWidth);

  void calculateDepletionWidths(double acceptorConcentration, double donorConcentration, 
				double& width, double& pWidth, double& nWidth);

};


#endif

