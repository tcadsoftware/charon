

#include <cmath>
#include "Charon_MMS_AnalyticFunctions.hpp"


namespace charon {

//**********************************************************************
// MMS_DD_RDH_1
//**********************************************************************
double MMS_DD_RDH_AnalyticFunction::potential(double const& x) const
{
  using std::atan;

  return -vl * (atan(beta*(x - length*0.5))/atan(beta*length*0.5));
}

//**********************************************************************
double MMS_DD_RDH_1_AnalyticFunction::edensity(double const& x) const
{
  return nhat*std::exp(potential(x)/lambda);
}

//**********************************************************************
double MMS_DD_RDH_1_AnalyticFunction::hdensity(double const& x) const
{
  return phat*std::exp(-potential(x)/lambda);
}


double MMS_DD_RDH_AnalyticFunction::dop(double const& x) const
{
  using std::atan;
  using std::pow;

  // obtain q & eps
  charon::PhysicalConstants const& cpc = charon::PhysicalConstants::Instance();
  double q = cpc.q;   // electron elemental charge in [C]
  double eps = cpc.eps0*11.8;

   return  (edensity(x) - hdensity(x)) - (2*pow(beta,3)*eps*vl*(-length*0.5 + x))/
  (q*pow(1 + pow(beta,2)*pow(-length*0.5 + x,2),2)*atan((beta*length)*0.5));

}

//**********************************************************************
std::vector<double> MMS_DD_RDH_AnalyticFunction::doping(double const& x) const
{
  std::vector<double> dopValue(2, 0.0);

  if (dop(x) < 0.0)
       {// negative doping = "Acceptor"
        dopValue[0] = std::abs(dop(x));
        dopValue[1] = 0.;
       }
  else
       {// positive doping = "Donor"
        dopValue[0] = 0.;
        dopValue[1] = dop(x);
       }

  return dopValue;
}

//**********************************************************************
// MMS_DD_RDH_2
//**********************************************************************
double MMS_DD_RDH_2_AnalyticFunction::edensity(double const& x) const
{
  using std::atan;
  using std::exp;

  return 7.69822e10*exp(-10.1075*atan(1.0e6*(-0.000025+x)));
}

//**********************************************************************
double MMS_DD_RDH_2_AnalyticFunction::hdensity(double const& x) const
{
  using std::atan;
  using std::exp;

  return 5.82759e10*exp(10.1075*atan(1.0e6*(-0.000025+x)));
}

//**********************************************************************
double MMS_DD_RDH_2_AnalyticFunction::bigN(double const& x) const
{
  return N0*exp(-4.0*x/length);
}

//**********************************************************************
double MMS_DD_RDH_2_AnalyticFunction::Nd(double const& x) const
{
  return 0.5*(bigN(x)+dop(x));
}

//**********************************************************************
double MMS_DD_RDH_2_AnalyticFunction::Na(double const& x) const
{
  return 0.5*(bigN(x)-dop(x));
}

}
