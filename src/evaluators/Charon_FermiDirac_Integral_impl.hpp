
#include <iostream>
#include <sstream>
#include <cmath>

#include "Charon_FermiDirac_Integral_decl.hpp"

#ifndef CHARON_FERMIDIRACINTEGRAL_IMPL_HPP
#define CHARON_FERMIDIRACINTEGRAL_IMPL_HPP

// Empty dtor required by ABC.
template <typename FDITemp>
charon::FermiDiracIntegralAlgorithm<FDITemp>::~FermiDiracIntegralAlgorithm()
{
}

// Default constructor not available for use
template <typename FDITemp>
charon::FermiDiracIntegral<FDITemp>::FermiDiracIntegral() :
  empty_name("NONE")
{
}

// Returns the name of the algorithm associated with this wrapper class
template <typename FDITemp>
std::string const& charon::FermiDiracIntegral<FDITemp>::algorithmName()
{
  if (fd_algo)
    return fd_algo->name();
  else
    return empty_name;
}

// The constructor that should be utilized for instantiating objects of this class
template <typename FDITemp>
charon::FermiDiracIntegral<FDITemp>::FermiDiracIntegral(enum fermidirac_type integ_type,
                                                        std::string algorithm_name,
                                                        double order)
{

  std::ostringstream err_msg;

  switch(integ_type)
  {
  case forward_PlusOneHalf:
    if (algorithm_name == "halenpulfrey")
    {
      fd_algo = new HalenPulfrey_PlusOneHalf_FIA<FDITemp>();
    }
    else if (algorithm_name == "bednarczyk" || algorithm_name == "") // DEFAULT
    {
      fd_algo = new Bednarczyk_PlusOneHalf_FIA<FDITemp>();
    }
    else
    {
      err_msg.str("");
      err_msg << "Unknown integral type \"" << algorithm_name << "\"requested "
              << "in charon::FermiDiracIntegral<FDITemp>::FermiDiracIntegral";
      throw std::runtime_error(err_msg.str());
    }
    break;

  case forward_MinusOneHalf:
    fd_algo = new HalenPulfrey_MinusOneHalf_FIA<FDITemp>();
    break;

  case inverse_PlusOneHalf:
    if (algorithm_name == "nilsson")
    {
      fd_algo = new Nilsson_InvPlusOneHalf_FIA<FDITemp>();
    }
    else if (algorithm_name == "aguilera")
    {
      fd_algo = new Aguilera_InvPlusOneHalf_FIA<FDITemp>();
    }
    else if (algorithm_name == "joycedixonmyers" || algorithm_name == "") // DEFAULT
    {
      fd_algo = new JoyceDixonMyers_InvPlusOneHalf_FIA<FDITemp>();
    }
    else if (algorithm_name == "joycedixon")
    {
      fd_algo = new JoyceDixon_InvPlusOneHalf_FIA<FDITemp>();
    }
    else
    {
      err_msg.str("");
      err_msg << "Unknown algorithm name for Fermi-Dirac integral";
      throw std::runtime_error(err_msg.str());
    }
    break;

  case forward_AnyOrder:
    fd_algo = new Aymerich_AnyOrder_FIA<FDITemp>(order);
    break;

  default:
    err_msg.str("");
    err_msg << "Unknown integral type requested in "
            << "charon::FermiDiracIntegral<FDITemp>::FermiDiracIntegral";
    throw std::runtime_error(err_msg.str());
  }

}

// Operator which simply invokes the specific algorithm operator
template <typename FDITemp>
typename FDITemp::ScalarT
charon::FermiDiracIntegral<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  return (*fd_algo)(arg);
}

// Unit-testing interface for algorithms. This should invoke the
// appropriate unit test member functions in the concrete algorithms
// class.
template <typename FDITemp>
void charon::FermiDiracIntegral<FDITemp>::unitTest_()
{
  std::cout << "Testing " << fd_algo->name() << ": " << std::endl;
  fd_algo->unitTest_();
}

// Destructor. Delete the allocated algorithm associated with the object
template <typename FDITemp>
charon::FermiDiracIntegral<FDITemp>::~FermiDiracIntegral()
{
  delete fd_algo;
}

//***************************************************************************
// Aymerich approximation for the Fermi-Dirac integral of order +1/2
template <typename FDITemp>
charon::Aymerich_AnyOrder_FIA<FDITemp>::Aymerich_AnyOrder_FIA(double order) :
  fd_order(order),
  algo_name("Aymerich F_j")
{
  a = std::sqrt(1.0 + (15.0/4.0)*(order+1.0) + (1.0/40.0)*(order+1.0)*(order+1.0));
  b = 1.8 + 0.61*order;
  c = 2.0 + (2.0 - std::sqrt(2.0))*std::pow(2.0, order);
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::Aymerich_AnyOrder_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  double v1 = (fd_order+1.0)*std::pow(2.0, fd_order+1.0);

  typename FDITemp::ScalarT v2 = std::pow(std::abs(arg-b), c);
  typename FDITemp::ScalarT v3 = std::pow(v2+std::pow(a,c), 1.0/c);

  typename FDITemp::ScalarT v4 = v1/std::pow(b+arg+v3, fd_order+1.0);

  typename FDITemp::ScalarT v5 = std::exp(-arg)/tgamma(fd_order+1.0);

  return (2.0/std::sqrt(M_PI))*(1.0/(v4+v5));
}

//***************************************************************************
// Halen-Pulfrey approximation for the Fermi-Dirac integral of order +1/2
template <typename FDITemp>
charon::HalenPulfrey_PlusOneHalf_FIA<FDITemp>::HalenPulfrey_PlusOneHalf_FIA() :
  algo_name("Halen-Pulfrey F_(+1/2)")
{
  a1[0] =  0.752253;
  a1[1] =  0.928195;
  a1[2] =  0.680839;
  a1[3] =  25.7829;
  a1[4] = -553.636;
  a1[5] =  3531.43;
  a1[6] = -3254.65;

  a2[0] =  0.765147;
  a2[1] =  0.604911;
  a2[2] =  0.189885;
  a2[3] =  0.020307;
  a2[4] = -0.004380;
  a2[5] = -0.000366;
  a2[6] =  0.000133;

  a3[0] =  0.777114;
  a3[1] =  0.581307;
  a3[2] =  0.206132;
  a3[3] =  0.017680;
  a3[4] = -0.006549;
  a3[5] =  0.000784;
  a3[6] = -0.000036;

  // For negative arguments
  n1[0] = 1.000000;
  n1[1] = 0.353568;
  n1[2] = 0.192439;
  n1[3] = 0.122973;
  n1[4] = 0.077134;
  n1[5] = 0.036228;
  n1[6] = 0.008346;
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::HalenPulfrey_PlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  typename FDITemp::ScalarT value = 0.0;

  if (arg <= 0.0)
  {
    for (int i=0; i < 7; i+=2)
    {
      value += n1[i]*std::exp((i+1.0)*arg);
    }
    for (int i=1; i < 7; i+=2)
    {
      value -= n1[i]*std::exp((i+1.0)*arg);
    }
  }
  else if (arg >= 4.0)
  {
    typename FDITemp::ScalarT mult = std::pow(arg, 1.5);
    for (int i=0; i < 7; ++i)
    {
      value += a1[i] / std::pow(arg, 2.0*i);
    }

    value = mult*value;
  }
  else if (arg <= 2.0)
  {
    for (int i=0; i < 7; ++i)
    {
      double pow_ = i;
      value += a2[i] * std::pow(arg, pow_);
    }
  }
  else // 2.0 < arg < 4.0
  {
    for (int i=0; i < 7; ++i)
    {
      double pow_ = i;
      value += a3[i] * std::pow(arg, pow_);
    }
  }

  return value;
}

//***************************************************************************
// Bednarczyk approximation for the Fermi-Dirac integral of order +1/2
template <typename FDITemp>
charon::Bednarczyk_PlusOneHalf_FIA<FDITemp>::Bednarczyk_PlusOneHalf_FIA() :
  algo_name("Bednarczyk F_(1/2)")
{
  sqrtPi = std::sqrt(M_PI);
  rFmult = 2.0/sqrtPi;
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::Bednarczyk_PlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  typename FDITemp::ScalarT aeta = std::pow(arg, 4.0) +
    33.6*arg*(1.0 - 0.68 * std::exp(-0.17*(arg+1.0)*(arg+1.0))) + 50.0;

  return rFmult*(0.5*sqrtPi*1.0/(0.75*sqrtPi*std::pow(aeta, -0.375)+std::exp(-arg)));
}

//***************************************************************************
// Halen-Pulfrey approximation for the Fermi-Dirac integral of order -1/2
template <typename FDITemp>
charon::HalenPulfrey_MinusOneHalf_FIA<FDITemp>::HalenPulfrey_MinusOneHalf_FIA() :
  algo_name("Halen-Pulfrey F_(-1/2)")
{
  a1[0] =  1.12837;
  a1[1] = -0.470698;
  a1[2] = -0.453108;
  a1[3] = -228.975;
  a1[4] =  8303.5;
  a1[5] = -118124.0;
  a1[6] =  632895.0;

  a2[0] =  0.604856;
  a2[1] =  0.380080;
  a2[2] =  0.059320;
  a2[3] = -0.014526;
  a2[4] = -0.004222;
  a2[5] =  0.001335;
  a2[6] =  0.000291;
  a2[7] = -0.000159;
  a2[8] =  0.000018;

  a3[0] =  0.638086;
  a3[1] =  0.292266;
  a3[2] =  0.159486;
  a3[3] = -0.077691;
  a3[4] =  0.018650;
  a3[5] = -0.002736;
  a3[6] =  0.000249;
  a3[7] = -0.000013;
  a3[8] =  2.9814e-07;

  // For negative arguments
  n1[0] = 0.999909;
  n1[1] = 0.706781;
  n1[2] = 0.572752;
  n1[3] = 0.466318;
  n1[4] = 0.324511;
  n1[5] = 0.152889;
  n1[6] = 0.033673;
}

template <typename FDITemp>
typename FDITemp::ScalarT charon::HalenPulfrey_MinusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  typename FDITemp::ScalarT value = 0.0;

  if (arg <= 0.0)
  {
    double sign = 1.0;
    double index;
    for (int i=0; i < 7; ++i)
    {
      index = i;
      value += sign*n1[i]*std::exp((index+1.0)*arg);
      sign *= -1.0;
    }
  }
  else if (arg >= 5.0)
  {
    typename FDITemp::ScalarT mult = std::sqrt(arg);
    for (int i=0; i < 7; ++i)
    {
      value += a1[i] / std::pow(arg, 2.0*i);
    }
    value = mult * value;
  }
  else if (arg < 2.5)
  {
    for (int i=0; i < 9; ++i)
    {
      double pow_ = i;
      value += a2[i] * std::pow(arg, pow_);
    }
  }
  else // 2.5 <= arg < 5.0
  {
    for (int i=0; i < 9; ++i)
    {
      double pow_ = i;
      value += a3[i] * std::pow(arg, pow_);
    }
  }

  return value;
}

//***************************************************************************
// Nilsson approximate inverse function of the Fermi-Dirac integral of
// order +1/2.
template <typename FDITemp>
charon::Nilsson_InvPlusOneHalf_FIA<FDITemp>::Nilsson_InvPlusOneHalf_FIA() :
  algo_name("Nilsson F^(-1)_(+1/2)")
{
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::Nilsson_InvPlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{

  typename FDITemp::ScalarT value;

  if (arg <= 0.0)
  {
    std::ostringstream err_msg;
    err_msg << "ERROR: An argument <= 0.0 for the Nilsson inverse Fermi "
            << "integral is not valid";
    throw std::runtime_error(err_msg.str());
  }
  else if (arg == 1.0) // singularity at 1.0, use the limit as arg->1.0
  {
    value = -0.5;
  }
  else
  {
    value = std::log(arg)/(1.0-arg*arg);
  }

  typename FDITemp::ScalarT v = std::pow(0.75*std::sqrt(M_PI)*arg, 2.0/3.0);
  return (value + v/(1.0+std::pow(0.24+1.08*v, -2.0)));
}

//***************************************************************************
// Aguilera approximate inverse function of the Fermi-Dirac integral of
// order +1/2.
template <typename FDITemp>
charon::Aguilera_InvPlusOneHalf_FIA<FDITemp>::Aguilera_InvPlusOneHalf_FIA() :
  algo_name("Aguilera F^(-1)_(+1/2)")
{
  k[0] = 4.8966851;
  k[1] = 3.3105795;
  k[2] = 73.6264033;
  k[3] = 0.1333760;
  k[4] = -21.0508644;
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::Aguilera_InvPlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{

  if (arg <= 0.0)
  {
    std::ostringstream err_msg;
    err_msg << "ERROR: An argument <= 0.0 for the Aguilera inverse Fermi "
            << "integral is not valid";
    throw std::runtime_error(err_msg.str());
  }

  return (std::log(arg) +
          k[0] * std::log(k[1]*arg + k[2]) +
          (k[3]*arg + k[4]));
}

//***************************************************************************
// Joyce-Dixon approximate inverse function of the Fermi-Dirac integral
// of order +1/2. There are refinements added by S.M. Myers to the
// algorithm for arguments greater than 7.5.
template <typename FDITemp>
charon::JoyceDixonMyers_InvPlusOneHalf_FIA<FDITemp>::JoyceDixonMyers_InvPlusOneHalf_FIA() :
  algo_name("Joyce-Dixon-Myers F^(-1)_(+1/2)")
{
  a[0] = 1.0/std::sqrt(8.0);
  a[1] = 0.1875 - std::sqrt(3.0)/9.0;
  a[2] = 0.125 + (5.0/48.0)*std::sqrt(2.0) - std::sqrt(6.0)/9.0;
  a[3] = 1585.0/6912.0 + (5.0/32.0)*std::sqrt(2.0) - (5.0/24.0)*std::sqrt(3.0)
    - 0.04*std::sqrt(5.0);

  a[5] = 4.0/3.0;
  a[4] = std::pow(0.75*std::sqrt(M_PI), a[5]);
  a[6] = M_PI*M_PI/6.0;
  a[7] = 1.0/3.0;

  x10 = 7.5;

  double tmp1 = x10*x10;
  y10 = std::log(x10) + a[0] * x10 + a[1] * (tmp1) + a[2]*(tmp1*x10) + a[3]*(tmp1 * tmp1);
  yp10 = 1.0/x10 + a[0] + a[1] * 2.0 * x10 + a[2] * 3.0 * (x10 * x10) + a[3] * 4.0 * (x10 * x10 * x10);

  x20 = 8.5;

  y20 = std::sqrt(a[4] * std::pow(x20, a[5]) - a[6]);
  yp20 = 0.5 / std::sqrt(a[4] * pow(x20, a[5]) - a[6]) * a[5] * a[4] * pow(x20, a[7]);

  double delx = 0.5;
  double dely = y20 - y10;
  c1 = dely * 0.5 / (delx * delx) - yp10 * 0.75 / delx - yp20 * 0.25 / delx;
  c2 = dely * 0.5 / (delx * delx) - yp20 * 0.75 / delx - yp10 * 0.25 / delx;
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::JoyceDixonMyers_InvPlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  typename FDITemp::ScalarT result;

  if (arg <= 0.0)
  {
    std::ostringstream err_msg;
    err_msg << "ERROR: An argument <= 0.0 for the Joyce-Dixon inverse Fermi "
            << "integral is not valid";
    throw std::runtime_error(err_msg.str());
  }
  else if (arg <= 7.5)
  {
    typename FDITemp::ScalarT temp_val = arg*arg;

    result = std::log(arg) + a[0]*arg + a[1]*temp_val +
      a[2]*temp_val*arg + a[3]*temp_val*temp_val;
  }
  // Added as found in code from Myers.
  else if (arg > 7.5 && arg <= 8.0)
  {
    typename FDITemp::ScalarT diff = arg - 7.5;
    result = y10 + yp10*diff + c1*(diff*diff);
  }
  else if (arg > 8.0 && arg < 8.5)
  {
    typename FDITemp::ScalarT diff = 8.5 - arg;
    result = y20 - yp20*diff - c2*diff*diff;
  }
  else // arg >= 8.5
  {
    result = std::sqrt(a[4] * std::pow(arg, a[5]) - a[6]);
  }

  return result;
}

//***************************************************************************
// Joyce-Dixon approximate inverse function of the Fermi-Dirac integral
// of order +1/2. There are refinements added by S.M. Myers to the
// algorithm for arguments greater than 7.5.
template <typename FDITemp>
charon::JoyceDixon_InvPlusOneHalf_FIA<FDITemp>::JoyceDixon_InvPlusOneHalf_FIA() :
  algo_name("Joyce-Dixon F^(-1)_(+1/2)")
{
  a[0] = 1.0/std::sqrt(8.0);
  a[1] = 0.1875 - std::sqrt(3.0)/9.0;
  a[2] = 0.125 + (5.0/48.0)*std::sqrt(2.0) - std::sqrt(6.0)/9.0;
  a[3] = 1585.0/6912.0 + (5.0/32.0)*std::sqrt(2.0) - (5.0/24.0)*std::sqrt(3.0)
    - 0.04*std::sqrt(5.0);
}

template <typename FDITemp>
typename FDITemp::ScalarT
charon::JoyceDixon_InvPlusOneHalf_FIA<FDITemp>::operator()(typename FDITemp::ScalarT arg)
{
  typename FDITemp::ScalarT result;

  if (arg <= 0.0)
  {
    std::ostringstream err_msg;
    err_msg << "ERROR: An argument <= 0.0 for the Joyce-Dixon inverse Fermi "
            << "integral is not valid";
    throw std::runtime_error(err_msg.str());
  }

  typename FDITemp::ScalarT temp_val = arg*arg;

  result = std::log(arg) + a[0]*arg + a[1]*temp_val +
    a[2]*temp_val*arg + a[3]*temp_val*temp_val;

  return result;
}

#endif // CHARON_FERMIDIRACINTEGRAL_IMPL_HPP
