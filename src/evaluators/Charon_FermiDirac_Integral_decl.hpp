
#include <stdexcept>

#ifndef CHARON_FERMIDIRACINTEGRAL_DECL_HPP
#define CHARON_FERMIDIRACINTEGRAL_DECL_HPP

namespace charon {

/**
 * @brief Abstract base class defining an interface for all Fermi-Dirac
 * integral algorithms.
 *
 * The developer can create their own Fermi-Dirac integral algorithm
 * class by inheriting this abstract base class.
 */
template <typename FDITemp>
class FermiDiracIntegralAlgorithm
{
public:

  /**
   * @brief this is the operator for invoking the algorithm.
   */
  virtual typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg) = 0;

  /**
   * @brief a pure virtual dtor for this ABC.
   */
  virtual ~FermiDiracIntegralAlgorithm() = 0;

  /**
   * @brief generally the class should have it's own built-in unit test
   * member function.
   */
  virtual void unitTest_() = 0;

  virtual std::string const& name() = 0;

protected:

  /**
   * @brief protected default ctor. Can't instantiate an object of this
   * class.
   */
  FermiDiracIntegralAlgorithm() {;}

};

/**
 * @brief Class for evaluating Fermi-Dirac integrals.
 *
 * This is a class that instantiates concrete objects in order to get
 * approximate values of the Fermi-Dirac integrals of order 1/2 and
 * -1/2. Note that for the purposes of efficiency the user should
 * instantiate an object of this class, of the desired Fermi-Dirac
 * integral type, and use it throughout the code as much as possible.
 */
template <typename FDITemp>
class FermiDiracIntegral
{
public:

  enum fermidirac_type
  {
    forward_PlusOneHalf,
    forward_MinusOneHalf,
    forward_AnyOrder,
    inverse_PlusOneHalf
  };

  FermiDiracIntegral(fermidirac_type integ_type,
                     std::string algorithm_name="",
                     double order=0.0 );

  void unitTest_();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  ~FermiDiracIntegral();

  std::string const& algorithmName();

private:

  // Can't instantiate using the default-argument version. Must supply
  // at least the order for the instance.
  FermiDiracIntegral();

  FermiDiracIntegralAlgorithm<FDITemp>* fd_algo;

  std::string const empty_name;

};

/**
 * @brief A concrete class using the approximation documented in the
 * Halen-Pulfrey paper to evaluate the Fermi-Dirac integral of order
 * 1/2.
 *
 * Implemented as documented in:
 *   P. Van Halen, D. L. Pulfrey, "Accurate, short series
 *   approximations to Fermi-Dirac integrals of order -1/2, 1/2, 3/2,
 *   2, 5/2, 3, and 7/2," J. of Appl. Phys., 57, 5271 (1985)
 * Note that there is an important erratum for the article as well in the
 * same journal, [59, 2264 (1986)].
 *
 * This may be more numerically intensive, but more accureate, than
 * other algorithms.
 */
template <typename FDITemp>
class HalenPulfrey_PlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  HalenPulfrey_PlusOneHalf_FIA();

  /**
   * @brief default method of invoking this algorithm
   */
  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  std::string const& name() {return algo_name;}

protected:

  /**
   * @brief invoke built-in unit test.
   */
  void unitTest_(){;}

private:
  double a1[7];
  double a2[7];
  double a3[7];
  double n1[7];

  std::string const algo_name;
};

/**
 * @brief A concreate class using the approximation documented in the
 * Aymerich-Humet paper to evaluate the Fermi-Dirac integral of order j.
 *
 * Implemented as documented in:
 *   X. Aymerich-Humet, F. Serra-Mestres, J. Millan, "A generalized
 *   approximation of the Fermi-Dirac integrals," J. of Appl. Phys., 54,
 *   2850 (May 1983)
 *
 * General approximation valid for any order of fermi-dirac integral,
 * but with increasing error for high values of j.
 */
template <typename FDITemp>
class Aymerich_AnyOrder_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  Aymerich_AnyOrder_FIA(double order=0.5);

  /**
   * @brief default method of invoking this algorithm
   */
  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  std::string const& name() {return algo_name;}

protected:

  /**
   * @brief invoke built-in unit test.
   */
  void unitTest_(){;}

private:
  double a, b, c;
  double fd_order;

  std::string const algo_name;
};

/**
 * @brief A concreate class using the approximation documented in the
 * Bednarczyk paper to evaluate the Fermi-Dirac integral of order 1/2.
 *
 * Implemented as documented in:
 *
 * D. Bednarczyk, J. Bednarczyk, "The Approximation of the Fermi-Dirac
 * Integral F_(1/2)(eta)," Physics Letters, 64A, No. 4, Jan. 1978
 */
template <typename FDITemp>
class Bednarczyk_PlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  Bednarczyk_PlusOneHalf_FIA();

  /**
   * @brief default method of invoking this algorithm
   */
  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  std::string const& name() {return algo_name;}

protected:

  /**
   * @brief invoke built-in unit test.
   */
  void unitTest_(){;}

private:
  double sqrtPi;
  double rFmult;

  std::string const algo_name;
};

/**
 * \brief A concrete class using the approximation documented in the
 * Halen-Pulfrey paper to evaluate the Fermi-Dirac integral of order
 * -1/2.
 *
 * Implemented as documented in:
 *   P. Van Halen, D. L. Pulfrey, "Accurate, short series
 *   approximations to Fermi-Dirac integrals of order -1/2, 1/2, 3/2,
 *   2, 5/2, 3, and 7/2," J. of Appl. Phys., 57, 5271 (1985)
 * Note that there is an important erratum for the article as well in the
 * same journal, [59, 2264 (1986)].
 *
 * This may be more numerically intensive, but more accureate, than
 * other algorithms.
 */
template <typename FDITemp>
class HalenPulfrey_MinusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  HalenPulfrey_MinusOneHalf_FIA();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  std::string const& name() {return algo_name;}

protected:
  void unitTest_(){;}

private:
  double a1[7];
  double a2[9];
  double a3[9];
  double n1[7];

  std::string const algo_name;
};

/**
 * @brief A concrete class used to approximate the inverse Fermi-Dirac
 * integral of order 1/2.
 *
 * Implemented as documented in:
 *    N.G. Nilsson, "An Accurate Approximation of the Generalized
 *    Einstein Relation for Degenerate Semiconductors,"
 *    Phys. Stat. Solidi (a), 19, pp. K75-K78, 1973
 */
template <typename FDITemp>
class Nilsson_InvPlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  Nilsson_InvPlusOneHalf_FIA();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  void unitTest_(){;}

  std::string const& name() {return algo_name;}

private:
  std::string const algo_name;
};


/**
 * @brief A concreate class used to approximate the inverse Fermi-Dirac
 * integral of order 1/2.
 *
 * Implemented as documented in:
 *    V.C. Aguilera-Navarro, G.A. Estevez, A. Kostecki, "A note on the
 *    Fermi-Dirac integral function," J. Appl. Physics, 63, 2848 (April
 *    1988).
 */
template <typename FDITemp>
class Aguilera_InvPlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  Aguilera_InvPlusOneHalf_FIA();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  void unitTest_(){;}

  std::string const& name() {return algo_name;}

private:
  double k[5];

  std::string const algo_name;
};

/**
 * @brief A concreate class used to approximate the inverse Fermi-Dirac
 * integral of order 1/2.
 *
 * Implemented as documented in:
 *   W.B. Joyce, R.W. Dixon, "Analytic approximations for the Fermi
 *   Energy of an ideal Fermi gas," Applied Physics Letters, 31, 354
 *   (Sept. 1977)
 *
 * This is the fastest (over the entire range) and most accurate of the
 * inverse Fermi algorithm.
 *
 * Note that there are refinements to the algorith as documented in the
 * original reference for arguments greater than 7.5 that originated
 * with a code written by S.M. Myers.
 */
template <typename FDITemp>
class JoyceDixonMyers_InvPlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  JoyceDixonMyers_InvPlusOneHalf_FIA();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  void unitTest_(){;}

  std::string const& name() {return algo_name;}

private:
  double a[8];
  double x10, x20;
  double y10, y20;
  double yp10, yp20;
  double c1, c2;

  std::string const algo_name;
};

/**
 * @brief A concreate class used to approximate the inverse Fermi-Dirac
 * integral of order 1/2.
 *
 * Implemented as documented in:
 *   W.B. Joyce, R.W. Dixon, "Analytic approximations for the Fermi
 *   Energy of an ideal Fermi gas," Applied Physics Letters, 31, 354
 *   (Sept. 1977)
 */
template <typename FDITemp>
class JoyceDixon_InvPlusOneHalf_FIA : public FermiDiracIntegralAlgorithm<FDITemp>
{
public:
  JoyceDixon_InvPlusOneHalf_FIA();

  typename FDITemp::ScalarT operator()(typename FDITemp::ScalarT arg);

  void unitTest_(){;}

  std::string const& name() {return algo_name;}

private:
  double a[8];
  double x10, x20;
  double y10, y20;
  double yp10, yp20;
  double c1, c2;

  std::string const algo_name;
};

} // END "namespace charon"

#endif // CHARON_FERMIINTEGRAL_DECL_HPP
