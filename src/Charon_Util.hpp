#ifndef CHARON_UTIL_HPP
#define CHARON_UTIL_HPP

#include <vector>
#include <map>
#include "Teuchos_ParameterList.hpp"

namespace charon {

  // a3*x*x*x + a2*x*x + a1*x + a0 = 0 (a3 != 0)
  void cubicsolve(double a3, double a2, double a1,
                  double a0, double& x1_re, double& x1_im,
                  double& x2_re, double& x2_im,
                  double& x3_re, double& x3_im);

  // a4*x*x*x*x + a3*x*x*x + a2*x*x + a1*x + a0 = 0 (a4 != 0)
  void quarticsolve_salzer(double a4, double a3, double a2,
                           double a1, double a0, double& r1,
                           double& r2, double& r3, double& r4);

  void findDopingPoints(const std::vector<double>& v, double dop,
                        std::pair<double,double>& points);

  double interpolateIonizEn(const std::map<double,double>& IonizEn,
                            const std::pair<double,double>& dop_bounds,
                            double dop);

  void expandIonizEnParams(Teuchos::ParameterList& incmpl_ioniz);
}


#endif
