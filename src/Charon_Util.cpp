#include "Charon_Util.hpp"

#include <fstream>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "Teuchos_TestForException.hpp"


namespace {
  class IsGreaterThan {
  public:
    IsGreaterThan(double x = 0.0) : val(x) {}
    bool operator() (double x) const {
      return (x > val);
    }
  private:
    double val;
  };


  class IsLessThan {
  public:
    IsLessThan(double x = 0.0) : val(x) {}
    bool operator() (double x) const {
      return (x < val);
    }
  private:
    double val;
  };
}

namespace charon {

  // a3*x*x*x + a2*x*x + a1*x + a0 = 0 (a3 != 0)
  void cubicsolve(double a3, double a2, double a1, double a0,
                  double& x1_re, double& x1_im,
                  double& x2_re, double& x2_im,
                  double& x3_re, double& x3_im) {
    // x*x*x + a1*x*x + a2*x + a3 = 0
    // x*x*x + a*x*x + b*x + c = 0
    assert(a3 != 0);
    double a = a2/a3;
    double b = a1/a3;
    double c = a0/a3;

    assert(a3 != 0);

    double Q = (3*b - a*a)/9;
    double R = (9*a*b - 27*c - 2*a*a*a)/54;
    double S = cbrt(R + sqrt(Q*Q*Q + R*R));
    double T = cbrt(R - sqrt(Q*Q*Q + R*R));
    double D = Q*Q*Q + R*R;

    if(D > 0) {
      // one real and two complex conjugate roots
      x1_re = S + T - a/3; x1_im = 0;
      x3_re = x2_re = -0.5*(S + T) - a/3;
      x2_im = 0.5*sqrt(3)*(S - T);
      x3_im = -x2_im;
    } else if(D == 0) {
      // all roots are real and at least two are equal
      S = T = cbrt(R);
      x1_im = x2_im = x3_im = 0;
      x1_re = S + T - a/3;
      x2_re = x3_re = -0.5*(S + T) - a/3;
    } else {
      // all roots are real and unequal
      double teta = acos(R/sqrt(-Q*Q*Q));
      x1_im = x2_im = x3_im = 0;
      x1_re = x2_re = x3_re = -a/3;
      double tmp = 2*sqrt(-Q);
      x1_re += tmp*cos(teta/3);
      x2_re += tmp*cos(teta/3 + 2*M_PI/3);
      x3_re += tmp*cos(teta/3 + 4*M_PI/3);
    }
  }



  // based on "A Note on the Solution of Quartic Equations" by Herbert E. Salzer
  // a4*x*x*x*x + a3*x*x*x + a2*x*x + a1*x + a0 = 0 (a4 != 0)
  void quarticsolve_salzer(double a4, double a3, double a2, double a1, double a0,
                           double& r1, double& r2,
                           double& r3, double& r4) {
    // x*x*x*x + a*x*x*x + b*x*x + c*x + d = 0
    assert(a4 != 0);
    double a = a3/a4;
    double b = a2/a4;
    double c = a1/a4;
    double d = a0/a4;
    r1 = 0, r2 = 0, r3 = 0, r4 = 0;

    // find the solution of the cubic equation
    // x*x*x - b*x*x + (a*c - 4*d)*x + d*(4*b - a*a) - c*c
    assert((d*(4*b-a*a)-c*c) != 0);

    double r1_re, r1_im, r2_re, r2_im, r3_re, r3_im;
    cubicsolve(1,-b,a*c-4*d,d*(4*b-a*a)-c*c,r1_re,r1_im,r2_re,r2_im,
               r3_re, r3_im);

    double x1 = r1_re;
    double m2 = 0.25*a*a - b + x1;
    double m = 0, n = 0;
    if(m2 > 0) {
      m = sqrt(m2);
      n = (a*x1 - 2*c)/4/m;
    } else if(m2 == 0) {
      m = 0;
      n = sqrt(0.25*x1*x1 -d);
    } else {
      // imaginary roots only
      r1 = r2 = r3 = r4 = 0;
    }
    double alpha = 0.5*a*a - x1 - b;
    double beta = 4*n - a*m;
    if(alpha + beta >= 0) {
      // one pair of real roots
      double gamma = sqrt(alpha + beta);
      r1 = -0.25*a + 0.5*m + 0.5*gamma;
      r2 = -0.25*a + 0.5*m - 0.5*gamma;
    }
    if(alpha-beta >= 0) {
      // another pair of real roots
      double delta = sqrt(alpha - beta);
      r3 = -0.25*a - 0.5*m + 0.5*delta;
      r4 = -0.25*a - 0.5*m - 0.5*delta;
    }
  }


  void findDopingPoints(const std::vector<double>& v, double dop,
                        std::pair<double,double>& points) {
    // the assumption is that v vector is ordered (ascending)
    if(std::find(v.begin(), v.end(), dop) != v.end()) {
      points.first = dop; points.second = dop;
    } else if(dop < v.front()) {
      points.first = v.front(); points.second = v.front();
    } else if(dop > v.back())  {
      points.first = v.back(); points.second = v.back();
    } else {
      std::vector<double>::const_iterator it_low =
        find_if(v.begin(), v.end(), IsLessThan(dop));
      std::vector<double>::const_iterator it_high =
        find_if(v.begin(), v.end(), IsGreaterThan(dop));
      points.first = *it_low; points.second = *it_high;
    }
  }


  double interpolateIonizEn(const std::map<double,double>& IonizEn,
                            const std::pair<double,double>& dop_bounds, double dop) {
    double en = 0;
    if(dop_bounds.first != dop_bounds.second) {
      // linearly interpolate the tabulated data on a loglin scale
      double dop1 = dop_bounds.first, dop2 = dop_bounds.second;
      const double& en1 = IonizEn.at(dop1);
      const double& en2 = IonizEn.at(dop2);
      double slope = (en2 - en1)/(log(dop2) - log(dop1));
      en = en1 + slope*(log(dop) - log(dop1));
    } else {
      en = IonizEn.at(dop_bounds.first);
    }
    return en;
  }


  void expandIonizEnParams(Teuchos::ParameterList& incmpl_ioniz) {
    const bool with_IonizAcc =
      (incmpl_ioniz.sublist("Acceptor").numParams() == 0) ? false : true;
    const bool with_IonizDon =
      (incmpl_ioniz.sublist("Donor").numParams() == 0) ? false : true;

    if(with_IonizAcc) {
      if(incmpl_ioniz.sublist("Acceptor").isParameter("Ionization Energy") &&
         incmpl_ioniz.sublist("Acceptor").isParameter("AccIncmplIoniz File") )
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "Either a value or a file must be specified for Acceptor Ionization Energy, not both" << std::endl);
      if(incmpl_ioniz.sublist("Acceptor").isParameter("AccIncmplIoniz File")) {
        // read ionization energy from file
        std::string ionizEnFile =
          incmpl_ioniz.sublist("Acceptor").get<std::string>("AccIncmplIoniz File");
        if(ionizEnFile == "")
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
             << "AccIncmplIoniz File name cannot be empty !" << std::endl);
        std::ifstream info(ionizEnFile.c_str());
        if(!info)
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
             << "Cannot read AccIncmplIoniz File '" << ionizEnFile << "'" << std::endl);
        Teuchos::RCP<std::vector<double> > accConc =
          Teuchos::rcp(new std::vector<double>);
        Teuchos::RCP<std::map<double,double> > accIonizEn =
          Teuchos::rcp(new std::map<double,double>);
        std::string line;
        while(std::getline(info,line)) {
          std::istringstream iss(line);
          double dop=0, ioniz_en=0;
          if(!(iss >> dop >> ioniz_en)) { break; }
          if(dop <= 0) dop = 1e-100;
          accConc->push_back(dop);
          accIonizEn->insert(std::pair<double,double>(dop,ioniz_en));
        }
        std::sort(accConc->begin(),accConc->end());
        incmpl_ioniz.sublist("Acceptor").sublist("AccIncmplIonizData").
          set("accConc",accConc);
        incmpl_ioniz.sublist("Acceptor").sublist("AccIncmplIonizData").
          set("accIonizEn",accIonizEn);
      }
    }

    if(with_IonizDon) {
      if(incmpl_ioniz.sublist("Donor").isParameter("Ionization Energy") &&
         incmpl_ioniz.sublist("Donor").isParameter("DonIncmplIoniz File") )
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
           << "Either a value or a file must be specified for Donor Ionization Energy, not both" << std::endl);
      if(incmpl_ioniz.sublist("Donor").isParameter("DonIncmplIoniz File")) {
        // read ionization energy from file
        std::string ionizEnFile =
          incmpl_ioniz.sublist("Donor").get<std::string>("DonIncmplIoniz File");
        if(ionizEnFile == "")
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
             << "DonIncmplIoniz File name cannot be empty !" << std::endl);
        std::ifstream info(ionizEnFile.c_str());
        if(!info)
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, std::endl
             << "Cannot read DonIncmplIoniz File '" << ionizEnFile << "'" << std::endl);
        Teuchos::RCP<std::vector<double> > donConc =
          Teuchos::rcp(new std::vector<double>);
        Teuchos::RCP<std::map<double,double> > donIonizEn =
          Teuchos::rcp(new std::map<double,double>);
        std::string line;
        while(std::getline(info,line)) {
          std::istringstream iss(line);
          double dop=0, ioniz_en=0;
          if(!(iss >> dop >> ioniz_en)) { break; }
          if(dop <= 0) dop = 1e-100;
          donConc->push_back(dop);
          donIonizEn->insert(std::pair<double,double>(dop,ioniz_en));
        }
        std::sort(donConc->begin(),donConc->end());
        incmpl_ioniz.sublist("Donor").sublist("DonIncmplIonizData").
          set("donConc",donConc);
        incmpl_ioniz.sublist("Donor").sublist("DonIncmplIonizData").
          set("donIonizEn",donIonizEn);
      }
    }
  }



}
