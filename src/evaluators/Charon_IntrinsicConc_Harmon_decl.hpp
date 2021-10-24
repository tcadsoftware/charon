
#ifndef CHARON_INTRINSICCONC_HARMON_DECL_HPP
#define CHARON_INTRINSICCONC_HARMON_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_FermiDirac_Integral.hpp"
#include "Charon_Scaling_Parameters.hpp"

using panzer::Cell;
using panzer::Point;

namespace charon {

template<typename EvalT, typename Traits>
class IntrinsicConc_Harmon
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    IntrinsicConc_Harmon(
      const Teuchos::ParameterList& p);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  // initialize the model parameters
  void initialize(const Teuchos::ParameterList& plist);

  // evaluate bgn values from a table read from a file
  void evaluateBGNFromFile(const ScalarT& dop, ScalarT& dEc, ScalarT& dEv);

  // find the lower index of a sorted vector for a given value
  int binarySearch(const ScalarT& dop);

  // output
  PHX::MDField<ScalarT,Cell,Point> nie;
  PHX::MDField<ScalarT,Cell,Point> effEg;
  PHX::MDField<ScalarT,Cell,Point> effChi;

  // input
  PHX::MDField<const ScalarT,Cell,Point> latt_temp;
  PHX::MDField<const ScalarT,Cell,Point> Eg;
  PHX::MDField<const ScalarT,Cell,Point> Chi;
  PHX::MDField<const ScalarT,Cell,Point> doping;
  PHX::MDField<const ScalarT,Cell,Point> elec_effdos;
  PHX::MDField<const ScalarT,Cell,Point> hole_effdos;

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0; // conc. scaling, [cm^-3]
  double T0; // temperature scaling, [K]

  int num_points;

  // material parameters
  double An, Ap;  // the A coefficient for n-type and p-type

  // include bgn when true
  bool includeBGN;

  // include the FD correction when true
  bool enableFD;

  // read BGN from a file when true
  bool bgnFromFile;

  // use a struct to store the doping-dependent bgn values
  struct dopBGNStruct
  {
    double dop;
    double dEc;
    double dEv;

    inline bool operator < (const dopBGNStruct &dbs) const
    { return ( dop < dbs.dop); }

    inline bool operator == (const dopBGNStruct &dbs) const
    { return ( dop == dbs.dop); }
  };

  std::vector<dopBGNStruct> dopBGNVec;

  Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

  Teuchos::RCP<charon::FermiDiracIntegral<EvalT> > inverseFermiIntegral;

}; // end of class IntrinsicConc_Harmon


}

#endif
