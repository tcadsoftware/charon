
#ifndef CHARON_BULKFIXCHARGE_FUNCTION_DECL_HPP
#define CHARON_BULKFIXCHARGE_FUNCTION_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Scaling_Parameters.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ParameterLibraryUtilities.hpp"
#include <vector>

using panzer::Cell;
using panzer::Point;

namespace charon {

  class uniformBulkFixQParams
  {
    public:
        void parseUniform (const Teuchos::ParameterList& plist);
        double value;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        bool varyingFixCharge;
  };


template<typename EvalT, typename Traits>
class BulkFixCharge_Function
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    BulkFixCharge_Function(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);

  private:

    using ScalarT = typename EvalT::ScalarT;

private:

  /**
   * @brief Evaluate total bulk fixed charge (due to various functions) at a given (x,y,z)
   */
  double evaluateBulkFixCharge(const double& x, const double& y, const double& z);

  /**
   * @brief Evaluate uniform/constant bulk fixed charge at given (x,y,z)
   */
  double evalUniformBulkFixQ(const double & x, const double & y,
    const double & z, const uniformBulkFixQParams & umfp);

  /**
   * @brief Evaluate other spatially dependent fixed charge profile 
   */
  
  // output fields
  PHX::MDField<ScalarT,Cell,Point> fixed_charge;    
  PHX::MDField<ScalarT,Cell,Point> fixed_charge_basis;

  // for IPs
  int int_rule_degree;
  std::size_t int_rule_index;
  int num_ips;
  int num_dims;
  int num_basis;

  // for basis points
  std::string basis_name;
  std::size_t basis_index;

  double chargeDensity;  // in [cm^{-3}]

  // scaling parameters
  Teuchos::RCP<charon::Scaling_Parameters> scaleParams;
  double C0;  // conc. scaling, [cm^-3]
 
  Teuchos::ParameterList bulkFixQParamList;

  std::vector<uniformBulkFixQParams> umfp_vec;

  //Homotopy stuff
  Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > user_value;
  bool varyingFixCharge;

}; // end of class BulkFixCharge_Function


}

#endif
