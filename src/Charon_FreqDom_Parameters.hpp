// FDparameter class
// its constructor, for this Panzer example, takes the input parameter list
// each of its objects has values/data:
//  truncation order
//  truncation scheme
//  fundamental harmonics
//  number of fundamental harmonics
//  (sorted possible) multi-indices (wrt the truncation scheme)
//  remapped fundamental harmonics 
//    (an array in the same ordering as the fundamental harmonics
//     e.g., fund = {1GHz, 1kHz, 2.3MHz} => remapped fund = {5, 2, 3})
//  number of time collocation points
//  
#ifndef CHARON_FREQDOM_PARAMETERS_HPP
#define CHARON_FREQDOM_PARAMETERS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include <cmath> // used to calculate sin(), cos(), and powers

class FreqDomParameters {

public:

    FreqDomParameters() = default;

    FreqDomParameters(Teuchos::RCP<Teuchos::ParameterList> freqdom_pl);

   ~FreqDomParameters(){ }

    void calc_NumFundamentalHarmonics();
    void calc_PossibleMultiIndices();
    void calc_TruncatedHarmonicBasis();
    void calc_RemappedFundamentalHarmonics();
    void calc_UnRemappedHarmonics();
    void calc_RemappedHarmonics();
    void calc_NumTotalHarmonics();

    void calc_NumTimeCollocationPoints();
    void calc_TimeCollocationPoints();

    void calc_TimeDomainDOF_InterwovenCoefficients();

    void calc_CosQuadratureWeightsForHarmonicNumber();
    void calc_SinQuadratureWeightsForHarmonicNumber();

    // get functions
    bool queryIsSmallSignal();
    Teuchos::Array<double> getHarmonicsVector();

    int getTruncationOrder();
    std::string getTruncationScheme();
    Teuchos::Array<double> getFundamentalHarmonics();
  
    int getNumFundamentalHarmonics();
    std::vector<std::vector<int> > getPossibleMultiIndices();
    Teuchos::RCP<std::vector<double> > getTruncatedHarmonicBasis();
    Teuchos::RCP<std::vector<double> > getRemappedFundamentalHarmonics(); 
    Teuchos::RCP<std::vector<double> > getUnRemappedHarmonics();
    Teuchos::RCP<std::vector<double> > getRemappedHarmonics();
    int getNumTotalHarmonics();  
    int getNumTimeCollocationPoints();
    std::vector<double> getTimeCollocationPoints();

    std::vector<std::vector<double> > getDofInterwovenCoeffsAtTimeCollPt();

    std::vector<std::vector<double> > getCosQuadratureWeightsForHarmonicNumber();
    std::vector<std::vector<double> > getSinQuadratureWeightsForHarmonicNumber();

private:

    bool isSmallSignal;

    int num_time_collocation_points;
    Teuchos::Array<double> harmonics_vector;

    // these values are set in the constructor
    int truncation_order;
    std::string truncation_scheme;
    double hybrid_exponent;
    Teuchos::Array<double> fundamental_harmonics;
    int specified_num_time_coll_points;
  
    // these values are calculated by a calc_ function
    std::vector<double> time_collocation_points;
    int num_fundamental_harmonics;
    std::vector<std::vector<int> > possible_multi_indices;
    Teuchos::RCP<std::vector<double> > truncated_harmonic_basis;
    Teuchos::RCP<std::vector<double> > remapped_fundamental_harmonics; 
    Teuchos::RCP<std::vector<double> > unremapped_harmonics;
    Teuchos::RCP<std::vector<double> > remapped_harmonics;
    int num_total_harmonics;

    std::vector<std::vector<double> > dof_interwoven_coeffs_at_time_coll_pt;
    std::vector<std::vector<double> > dof_cos_coeffs_at_time_coll_pt;
    std::vector<std::vector<double> > dof_sin_coeffs_at_time_coll_pt;

    std::vector<std::vector<double> > cos_quadrature_weights_for_harmonic_number;
    std::vector<std::vector<double> > sin_quadrature_weights_for_harmonic_number;

    bool printDebug = false;
};

#endif
