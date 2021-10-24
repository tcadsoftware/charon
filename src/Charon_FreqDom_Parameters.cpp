#include "Charon_FreqDom_Parameters.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include <set>

///////////////////////////////////////////////////////////////
// depending on where a FreqDomParameters object is instantiated, we may need to compute more or fewer internal paramers
FreqDomParameters::FreqDomParameters(Teuchos::RCP<Teuchos::ParameterList> freqdom_pl){

  this->isSmallSignal = freqdom_pl->get<bool>("Enable Small Signal Analysis");
  std::cout << "We are performing a ";
  if(this->isSmallSignal) std::cout << "SMALL SIGNAL ANALYSIS." << std::endl;
  if(!(this->isSmallSignal)) std::cout << "LARGE SIGNAL ANALYSIS." << std::endl;

  // when in doubt of the order in which calc_ functions can be called,
  // refer to the order of the funcation declarations

  this->truncation_order = freqdom_pl->get<int>("Truncation Order");
  this->truncation_scheme = freqdom_pl->get<std::string>("Truncation Scheme");
  this->hybrid_exponent = freqdom_pl->get<double>("Hybrid Exponent");
  this->fundamental_harmonics = freqdom_pl->get<Teuchos::Array<double> >("Fundamental Harmonics");
  this->specified_num_time_coll_points = freqdom_pl->get<int>("Number of Time Collocation Points");
  this->calc_NumFundamentalHarmonics();

  if(!freqdom_pl->isParameter("Remapped Fundamental Harmonics")){
    std::cout << "We were not provided with remapped harmonics, so we calculate them!" << std::endl;
    // these can only be performed after calc_PossibleMultiIndices();
    this->calc_RemappedFundamentalHarmonics(); // perform a frequency remapping
    this->calc_PossibleMultiIndices();
    this->calc_UnRemappedHarmonics();
    this->calc_RemappedHarmonics();
  }

  else if(freqdom_pl->isParameter("Remapped Fundamental Harmonics")){
    std::cout << "We were provided with remapped harmonics, so we're using them!" << std::endl;
    this->remapped_fundamental_harmonics = Teuchos::rcp(new std::vector<double>());
    for(auto harmonic : freqdom_pl->get<Teuchos::Array<double> >("Remapped Fundamental Harmonics"))
      this->remapped_fundamental_harmonics->emplace_back(harmonic);
    this->calc_PossibleMultiIndices();
    this->calc_UnRemappedHarmonics();
    // skip calculating the remapped fundamental harmonics via calc_RemappedFundamentalHarmonics(),
    // since they're provided these can only be performed after calc_PossibleMultiIndices();
    this->calc_RemappedHarmonics(); // uses remapped fundamental harmonics
  }

  // these can only be performed after calc_PossibleMultiIndices();
  this->calc_TruncatedHarmonicBasis();

  // this can only be performed after calc_RemappedHarmonics();
  this->calc_NumTotalHarmonics(); 
  std::cout << "The total num of harmonics is: " << std::to_string(this->num_total_harmonics) << std::endl;

  // can only be performed after calc_RemappedFundamentalHarmonics
  this->calc_NumTimeCollocationPoints();
  this->calc_TimeCollocationPoints();

  // can only be performed after calc_RemappedHarmonics, calc_Time_CollocationPoint
  this->calc_TimeDomainDOF_InterwovenCoefficients();
  this->calc_CosQuadratureWeightsForHarmonicNumber();
  this->calc_SinQuadratureWeightsForHarmonicNumber();

}


///////////////////////////////////////////////////////////////
// calculate the number of fundamental harmonics from 
void FreqDomParameters::calc_NumFundamentalHarmonics(){ 
  this->num_fundamental_harmonics = fundamental_harmonics.size();
}


// take the fundamental harmonics array and determine the possible multi-indices
// which are the coefficients of the fundamental harmonics which form 
// linear combinations lying within the bounds of the truncation scheme
void FreqDomParameters::calc_PossibleMultiIndices(){ 

  // Description of this function
  // (1) create possible index values vector, 
  // (2) create possible_multi_index_vector, 
  // (3) slice the multi_index vector by the truncation scheme, obtaining the (admissible) multi_index vector
  // (4) order the multi_index_vector with respect to the fundamental_harmonics

  // (1) create a vector of possible index values
  std::vector<int> possible_index_values;
  possible_index_values.push_back(0);
  for (int j = 0 ; j < (this->truncation_order) ; j++){
    possible_index_values.push_back(j+1);
    possible_index_values.push_back(-j-1);
  }
  /*
  std::cout << "The list of possible index values are: ";
  for(auto number : possible_index_values) std::cout << number << " ";
  std::cout << std::endl;
  */

  // (2) create the possible_multi_index vector
  // this includes all vectors in a Box truncation scheme
  std::vector<std::vector<int> > poss_multi_indices;
  std::vector<int> partial_multi_index;
  for (auto value : possible_index_values){    
    partial_multi_index.clear();
    partial_multi_index.push_back(value);
    poss_multi_indices.push_back(partial_multi_index);
  }

  std::vector<int> temp_partial_multi_index;
  for (int dim = 1 ; dim < (this->num_fundamental_harmonics) ; dim++){

    unsigned num_old_vectors = poss_multi_indices.size();
    for (unsigned row = 0 ; row < num_old_vectors ; row++){
      for (unsigned nz_entry = 1 ; nz_entry < possible_index_values.size() ; nz_entry++){
        temp_partial_multi_index.clear();
        temp_partial_multi_index = poss_multi_indices[row];
        temp_partial_multi_index.push_back(possible_index_values[nz_entry]);
        poss_multi_indices.push_back(temp_partial_multi_index);
      }
      poss_multi_indices[row].push_back(0);
    }
  }

  // print to check
  /*
  for(int row = 0 ; row < poss_multi_indices.size() ; row++){
    std::cout << "Multi index number " << std::to_string(row) << " is: ";
    for(auto entry : poss_multi_indices[row]){
      std::cout << std::to_string(entry) << " ";
    }
    std::cout << std::endl;
  }
  */

  // (3) slicing away the non-admissible multi-indices according to the truncation scheme
  // essentially, the other truncation schemes provide a different row_norm to compare to the cut_off

  for(auto row = poss_multi_indices.begin() ; row != poss_multi_indices.end() ; ){
    double row_norm = 0;
    double cut_off = this->truncation_order;
    double linear_combination = 0;

    // this computation of the row norm is for the Diamond truncation scheme
    if(this->truncation_scheme == "Diamond"){
      for(auto entry : *row) 
        row_norm = row_norm + abs(entry);
    }

    // this computation of the row norm is for the Box truncation scheme
    if(this->truncation_scheme == "Box"){
      // add nothing, since we admit every multi-index
    }

    // this computation of the row norm is for the Box truncation scheme
    if(this->truncation_scheme == "Hybrid"){
      for(auto entry : *row) 
        row_norm = row_norm + std::pow(abs(entry), this->hybrid_exponent);
    }

    for(int i = 0 ; i < (this->num_fundamental_harmonics) ; i++)
      linear_combination += (*row)[i]*(*(this->remapped_fundamental_harmonics))[i];

    if(row_norm > cut_off || linear_combination <0)
      poss_multi_indices.erase(row);
    else ++row;
  }

  // print to check
  /*
  for(int row = 0 ; row < poss_multi_indices.size() ; row++){
    double row_norm = 0;
    for(auto entry : poss_multi_indices[row]){
        row_norm = row_norm + abs(entry);
    }
    std::cout << "With multi-index norm = " << std::to_string(row_norm) << ", we admit: (";
    for(auto entry : poss_multi_indices[row]){
        std::cout << std::to_string(entry) << ",";
    }
    std::cout << ") and it is " << std::to_string(row_norm <= cut_off) << " row_norm <= cut_off." << std::endl;
  }
  */

  // (4) order the admissible multi-indices according to the fundamental harmonics
  // if the inputted frequencies are already in increasing order, then this is simply the left-to-right lexicographical ordering
  // which can be achieved by a simple sort(poss_multi_indices.begin(), poss_multi_indices.end());
  // but this requires the user to input in increasing order 
  // if we order the fundamental frequencies after input, book-keeping needs to be done to re-interpret each multi-index

  /* using sort(), but need to pass a parameter (fundamental harmonics) into the comparator function
  // using struct to implement operator() approach

  // create a comparison function which depends on the fundamental harmonics  
  struct compareFrequencyLinComb {

    compareFrequencyLinComb(Teuchos::Array<double> fundamental_harmonics) {
      this->fundamental_harmonics = fundamental_harmonics; }

    bool operator () (std::vector<int> multi_index_a, std::vector<int> multi_index_b) {

      double linear_combination = 0;
      for(int i = 0 ; i < fundamental_harmonics.size() ; i++) 
        linear_combination += (multi_index_a[i]-multi_index_b[i])*fundamental_harmonics[i];
      return linear_combination < 0;
      }

    Teuchos::Array<double> fundamental_harmonics;
  };

  sort(poss_multi_indices.begin(), poss_multi_indices.end(), compareFrequencyLinComb(fundamental_harmonics));
  */

  // lambda expression approach (cleaner?)
  sort(poss_multi_indices.begin(), poss_multi_indices.end(),
      [&](std::vector<int> multi_index_a, std::vector<int> multi_index_b)-> bool {
            double linear_combination = 0;
            for(int i = 0 ; i < (this->fundamental_harmonics).size() ; i++)
              linear_combination += 
                (multi_index_a[i]-multi_index_b[i])*(*(this->remapped_fundamental_harmonics))[i];
            return linear_combination < 0;
         }
      );

  /*
  for(int row = 0 ; row < poss_multi_indices.size() ; row++){
    std::cout << "The row number " << std::to_string(row) << " multi-index is: ";
    for(auto entry : poss_multi_indices[row]){
      std::cout << std::to_string(entry) << " ";
    }
    std::cout << std::endl;
  }
  */

  this->possible_multi_indices = poss_multi_indices;
}

// the following calculates: truncated_harmonic_basis
// the following requires: possible_muti_indices, fundamental_harmonics
void FreqDomParameters::calc_TruncatedHarmonicBasis(){

  Teuchos::RCP<std::vector<double> > trunc_harmonic_basis = Teuchos::rcp(new std::vector<double>);
  for(auto row = (this->possible_multi_indices).begin() ; row != (this->possible_multi_indices).end() ; row++){
    // this computation of the row norm is for the Diamond truncation scheme
    double linear_combination = 0;
    for(int i = 0 ; i < (this->num_fundamental_harmonics) ; i++)
      linear_combination += (*row)[i]*(this->fundamental_harmonics)[i];
    //std::cout << "The truncated frequency basis contains: " << std::to_string(linear_combination) << std::endl;
    trunc_harmonic_basis->push_back(linear_combination);
  }
  //std::cout << "This concludes the truncated frequency basis." << std::endl;

  this->truncated_harmonic_basis = trunc_harmonic_basis;
}


// take the array of fundamental harmonics and return remapped fundamental harmonics
// this exploits p-adic expansions. the idea is, given the fundamental frequencies vector (\omega_1, \ldots , \omega_M),
// the linear combinations \sum_i (k_i \omega_i) are distinct frequencies; a remapping of \omega_i to \eta_i
// should form a bijection onto the linear combinations \sum_i (k_i \eta_i)
// requires: fundamental_harmonics, num_fundamental_harmonics
void FreqDomParameters::calc_RemappedFundamentalHarmonics(){

  //std::cout << "The remapped fundamental harmonics are among: " ;
  std::vector<double> target(this->num_fundamental_harmonics, 0);

  for(int k = 0 ; k < (this->num_fundamental_harmonics) ; k++){
    if(this->truncation_scheme == "Box") 
      target[k] = (double)std::pow(2 * (this->truncation_order) + 1 , k);
    if(this->truncation_scheme == "Diamond")
      // TODO: this isn't very bijective...
      target[k] = (double)((this->num_fundamental_harmonics) * k + std::pow(this->truncation_order , (double)k) + 1);
    if(this->truncation_scheme == "Hybrid")
      // TODO: for now, just use the Box remapping (guaranteed to be bijective, but not densely packing) 
      target[k] = (double)std::pow(2 * (this->truncation_order) + 1, k);
      //std::cout << std::to_string(target[k]) << " ";
  }
  //std::cout << std::endl;

  std::vector<std::vector<double> > 
    harmonics_and_indices(this->num_fundamental_harmonics, std::vector<double>(2) );

  for(int k = 0 ; k < (this->num_fundamental_harmonics) ; k++)
    harmonics_and_indices[k] = {(this->fundamental_harmonics)[k] , (double)k};

  //std::cout << "The read-in fundamental frequency list is: " << std::endl;
  //for(unsigned k = 0 ; k < harmonics_and_indices.size() ; k++)
  //    std::cout << std::to_string(harmonics_and_indices[k][0]) << std::endl;
  // end print

  // sorting harmonics_and_indices lexicographically will store a permutation in the second column
  sort(harmonics_and_indices.begin() , harmonics_and_indices.end());

  // use the inverse of the permutation stored in the second column of the sorted harmonics_and_indices
  Teuchos::RCP<std::vector<double> > 
    remapped_fund_harmonics = Teuchos::rcp(
      new std::vector<double>(this->num_fundamental_harmonics, 0) );
  for(int k = 0 ; k < (this->num_fundamental_harmonics) ; k++)
      (*remapped_fund_harmonics)[ harmonics_and_indices[k][1] ] = target[k];

  // print the remapped frequency list (verify that it is in the same ordering as the read-in list)
  for(auto freq = (*remapped_fund_harmonics).begin() ; freq != (*remapped_fund_harmonics).end() ; freq++)
      std::cout << "The remapped fundamental frequency list contains: " << std::to_string(*freq) << std::endl;

  this->remapped_fundamental_harmonics = remapped_fund_harmonics;
}

// calculate the un-remapped harmonics from the multi-indices
void FreqDomParameters::calc_UnRemappedHarmonics(){

  // (5) use those admissible indices (a vector of coefficients) to calculate all the harmonics of the truncated_harmonic_basis
  // computed by looping over the admissible indices set, dotting each multi-index with the fundamental harmonics vector
  // doing this in the lexicographical order should ensure the resulting list is in the real order

  Teuchos::RCP<std::vector<double> > harmonics = 
    Teuchos::rcp(new std::vector<double>( (this->possible_multi_indices).size() ) );

  for(unsigned k = 0 ; k < (this->possible_multi_indices).size() ; k++){
    double lin_comb = 0;
    for(int l = 0 ; l < (this->num_fundamental_harmonics) ; l++)
      lin_comb += (this->fundamental_harmonics)[l] 
                  * (this->possible_multi_indices)[k][l];

    (*harmonics)[k] = lin_comb;

    // print un-remapped frequencies list
    std::cout << "The un-remapped harmonic with a multi-index of: (";
    std::cout << std::to_string((this->possible_multi_indices)[k][0]);
    for(auto i = 1 ; i < (this->num_fundamental_harmonics) ; i++){
      std::cout << "," << std::to_string((this->possible_multi_indices)[k][i]);
    }
    std::cout << ") has value " << std::to_string(lin_comb) << std::endl;
  }

  // dont simply do 'this->unremapped_harmonics = harmonics;'
  // because there may be duplicated harmonics in the truncation scheme!
  // for example, if omega = {2 Hz,3 Hz}, then (2,-1) and (-1,2) correspond to 1 Hz

  // remove duplicated harmonics
  std::set<double> set_of_harmonics;
  unsigned size = harmonics->size();
  for( unsigned i = 0; i < size; ++i )
    set_of_harmonics.insert( (*(harmonics))[i] );

  harmonics->assign( set_of_harmonics.begin(), set_of_harmonics.end() );
  this->unremapped_harmonics = Teuchos::rcp(new std::vector<double>());
  for(auto h : *harmonics)
    this->unremapped_harmonics->emplace_back(double(h));

  std::cout << "Thus, the un-remapped harmonics we consider are: (";
  for(auto i : *(this->unremapped_harmonics))
    std::cout << std::to_string(i) << ",";
  std::cout << ")" << std::endl;
}

// calculate the remapped harmonics from the multi-indices
void FreqDomParameters::calc_RemappedHarmonics(){

  // (5) use those admissible indices (a vector of coefficients) to calculate all the harmonics of the truncated_harmonic_basis
  // computed by looping over the admissible indices set, dotting each multi-index with the fundamental harmonics vector
  // doing this in the lexicographical order should ensure the resulting list is in the real order

  Teuchos::RCP<std::vector<double> > harmonics = 
    Teuchos::rcp(new std::vector<double>( (this->possible_multi_indices).size() ) );

  for(unsigned k = 0 ; k < (this->possible_multi_indices).size() ; k++){
    double lin_comb = 0;
    for(int l = 0 ; l < (this->num_fundamental_harmonics) ; l++)
      lin_comb += (*(this->remapped_fundamental_harmonics))[l] 
                  * (this->possible_multi_indices)[k][l];

    (*harmonics)[k] = lin_comb;

    // print remapped frequencies list
    std::cout << "The remapped harmonic with a multi-index of: (";
    std::cout << std::to_string((this->possible_multi_indices)[k][0]);
    for(auto i = 1 ; i < (this->num_fundamental_harmonics) ; i++){
      std::cout << "," << std::to_string((this->possible_multi_indices)[k][i]);
    }
    std::cout << ") has value " << std::to_string(lin_comb) << std::endl;
  }

  // dont simply do 'this->remapped_harmonics = harmonics;'
  // because there may be duplicated harmonics in the truncation scheme!
  // for example, if omega = {2 Hz,3 Hz}, then (2,-1) and (-1,2) correspond to 1 Hz

  // remove duplicated harmonics
  std::set<double> set_of_harmonics;
  unsigned size = harmonics->size();
  for( unsigned i = 0; i < size; ++i )
    set_of_harmonics.insert( (*(harmonics))[i] );

  harmonics->assign( set_of_harmonics.begin(), set_of_harmonics.end() );
  this->remapped_harmonics = Teuchos::rcp(new std::vector<double>());
  for(auto h : *harmonics)
    this->remapped_harmonics->emplace_back(double(h));

  std::cout << "Thus, the harmonics we consider are: (";
  for(auto i : *(this->remapped_harmonics))
    std::cout << std::to_string(i) << ",";
  std::cout << ")" << std::endl;
}

// calculate the total number of harmonics (remapped or not)
void FreqDomParameters::calc_NumTotalHarmonics(){ 
  this->num_total_harmonics = (this->remapped_harmonics)->size();
}

// calculate the number of time collocaiton points
void FreqDomParameters::calc_NumTimeCollocationPoints(){

  // small signal
  // time is essentially an index for the harmonic "block" index
  if(this->isSmallSignal)
    this->num_time_collocation_points = 2*(this->num_total_harmonics);

  // large signal
  if(!(this->isSmallSignal)) 
    this->num_time_collocation_points = 2 * (*(std::max_element( (*(this->remapped_harmonics)).begin(), (*(this->remapped_harmonics)).end() ))) + 3;

  if(this->specified_num_time_coll_points != 0){
    if(this->isSmallSignal)
      std::cout << "The number of time collocation points can't be manually specified for a small signal analysis." << std::endl;
    if(!(this->isSmallSignal)){
      std::cout << "The number of time collocation points was manually specified for a large signal analysis." << std::endl;
      if(this->specified_num_time_coll_points < this->num_time_collocation_points)
        std::cout << "However, the specified number is less than the Nyquist Sampling Theorem requirement. Defaulting to the NST minimum, instead." << std::endl;
      if(!(this->specified_num_time_coll_points < this->num_time_collocation_points)){
        std::cout << "It is specified to be at least that required by the Nyquist Sampling Theorem." << std::endl;
        this->num_time_collocation_points = this->specified_num_time_coll_points;
      }
    }
  }
  std::cout << "The number of time collocation points is: " << std::to_string(this->num_time_collocation_points) << std::endl;
}

// calculate the time collocation points themselves
void FreqDomParameters::calc_TimeCollocationPoints(){

  std::vector<double> points(this->num_time_collocation_points);

  for(double t = 0 ; t < (this->num_time_collocation_points) ; t++)
    points[t] = t/((double)(this->num_time_collocation_points - 1));	

  // print the time collocation points to check
  //std::cout << "The time collocation points are: " << std::endl;
  //for(int t = 0 ; t < (this->num_time_collocation_points) ; t++){
  //  std::cout << std::to_string(points[t]) << " ";
  //}
  //std::cout << std::endl;

  this->time_collocation_points = points;
}

// calculate the time domain dof (at the tiume collocation point) coefficients
void FreqDomParameters::calc_TimeDomainDOF_InterwovenCoefficients(){
  
  std::vector<std::vector<double> > 
    cos_sin_coeffs(this->num_time_collocation_points , 
                   std::vector<double>(this->num_total_harmonics) );

  for(int time = 0 ; time < (this->num_time_collocation_points) ; time++){
    std::vector<double> coeffs(2 * (this->num_total_harmonics), 0.0);
    for(int harmonic = 0 ; harmonic < (this->num_total_harmonics) ; harmonic++){

      // small signal, ansatz: U= [U0] + sum_k [U_c^k cos(2*\pi*\omega_k*t) + U_s^k sin(2*\pi*\omega_k*t)]
      // where U_c^k cos(2*\pi*\omega_k*t) and U_s^k sin(2*\pi*\omega_k*t), for each k, solve the eqn separately and with the respective boundary condition at this harmonic
      // in other words, sine and cosine coefficients are coupled only for a given harmonic
      if(this->isSmallSignal){ 
        if(time == 2*harmonic) coeffs[time] = 1.0;
        if(time == 2*harmonic+1) coeffs[time] = 1.0;
      }

      // large signal
      if(!(this->isSmallSignal)){
        coeffs[2*harmonic    ] = std::cos(2 * M_PI * (*(this->remapped_harmonics))[harmonic] 
                                          * (this->time_collocation_points)[time] );
        coeffs[2*harmonic + 1] = std::sin(2 * M_PI * (*(this->remapped_harmonics))[harmonic] 
                                          * (this->time_collocation_points)[time] );
      }

      //std::cout << "At time " << std::to_string((this->time_collocation_points)[time]) << " and harmonic " << std::to_string(harmonic);
      //std::cout << ", cos(2 pi h time) equals " << std::to_string(coeffs[2*harmonic]);
      //std::cout << " and sin(2 pi h time) equals " << std::to_string(coeffs[2*harmonic+1]) << std::endl;

    }
    cos_sin_coeffs[time] = coeffs;
  }

  this->dof_interwoven_coeffs_at_time_coll_pt = cos_sin_coeffs;
}


// calculate the cosine time domain quadrature rule weights (for the remapped harmonics!)
void FreqDomParameters::calc_CosQuadratureWeightsForHarmonicNumber(){

  std::vector<std::vector<double> > 
    weights(this->num_total_harmonics,
            std::vector<double>(this->num_time_collocation_points));

  for(int harmonic = 0 ; harmonic < (this->num_total_harmonics) ; harmonic++){
    std::vector<double> coefficients( this->num_time_collocation_points , 0.0);
    double pi = 4.0*std::atan(1.0);
    for(int time = 0 ; time < (this->num_time_collocation_points) ; time++){

      // small signal
      // include only the COSINE summand involving this harmonic
      if(this->isSmallSignal)
        if(time == 2*harmonic) coefficients[time] = 1.0;

      // large signal
      if(!(this->isSmallSignal)){
        if(harmonic > 0)
          coefficients[time] = (2.0/(double)(this->num_time_collocation_points - 1)) 
                               * std::cos( 2 * pi * (*(this->remapped_harmonics))[harmonic] 
                               * (this->time_collocation_points)[time] );
        if(harmonic == 0)
          coefficients[time] = (1.0/(double)(this->num_time_collocation_points - 1)) 
                               * std::cos( 2 * pi * (*(this->remapped_harmonics))[harmonic] 
                               * (this->time_collocation_points)[time] );
      }
    }

    // large signal, account for the endpoits of the trapezoidal rule
    if(!(this->isSmallSignal)){
      coefficients[0] *= 0.5;
      coefficients[(this->num_time_collocation_points) - 1] *= 0.5; 
    }

    weights[harmonic] = coefficients;
  }

  this->cos_quadrature_weights_for_harmonic_number = weights;
}

// calculate the sine time domain quadrature rule weights (for the remapped harmonics!)
void FreqDomParameters::calc_SinQuadratureWeightsForHarmonicNumber(){

  std::vector<std::vector<double> > weights(this->num_total_harmonics, 
                                            std::vector<double>(this->num_time_collocation_points));

  for(int harmonic = 0 ; harmonic < (this->num_total_harmonics) ; harmonic++){
    std::vector<double> coefficients( this->num_time_collocation_points , 0);
    double pi = 4.0*std::atan(1.0);
    for(int time = 0 ; time < (this->num_time_collocation_points) ; time++){

      // small signal (include only the SINE summand of this harmonic)
      if(this->isSmallSignal)
        if(time == 2*harmonic+1) coefficients[time] = 1.0;

      // large signal
      if(!(this->isSmallSignal)){
        if(harmonic > 0)
          coefficients[time] = (2.0/(double)(this->num_time_collocation_points - 1)) 
                               * std::sin( 2 * pi * (*(this->remapped_harmonics))[harmonic] 
                               * (this->time_collocation_points)[time] );
        if(harmonic == 0)
          coefficients[time] = (1.0/(double)(this->num_time_collocation_points - 1)) 
                               * std::sin( 2 * pi * (*(this->remapped_harmonics))[harmonic] 
                               * (this->time_collocation_points)[time] );
      }
    }

    // large signal, account for the endpoits of the trapezoidal rule
    if(!(this->isSmallSignal)){
      coefficients[0] *= 0.5;
      coefficients[(this->num_time_collocation_points) - 1] *= 0.5; 
    }

    weights[harmonic] = coefficients;
  }

  this->sin_quadrature_weights_for_harmonic_number = weights;
}

///////////////////////////////////////////////////////////////////////////////////////
// get functions
bool FreqDomParameters::queryIsSmallSignal(){
  return this->isSmallSignal; }

Teuchos::Array<double> FreqDomParameters::getHarmonicsVector(){
  return this->harmonics_vector; }

int FreqDomParameters::getTruncationOrder(){
  return this->truncation_order; }

std::string FreqDomParameters::getTruncationScheme(){
  return this->truncation_scheme; }

Teuchos::Array<double> FreqDomParameters::getFundamentalHarmonics(){
  return this->fundamental_harmonics; }

int FreqDomParameters::getNumFundamentalHarmonics(){
  return this->num_fundamental_harmonics; }

std::vector<std::vector<int> > FreqDomParameters::getPossibleMultiIndices(){ 
  return this->possible_multi_indices; }

Teuchos::RCP<std::vector<double> > 
  FreqDomParameters::getTruncatedHarmonicBasis(){ 
    return this->truncated_harmonic_basis; }

Teuchos::RCP<std::vector<double> > 
  FreqDomParameters::getRemappedFundamentalHarmonics(){ 
    return this->remapped_fundamental_harmonics; }

Teuchos::RCP<std::vector<double> > FreqDomParameters::getUnRemappedHarmonics(){
  return this->unremapped_harmonics; }

Teuchos::RCP<std::vector<double> > FreqDomParameters::getRemappedHarmonics(){
  return this->remapped_harmonics; }

int FreqDomParameters::getNumTotalHarmonics(){
  return this->num_total_harmonics; }

int FreqDomParameters::getNumTimeCollocationPoints(){
  return this->num_time_collocation_points; }

std::vector<double> FreqDomParameters::getTimeCollocationPoints(){ 
  return this->time_collocation_points; }

std::vector<std::vector<double> > 
  FreqDomParameters::getDofInterwovenCoeffsAtTimeCollPt(){ 
    return this->dof_interwoven_coeffs_at_time_coll_pt; }

std::vector<std::vector<double> > 
  FreqDomParameters::getCosQuadratureWeightsForHarmonicNumber(){
    return this->cos_quadrature_weights_for_harmonic_number; }
std::vector<std::vector<double> > 
  FreqDomParameters::getSinQuadratureWeightsForHarmonicNumber(){
    return this->sin_quadrature_weights_for_harmonic_number; }
