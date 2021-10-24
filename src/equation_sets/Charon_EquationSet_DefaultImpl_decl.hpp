//  idea: for now, instead of adding the new HB-useful constructor to Panzer_EquationSet_DefaultImpl*, create an independent fancy EquationSet in Charon which has the old style constructor along with the HB-style constructor. Then, have the Charon Equation Sets inherit from this fancy one.

//  idea: panzer equationset composite structure:
//  - take as input a vector of equation sets, composite parameters (from a separate CompositeParameters.cpp? with default Fourier implementation...)
//  - virtual functions: 
//      setupDOFs
//      buildAndRegisterCompositeEquationSetEvaluators


// the EquationSet_Composite should take care of all of the necessary DOF stuff, calling EquationSet_CompositeParameters as necessary
//   to check that they are created appropriately from the EquationSet_Constituents
// let the EquationSet_Constituents handle

#ifndef CHARON_EQUATIONSET_DEFAULTIMPL_DECL_HPP
#define CHARON_EQUATIONSET_DEFAULTIMPL_DECL_HPP

#include "Charon_Names.hpp"

#include "PanzerDiscFE_config.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"

// used to create default parameters
#include "Panzer_GlobalData.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace PHX {
  template<typename Traits> class FieldManager;
}

namespace charon {

  template <typename EvalT>
  class EquationSet_DefaultImpl : public panzer::EquationSet_DefaultImpl<EvalT> {

  public:

    // With DOF setup constructor, inherited
    using panzer::EquationSet_DefaultImpl<EvalT>::EquationSet_DefaultImpl;

    // Sans DOF setup constructor, defined in impl
    EquationSet_DefaultImpl(Teuchos::RCP<panzer::IntegrationRule> input_ir,
                            Teuchos::RCP<panzer::BasisIRLayout> input_basis,
			    Teuchos::ParameterList timedom_residual_pl);

    virtual ~EquationSet_DefaultImpl() {}

  protected:

    struct DOFDescriptor {
      DOFDescriptor() 
        : dofName("")
        , residualName(std::make_pair(false,""))
        , grad(std::make_pair(false,""))
        , curl(std::make_pair(false,""))
        , div(std::make_pair(false,""))
        , timeDerivative(std::make_pair(false,"")) {}

      std::string dofName;
      std::string basisType;
      int basisOrder;
      Teuchos::RCP<panzer::PureBasis> basis;
      int integrationOrder;
      Teuchos::RCP<panzer::IntegrationRule> intRule;
      std::pair<bool,std::string> residualName;
      std::string scatterName;
      std::pair<bool,std::string> grad;
      std::pair<bool,std::string> curl;
      std::pair<bool,std::string> div;
      std::pair<bool,std::string> timeDerivative;

      void print(std::ostream & os) const {
        os << "DOF Desc = \"" << dofName << "\": "
           << "Basis = (" << basisType << ", \"" << basisOrder << "\"), "
           << "Res = (" << residualName.first << ", \"" << residualName.second << "\"), "
           << "Grad = (" << grad.first << ", \"" << grad.second << "\"), "
           << "Curl = (" << curl.first << ", \"" << curl.second << "\"), "
           << "Div = (" << div.first << ", \"" << div.second << "\"), "
           << "Time = (" << timeDerivative.first << ", \"" << timeDerivative.second << "\")";
      }
    };

  private:

    //! For convenience, declare the DOFDescriptor iterator
    typedef typename std::map<std::string,DOFDescriptor>::const_iterator DescriptorIterator;

    //! For convenience, declare a basis iterator
    typedef typename std::map<std::string,std::pair<Teuchos::RCP<panzer::PureBasis>,Teuchos::RCP<std::vector<std::string> > > >::const_iterator BasisIterator;


  public:

    // frequency domain boolean
    bool isFreqDom = false;

    // for frequency domain simulations
    Teuchos::RCP<panzer::IntegrationRule> base_ir;
    Teuchos::RCP<panzer::BasisIRLayout> base_basis;    
    Teuchos::ParameterList m_timedom_residual_pl;

    Teuchos::RCP<charon::Names> base_names;

    // options
    std::string base_solveElectron;
    std::string base_solveHole;

    std::string base_supg_stab;
    std::string base_tau_e_type;
    std::string base_tau_h_type;
    std::string base_ls_type;

    bool base_haveSource;
    bool base_add_source_stab;
    bool base_addTrapCharge;
    bool base_addFixCharge;
    std::string base_drForce;

  };
}

#endif //CHARON_EQUATIONSET_DEFAULTIMPL_DECL_HPP
