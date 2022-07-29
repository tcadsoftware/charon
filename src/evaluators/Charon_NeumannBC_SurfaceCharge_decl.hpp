
#ifndef CHARON_NEUMANNBC_SURFACECHARGE_DECL_HPP
#define CHARON_NEUMANNBC_SURFACECHARGE_DECL_HPP

#include "Panzer_Dimension.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Panzer_ParameterLibrary.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Charon_Names.hpp"

using panzer::Cell;
using panzer::Point;


namespace charon {


/**
 *  \brief Neumann boundary condition evaluator for surfaces and
 *         interfaces.
 *
 *  Calculates Neumann boundary condition for charge specified by the user 
 *  in #/cm^2.
 *  Or the user can specify piezo and spontaneous polarization for 
 *  wurtzite materials.
 *  Using the equations given in J. Appl. Phys. Vol. 87, No. 1 the charge 
 *  density due to piezo and spontaneous polarization at an interface can be
 *  calculated.
 */

template<typename EvalT, typename Traits>
class NeumannBC_SurfaceCharge
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    NeumannBC_SurfaceCharge(const Teuchos::ParameterList& p);

    void evaluateFields(typename Traits::EvalData d);

    void postRegistrationSetup(typename Traits::SetupData d,
                               PHX::FieldManager<Traits>& fm);

  private:
    using ScalarT = typename EvalT::ScalarT;

    /**
     * \brief linearly interpolated piezoelectric constants
     *  
     */
    std::map<std::string,double> interp;

    //std::string paramName;
    
    /**
     *  \brief Produces a map of material values for use by piezo and
     *         polarization functions
     *
     *  \param[in]  arg1  Wurtzite material properties.
     *
     *  \returns A std::map of wurtzite material properties
     */
    auto getPolarvals(const std::string& mat) -> std::map<std::string, double>;

    /**
     *  \brief Calculates the piezo electric charge between given materials
     *  in C/m^2
     *
     *  \param[in] arg1 getPolarvals for top material in stack.
     *
     *  \param[in] arg2 getPolarvals for bottom material in stack.
     *
     *  \param[in] arg3 metal fraction i.e \f$ Al_xGa_{1-x}N \f$, optional
     *                  defaults to 0.3.
     *
     *  \returns charge in C/m^2
     *
     */
    auto piezo(const std::map<std::string,double>& top, 
          const std::map<std::string,double>& bottom,
          double xcomp) -> double;

    /**
     *  \brief Calculates the spontaneous polarization charge between given 
     *  materials in C/m^2
     *
     *  \note this function calls the piezo function which populates the
     *        interp map
     *
     *  \param[in] arg1 getPolarvals for top material in stack.
     *
     *  \param[in] arg2 getPolarvals for bottom material in stack.
     *
     *  \param[in] arg3 metal fraction i.e \f$ Al_xGa_{1-x}N \f$, optional
     *                  defaults to 0.3.
     *
     *  \returns charge in C/m^2
     *
     */
    auto polarization(const std::map<std::string,double>& top, 
          const std::map<std::string,double>& bottom,
          double xcomp) -> double;

    // output fields
    PHX::MDField<ScalarT,Cell,Point> surface_charge;   // scaled, no unit
    PHX::MDField<ScalarT,Cell,Point> surface_recomb;   // scaled, no unit

    // input fields
    PHX::MDField<const ScalarT,Cell,Point> edensity;      // scaled, no unit
    PHX::MDField<const ScalarT,Cell,Point> hdensity;      // scaled, no unit
    PHX::MDField<const ScalarT,Cell,Point> intrin_conc;   // scaled, no unit
    PHX::MDField<const ScalarT,Cell,Point> latt_temp;     // scaled, no unit

    double C0, X0, R0, T0;  // scaling parameters
    double varyingCharge;     // in unit of cm^{-2}
    Teuchos::RCP<panzer::ScalarParameterEntry<EvalT> > fixedCharge;     // in unit of cm^{-2}

    int num_ips, num_nodes; 
    std::string basis_name; 
    std::size_t basis_index;
 
    bool bFixCharge, bVaryingCharge, bSurfTrap, bSurfRecomb, bPolar;

    std::string fluxSurfCharge, fluxSurfRecomb; 
 
    std::vector<double> trapEnergy, trapDensity, eTrapVelocity, hTrapVelocity, trapEnWidth; 
    std::vector<std::string> trapType, trapEnDistr;
    std::vector<int> trapNL;
    void discretizeContDistribution(std::vector<ScalarT>* enLevels, 
				 std::vector<ScalarT>* norm_densities,
				 const std::string& enDistr, double Et, 
				 double enSigma, int NL);
    
    double eSurfRecombVel, hSurfRecombVel, surfRecombEnergy; 
  
    Teuchos::RCP<Teuchos::ParameterList> getValidParameters() const;

    void initSurfTrapParams(const Teuchos::ParameterList& trapPList);
 
}; // end of class NeumannBC_SurfaceCharge


}

#endif
