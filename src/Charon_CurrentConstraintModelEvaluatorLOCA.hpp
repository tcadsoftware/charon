#ifndef __Charon_CurrentConstraintModelEvaluatorLOCA_hpp__
#define __Charon_CurrentConstraintModelEvaluatorLOCA_hpp__

#include "Charon_CurrentConstraintList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

namespace charon
{
  /**
   *  \brief LOCA Current Constraint `ModelEvaluator`.
   *
   *  This ModelEvaluator intercepts the responses in the physics ME
   *  and replaces them with constraint equations required by LOCA.
   *
   *  This `ModelEvaluator` is used to apply a current constraint to a
   *  boundary.  Two different types of constraints are available:
   *  - Constant Current:\n
   *    In this case we're hooking up a current source to a boundary.  We
   *    therefore need to specify that current value.
   *  - Resistor Contact:\n
   *    In this case we're attaching a resistor to a boundary, and then hooking
   *    up a voltage source to that resistor.  We therefore need to specify the
   *    resistor value and the applied voltage.
   */
  template<typename Scalar>
  class CurrentConstraintModelEvaluatorLOCA
    :
    public Thyra::ModelEvaluatorDelegatorBase<Scalar>
  {
    public:
    
    /**
     *  \brief Physics Constructor.
     *
     *  This constructor saves the input as member data and then calls
     *  `initialize()`.
     *
     *  \param[in] physics          The physics `ModelEvaluator`.
     *  \param[in] comm             The MPI communicator.
     *  \param[in] constraints      The list of all the current constraints.
     *  \param[in] dimension        The spatial dimension of the simulation.
     *  \param[in] print_debug      If set to true, debug info will be dumped to screen.
     */
    CurrentConstraintModelEvaluatorLOCA(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>>& physics,
                                        MPI_Comm comm,
                                        const charon::CurrentConstraintList& constraints,
                                        const int& dimension,
					const Teuchos::RCP<Teuchos::ParameterList>& parameters,
                                        const bool print_debug = false);

    CurrentConstraintModelEvaluatorLOCA() = delete;    
    CurrentConstraintModelEvaluatorLOCA(const CurrentConstraintModelEvaluatorLOCA& src) = delete;
    
    /// Returns the current constraints
    const charon::CurrentConstraintList& constraints() const
    {return constraints_;}

    /// Copies value from Tpetra::MultiVector to Thyra::Simd.
    void assignValueTpetraToSpmd(const Teuchos::RCP<const ::Thyra::VectorBase<Scalar>>& source_tpetra_value,
                                 const Teuchos::RCP<::Thyra::VectorBase<Scalar>>& target_simd_value) const;

    /// Copies value from Thyra::Simd into Tpetra::MultiVector.
    void assignValueSpmdToTpetra(const Teuchos::RCP<const ::Thyra::MultiVectorBase<Scalar>>& source_simd_value,
                                 const Teuchos::RCP<::Thyra::MultiVectorBase<Scalar>>& target_tpetra_value) const;

    /// Returns the underlying physics model evalautor.
    Teuchos::RCP<::Thyra::ModelEvaluator<Scalar>> getInternalPhysicsME() const;

    // Thr following overrides are needed to swap Thyra::SpmdVectors
    // for locally replicated p, g and DgDp (allocated in
    // panzer::ModelEvaluator) with a corresponding
    // Tpetra::MultiVector versions. LOCA constraint enforcement
    // requires Tpetra::MV data structures.
    std::string description() const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_p_space(int l) const;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> get_g_space(int j) const;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
    Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
    Teuchos::RCP<Thyra::LinearOpBase<Scalar>> create_DgDp_op( int j, int l ) const;

  private:    

    void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar>&  inArgs,
                       const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
    
    /**
     *  \brief Get the voltages.
     *
     *  This grabs the voltages out of the constrained solution vector.
     *
     *  \param[in,out] voltages The voltages.
     *  \param[in]     xFull    The constrained solution vector.
     */
    // void getVoltages(std::vector<Scalar>& voltages,
    //                  Teuchos::RCP<const Thyra::ProductVectorBase<Scalar>> xFull) const;

    /// Underlying physics `ModelEvaluator`.
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> physics_;

    /// The MPI communicator to use with this `ModelEvaluator`.
    Teuchos::RCP<const Teuchos::Comm<int>> comm_;

    /// The list of all the current constraints.
    const charon::CurrentConstraintList constraints_;

    /// The spatial dimension of the simulation (2D or 3D);
    const int dimension_;

    Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> pSpace_;
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> gSpace_;

    Teuchos::RCP<Teuchos::ParameterList> parameters_;
    const bool print_debug_;
  };

} // namespace charon

#endif
