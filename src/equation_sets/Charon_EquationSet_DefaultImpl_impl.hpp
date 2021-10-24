#ifndef CHARON_EQUATIONSET_DEFAULTIMPL_IMPL_HPP
#define CHARON_EQUATIONSET_DEFAULTIMPL_IMPL_HPP

// Sans DOF setup constructor
template <typename EvalT>
charon::EquationSet_DefaultImpl<EvalT>::
EquationSet_DefaultImpl(Teuchos::RCP<panzer::IntegrationRule> input_ir,
                        Teuchos::RCP<panzer::BasisIRLayout> input_basis,
                        Teuchos::ParameterList timedom_residual_pl)
:
panzer::EquationSet_DefaultImpl<EvalT>::EquationSet_DefaultImpl(
	timedom_residual_pl.get<Teuchos::RCP<Teuchos::ParameterList>>("params"),
        timedom_residual_pl.get<int>("default_integration_order"),
	panzer::CellData(),
        timedom_residual_pl.get<Teuchos::RCP<panzer::GlobalData>>("global_data"),
        timedom_residual_pl.get<bool>("build_transient_support"))
{
    isFreqDom = true;

    base_ir = input_ir;
    base_basis = input_basis;
    m_timedom_residual_pl = timedom_residual_pl;

    base_names = timedom_residual_pl.get<Teuchos::RCP<charon::Names>>("Names");

    // Options
    base_solveElectron = timedom_residual_pl.get<std::string>("solveElectron");
    base_solveHole = timedom_residual_pl.get<std::string>("solveHole");

    if(timedom_residual_pl.isParameter("supg_stab"))
      base_supg_stab = timedom_residual_pl.get<std::string>("supg_stab");
    if(timedom_residual_pl.isParameter("tau_e_type"))
      base_tau_e_type = timedom_residual_pl.get<std::string>("tau_e_type");
    if(timedom_residual_pl.isParameter("tau_h_type"))
      base_tau_h_type = timedom_residual_pl.get<std::string>("tau_h_type");
    if(timedom_residual_pl.isParameter("ls_type"))
      base_ls_type = timedom_residual_pl.get<std::string>("ls_type");

    if(timedom_residual_pl.isParameter("haveSource"))
      base_haveSource = timedom_residual_pl.get<bool>("haveSource");
    if(timedom_residual_pl.isParameter("add_source_stab"))
      base_add_source_stab = timedom_residual_pl.get<bool>("add_source_stab");
    if(timedom_residual_pl.isParameter("addTrapCharge"))
      base_addTrapCharge = timedom_residual_pl.get<bool>("addTrapCharge");
    if(timedom_residual_pl.isParameter("addFixCharge"))
      base_addFixCharge = timedom_residual_pl.get<bool>("addFixCharge");
      
    // note that SGCVFEM Drift Diffusion equation set uses this string
    // whereas Drift Diffusion will not, even though this is set to a default value here
    base_drForce = timedom_residual_pl.get<std::string>("drForce","EffectiveField");
}

#endif // CHARON_EQUATIONSET_DEFAULTIMPL_IMPL_HPP
