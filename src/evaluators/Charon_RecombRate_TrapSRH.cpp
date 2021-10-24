
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_RecombRate_TrapSRH_decl.hpp"
#include "Charon_RecombRate_TrapSRH_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::RecombRate_TrapSRH)
// PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::FermiDiracIntrinsicDensity)

#endif
