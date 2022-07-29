
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_RecombRate_DynamicTraps_decl.hpp"
#include "Charon_RecombRate_DynamicTraps_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::RecombRate_DynamicTraps)
PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::DynamicTraps)
PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::Trap)

#endif
