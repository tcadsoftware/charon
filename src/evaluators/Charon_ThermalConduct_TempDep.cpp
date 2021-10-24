
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_ThermalConduct_TempDep_decl.hpp"
#include "Charon_ThermalConduct_TempDep_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::ThermalConduct_TempDep)

#endif
