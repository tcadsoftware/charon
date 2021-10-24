
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_SUPG_Tau_Linear_decl.hpp"
#include "Charon_SUPG_Tau_Linear_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::SUPG_Tau_Linear)

#endif
