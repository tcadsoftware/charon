
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_BC_OhmicContact_decl.hpp"
#include "Charon_BC_OhmicContact_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::BC_OhmicContact)

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::OhmicContact)

#endif
