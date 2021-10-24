
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_ClosureModel_Factory_decl.hpp"
#include "Charon_ClosureModel_Factory_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::ClosureModelFactory)

#endif
