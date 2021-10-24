
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_BCStrategy_Dirichlet_Sinusoid_decl.hpp"
#include "Charon_BCStrategy_Dirichlet_Sinusoid_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::BCStrategy_Dirichlet_Sinusoid)

#endif
