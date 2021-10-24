
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_BCStrategy_Dirichlet_Constant_decl.hpp"
#include "Charon_BCStrategy_Dirichlet_Constant_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::BCStrategy_Dirichlet_Constant)

#endif
