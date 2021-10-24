
#include "Charon_config.hpp"

#ifdef    HAVE_CHARON_EXPLICIT_INSTANTIATION

// Charon
#include "Charon_BCStrategy_Dirichlet_CurrentConstraint_decl.hpp"
#include "Charon_BCStrategy_Dirichlet_CurrentConstraint_impl.hpp"

// Panzer
#include "Panzer_ExplicitTemplateInstantiation.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(
  charon::BCStrategy_Dirichlet_CurrentConstraint)

#endif // HAVE_CHARON_EXPLICIT_INSTANTIATION

// end of Charon_BCStrategy_Dirichlet_CurrentConstraint.cpp
