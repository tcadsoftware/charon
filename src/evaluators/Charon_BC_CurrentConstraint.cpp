
#include "Charon_config.hpp"

#ifdef    HAVE_CHARON_EXPLICIT_INSTANTIATION

// Charon
#include "Charon_BC_CurrentConstraint_decl.hpp"
#include "Charon_BC_CurrentConstraint_impl.hpp"

// Panzer
#include "Panzer_ExplicitTemplateInstantiation.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::BC_CurrentConstraint)

#endif // HAVE_CHARON_EXPLICIT_INSTANTIATION

// end of Charon_BC_CurrentConstraint.cpp
