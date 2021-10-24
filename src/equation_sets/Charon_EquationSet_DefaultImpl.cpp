#include "PanzerDiscFE_config.hpp"

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"
#include "Panzer_Traits.hpp"

#include "Charon_EquationSet_DefaultImpl_decl.hpp"
#include "Charon_EquationSet_DefaultImpl_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::EquationSet_DefaultImpl)

#endif
