
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_Integrator_GradBasisDotGradDOF_decl.hpp"
#include "Charon_Integrator_GradBasisDotGradDOF_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Integrator_GradBasisDotGradDOF)

#endif
