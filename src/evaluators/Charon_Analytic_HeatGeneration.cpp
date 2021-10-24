
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_Analytic_HeatGeneration_decl.hpp"
#include "Charon_Analytic_HeatGeneration_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Analytic_HeatGeneration)

#endif
