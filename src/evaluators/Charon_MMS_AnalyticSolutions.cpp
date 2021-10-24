
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_MMS_AnalyticSolutions.hpp"
#include "Charon_MMS_AnalyticSolutions_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::MMS_NLP_GLH_1_AnalyticSolution)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::MMS_DD_RDH_1_AnalyticSolution)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::MMS_DD_RDH_2_AnalyticSolution)

#endif
