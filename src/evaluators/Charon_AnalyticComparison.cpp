
#include "Charon_config.hpp"

#if 1

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_AnalyticComparison.hpp"
#include "Charon_AnalyticComparison_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::AnalyticComparison)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::AnalyticComparison_L2Error)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::AnalyticComparison_RelError)
//PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::AnalyticSolution)
//PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::ManufacturedSolution)

#endif
