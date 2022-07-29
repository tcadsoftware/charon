#include "Charon_config.hpp"

#if 1

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_ManufacturedSolution.hpp"
#include "Charon_ManufacturedSolution_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::AnalyticSolution)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::DD_RDH_1_AnalyticSolution)

#endif
