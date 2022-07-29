#include "Charon_config.hpp"

#if 1

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_NormCalculation.hpp"
#include "Charon_NormCalculation_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Norm_L2Error)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Norm_L2)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Norm_H1Error)
PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::Norm_H1)

#endif

