
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_RecombRate_Empirical_Defect_decl.hpp"
#include "Charon_RecombRate_Empirical_Defect_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(charon::RecombRate_Empirical_Defect,panzer::BASIS)
PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(charon::RecombRate_Empirical_Defect,panzer::Point)

#endif
