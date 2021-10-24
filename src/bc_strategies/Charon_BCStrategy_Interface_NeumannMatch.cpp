
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_BCStrategy_Interface_NeumannMatch_decl.hpp"
#include "Charon_BCStrategy_Interface_NeumannMatch_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(charon::BCStrategy_Interface_NeumannMatch)

#endif
