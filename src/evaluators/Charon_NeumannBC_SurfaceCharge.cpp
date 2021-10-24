
#include "Charon_config.hpp"

#ifdef HAVE_CHARON_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Charon_NeumannBC_SurfaceCharge_decl.hpp"
#include "Charon_NeumannBC_SurfaceCharge_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(charon::NeumannBC_SurfaceCharge)

#endif
