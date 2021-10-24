

#include "Charon_Physical_Constants.hpp"

//========================================================================
charon::PhysicalConstants const& charon::PhysicalConstants::Instance()
{
  static PhysicalConstants instance;
  return instance;
}

//========================================================================
charon::PhysicalConstants::PhysicalConstants() :
  kb(8.617343e-05),
  q(1.602176487e-19),
  eps0(8.854187817e-12*0.01),

  m0(9.10938215e-31),
  h(6.62606957e-34),
  hbar(1.054571628e-34),
  pi(3.141592654)

// constants needed to run the NLP MMS simulations for code verification
//  kb(8.618333e-05),
//  q(1.602e-19),
//  eps0(8.853559e-12*0.01),

// constants found in Bill Wampler's he3 code
//  kb(1./11604.0),
//  q(1.60218e-19),
//  eps0(8.854187817e-12*0.01),
//  m0(0.910953e-30),
//  h(6.62617e-34),
//  hbar(h/(2*pi)),
//  pi(3.141592653589793)
{
}
