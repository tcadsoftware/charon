
#ifndef CHARON_PHYSICALCONSTANTS_HPP
#define CHARON_PHYSICALCONSTANTS_HPP

namespace charon {

/**
 * @brief Stores common physical constants
 *
 */
class PhysicalConstants {

public:
  static PhysicalConstants const& Instance();

private:
  PhysicalConstants();

public:
  //! Boltzmann constant in [eV/K]
  const double kb;

  //! Electron elemental charge in [C]
  const double q;

  //! Vacuum permittivity in [C/(V.cm)]
  const double eps0;

  //! Electron rest mass (kg)
  const double m0;

  //! Planck's constant (J-sec)
  const double h;

  //! Reduced Planck's constant (J-sec)
  const double hbar;

  //! pi constant (unitless)
  const double pi;

};

} // END "namespace charon"

#endif
