
INCLUDE(TribitsListHelpers)

SET( src_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Charon                 .       SS
  )

#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Charon Windows)
