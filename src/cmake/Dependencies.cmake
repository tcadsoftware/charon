SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Kokkos Panzer PanzerDiscFE PanzerAdaptersSTK STKMesh STKUtil STKTopology SEACASNemesis SEACASExodus Stratimikos Piro NOX Rythmos Teko Belos AztecOO)

SET(TEST_REQUIRED_DEP_PACKAGES SEACAS SEACASEpu SEACASExodiff SEACASNemslice SEACASNemspread)

IF (ENABLE_XYCE_CLUSTER)

  # The coupled Charon/Xyce cluster build requires additional Boost
  # libraries and the xyce libraries.
  SET(LIB_REQUIRED_DEP_TPLS BoostCharonLib XyceCharonLib)

  # Add TrilinosCouplings for Xyce
  SET(LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES} TrilinosCouplings)
ELSE()

  # Set to empty for non cluster builds
  SET(LIB_REQUIRED_DEP_TPLS)
ENDIF()

# Set some variables for specific platforms. This can be used, for
# example, to run certain tests on only that platform.
IF ($ENV{SNLCLUSTER} MATCHES "skybridge")
  SET(CHARON_CLUSTER_SKYBRIDGE TRUE)
ELSEIF ($ENV{SNLCLUSTER} MATCHES "linux_rh7")
  SET(CHARON_CLUSTER_RHEL7 TRUE)
ENDIF()

IF (ENABLE_HEAVY_TESTS)
  SET(MPI_EXEC_MAX_NUMPROCS_DEFAULT 16)
  MESSAGE(STATUS "WARNING!!! Heavy tests only enabled!")
ENDIF()

SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
