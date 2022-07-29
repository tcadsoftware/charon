SET(LIB_REQUIRED_DEP_PACKAGES Teuchos Kokkos Panzer PanzerDiscFE PanzerAdaptersSTK STKMesh STKUtil STKTopology SEACASNemesis SEACASExodus SEACASAprepro Stratimikos Piro NOX Rythmos Teko Belos AztecOO)

SET(TEST_REQUIRED_DEP_PACKAGES SEACAS SEACASEpu SEACASExodiff SEACASNemslice SEACASNemspread)

IF (ENABLE_XYCE_CLUSTER OR ENABLE_MIXED_MODE)

  # The coupled Charon/Xyce cluster build requires additional Boost
  # libraries and the xyce libraries.
  SET(LIB_REQUIRED_DEP_TPLS BoostCharonLib XyceCharonLib)

  # Add TrilinosCouplings for Xyce
  SET(LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES} TrilinosCouplings)
ELSE()

  # Set to empty for non cluster builds
  SET(LIB_REQUIRED_DEP_TPLS)
ENDIF()

# For excluding tests from skybridge. The weird logic is due to the
# fact that tests can only be excluded via Tribits EXCLUDE_IF_NOT_TRUE
# directive in the test. There is no EXCLUDE_IF_TRUE directive.
IF ($ENV{SNLCLUSTER} MATCHES "skybridge")
  SET(CHARON_CLUSTER_NOTSKYBRIDGE FALSE)
ELSE()
  SET(CHARON_CLUSTER_NOTSKYBRIDGE TRUE)
ENDIF()

SET(MPI_EXEC_MAX_NUMPROCS_DEFAULT 16)

SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
