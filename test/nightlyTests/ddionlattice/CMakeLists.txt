INCLUDE(TribitsAddAdvancedTest)

################################################################################
# Test the FEM-SUPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.supg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    tiox.ddionlatt.supg.yaml
    tiox.exo.4.0
    tiox.exo.4.1
    tiox.exo.4.2
    tiox.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_supg.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=tiox.ddionlatt.supg.yaml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE supg.log

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.supg  PROPERTY LABELS tiox.ddionlatt.supg  nightly)
ENDIF()

################################################################################
# Test the FEM-EFFPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.effpg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    tiox.ddionlatt.effpg.yaml
    tiox.exo.4.0
    tiox.exo.4.1
    tiox.exo.4.2
    tiox.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_effpg.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=tiox.ddionlatt.effpg.yaml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE effpg.log

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.effpg  PROPERTY LABELS
    tiox.ddionlatt.effpg  nightly)
ENDIF()

################################################################################
# Test the FEM-SymEFFPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.symeffpg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    tiox.ddionlatt.symeffpg.yaml
    tiox.exo.4.0
    tiox.exo.4.1
    tiox.exo.4.2
    tiox.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_symeffpg_rpi.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=tiox.ddionlatt.symeffpg.yaml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE symeffpg.log

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.symeffpg  PROPERTY LABELS
    tiox.ddionlatt.symeffpg  nightly)
ENDIF()

# NUM_MPI_PROCS has to be <= MPI_EXEC_MAX_NUMPROCS set in charon-cmake.sh script.
# See the email below from Roger on 08/19/2015.
#
# Hi Suzey,
#
# When you configure tribits, there is a deault for the maximum number of mpi processes per test added.
# This is to not overload the machine.  You just need to reconfigure with a default larger than 4.
# See this section in the tribits user guide:
# https://tribits.org/doc/TribitsBuildReference.html#configuring-with-mpi-support
#
# From the link above:
#
# -D MPI_EXEC_MAX_NUMPROCS=4
#
# (The maximum number of processes to allow when setting up and running MPI test and example executables.
# The default is set to '4' but should be set to the largest number that can be tolerated for the given machine.
# Tests with more processes than this are excluded from the test suite at configure time.)
#
# Roger
