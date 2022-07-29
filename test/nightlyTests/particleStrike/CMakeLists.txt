INCLUDE(TribitsAddAdvancedTest)

####################################################
# Diode particle Strike test
SET(testName charon_mp_particleStrike)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    diode-particleStrike.sg.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    diode-particleStrike.gold.exo
    diode-particleStrike.exodiff

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i diode-particleStrike.sg.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=diode-particleStrike.sg.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "diode-particleStrike.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         diode-particleStrike.exodiff
         diode-particleStrike.exo
         diode-particleStrike.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_particleStrike PROPERTY LABELS
    particleStrike nightly debugexclude)
ENDIF()


####################################################
# 3D Diode Particle Strike test
SET(testName charon_mp_3DparticleStrike)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    diode3D-particleStrike.sg.inp
    diode3D.exo.4.0
    diode3D.exo.4.1
    diode3D.exo.4.2
    diode3D.exo.4.3
    diode3D-particleStrike.gold.exo
    diode3D-particleStrike.exodiff

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i diode3D-particleStrike.sg.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=diode3D-particleStrike.sg.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "diode3D-particleStrike.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         diode3D-particleStrike.exodiff
         diode3D-particleStrike.exo
         diode3D-particleStrike.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_3DparticleStrike PROPERTY LABELS
    particleStrike nightly debugexclude)
ENDIF()
