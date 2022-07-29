INCLUDE(TribitsAddAdvancedTest)

##########################################################
# Test the DDIon equation set (only FEM-SUPG is supported)
##########################################################

SET(testName charon_mp_ddion)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    memristor.jou
    memristor.exo.4.0 memristor.exo.4.1 memristor.exo.4.2 memristor.exo.4.3
    memristor.nlp.inp
    memristor.dd.equ.inp
    memristor.ddion.trans.3V.inp
    memristor.ddion.iv.loca.inp
    memristor.ddion.sine.inp
    memristor.ddion.iv.loca.gold.exo
    memristor.ddion.sine.gold.exo
    memristor.ddion.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.nlp.inp

  # Nonlinear Poisson simulation
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.dd.equ.inp

  # Equilibrium DD simulation
  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.dd.equ.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddion.trans.3V.inp

  # Transient DDIon simulation
  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.ddion.trans.3V.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for next test
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "memristor.ddion.trans.3V.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddion.iv.loca.inp

  # Steady-state DDIon simulation (ION_DENSITY is frozen)
  # This simulation converges only in a serial run
#  TEST_9 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=memristor.ddion.iv.loca.inp.xml --current
#    NUM_MPI_PROCS 1
#    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Compare simulation results
#  TEST_10 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
#    NOEXESUFFIX
#    ARGS -f
#         memristor.ddion.exodiff
#         memristor.ddion.iv.loca.exo
#         memristor.ddion.iv.loca.gold.exo
#    NUM_MPI_PROCS 1

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddion.sine.inp

  # Transient DDIon simulation with a sinusoidal voltage source
  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.ddion.sine.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "memristor.ddion.sine.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
#  TEST_13 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
#    NOEXESUFFIX
#    ARGS -f
#         memristor.ddion.exodiff
#         memristor.ddion.sine.exo
#         memristor.ddion.sine.gold.exo
#    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ddion PROPERTY LABELS ddion nightly debugexclude)
ENDIF()
