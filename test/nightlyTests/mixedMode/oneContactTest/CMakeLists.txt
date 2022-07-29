
INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run manufactured and mixed mode problems, then compare
SET(testName charon_mp_mixedMode.oneContactTest)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistorC_resistorX.cir
    resistorC_resistorX.dd.inp
    resistor.dd.equ.exo.4.0
    resistor.dd.equ.exo.4.1
    resistor.dd.equ.exo.4.2
    resistor.dd.equ.exo.4.3
    resistorC_resistorX.dd.gold.exo
    pndiodeC_resistorX.cir
    pndiodeC_resistorX.dd.inp
    pndiode.dd.exo.4.0
    pndiode.dd.exo.4.1
    pndiode.dd.exo.4.2
    pndiode.dd.exo.4.3
    pndiodeC_resistorX.dd.gold.exo
    tol.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistorC_resistorX.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistorC_resistorX.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistorC_resistorX.dd.inp.log

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistorC_resistorX.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         resistorC_resistorX.dd.exo
         resistorC_resistorX.dd.gold.exo
    NUM_MPI_PROCS 1

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiodeC_resistorX.dd.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiodeC_resistorX.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiodeC_resistorX.dd.inp.log

  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiodeC_resistorX.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         pndiodeC_resistorX.dd.exo
         pndiodeC_resistorX.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mixedMode.oneContactTest PROPERTY LABELS mixedMode.oneContactTest
    mixedmode)
ENDIF()
