
INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run manufactured and mixed mode problems, then compare
SET(testName charon_mp_mixedMode.loca.twoContact)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistorC_linearX.cir
    #resistorC_linearX.dd.loca.inp
    resistorC_linearX.dd.loca.inp.xml
    resistor.dd.equ.exo.4.0
    resistor.dd.equ.exo.4.1
    resistor.dd.equ.exo.4.2
    resistor.dd.equ.exo.4.3
    resistorC_linearX.dd.loca.gold.exo
    tol.exodiff

  #TEST_N CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
  #  ARGS -i resistorC_linearX.dd.loca.inp

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistorC_linearX.dd.loca.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistorC_linearX.dd.loca.inp.log

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistorC_linearX.dd.loca.exo.4."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         resistorC_linearX.dd.loca.exo
         resistorC_linearX.dd.loca.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mixedMode.loca.twoContact PROPERTY LABELS mixedMode.loca.twoContact
    mixedmode)
ENDIF()
