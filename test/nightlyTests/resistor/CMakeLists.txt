INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run nlp
SET(testName charon_mp_resistor.nlp)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.nlp.inp
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor.nlp PROPERTY LABELS resistor.nlp
    nightly)
ENDIF()

####################################################
# Run fem-dd
SET(testName charon_mp_resistor.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.dd.inp
    resistor.nlp.exo.4.0
    resistor.nlp.exo.4.1
    resistor.nlp.exo.4.2
    resistor.nlp.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor.dd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor.dd.output
         --anodeCurrent   -8.01088248e-01
         --cathodeCurrent 8.01088230e-01
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor.dd PROPERTY LABELS resistor.dd
    nightly)
ENDIF()
