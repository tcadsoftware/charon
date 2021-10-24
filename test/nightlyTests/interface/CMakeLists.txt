INCLUDE(TribitsAddAdvancedTest)

####################################################
# Interface test: Simple
SET(testName charon_mp_interface)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    interface.yaml
    simple-ef-coupled.yaml
    interface.exo.4.0
    interface.exo.4.1
    interface.exo.4.2
    interface.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=interface.yaml --current --current-ho
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE interface.output

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_interface  PROPERTY LABELS interface
    nightly)
ENDIF()
