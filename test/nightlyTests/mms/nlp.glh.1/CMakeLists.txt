INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run nlp
SET(testName charon_mp_mms_nlp_glh_1)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mms-nlp-glh-1_diode.yaml
    mms-nlp-glh-1_solution.gold.exo
    mms-nlp-glh-1.exodiff
    mms-nlp-glh-1_diode.gen.4.0
    mms-nlp-glh-1_diode.gen.4.1
    mms-nlp-glh-1_diode.gen.4.2
    mms-nlp-glh-1_diode.gen.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=mms-nlp-glh-1_diode.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "mms-nlp-glh-1_solution.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         mms-nlp-glh-1.exodiff
         mms-nlp-glh-1_solution.exo
         mms-nlp-glh-1_solution.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mms_nlp_glh_1 PROPERTY LABELS mms_nlp_glh_1 nightly)
ENDIF()
