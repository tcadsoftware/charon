INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run nlp
SET(testName charon_mp_mms_nlp_glh_1)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mms-nlp-glh-1_diode.inp
    mms-nlp-glh-1_solution.gold.exo
    mms-nlp-glh-1.exodiff
    mms-nlp-glh-1_diode.gen.4.0
    mms-nlp-glh-1_diode.gen.4.1
    mms-nlp-glh-1_diode.gen.4.2
    mms-nlp-glh-1_diode.gen.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i mms-nlp-glh-1_diode.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=mms-nlp-glh-1_diode.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "mms-nlp-glh-1_solution.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
