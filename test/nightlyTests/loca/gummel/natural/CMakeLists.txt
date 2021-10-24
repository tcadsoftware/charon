INCLUDE(TribitsAddAdvancedTest)

####################################################
# SUPG-FEM test: 1st order current calculation
SET(testName charon_mp_pndiode.loca.gummel.natural)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.gen.4.0
    pndiode.gen.4.1
    pndiode.gen.4.2
    pndiode.gen.4.3
    pndiode.nlp.yaml
    pndiode.nlp-result.gold.exo
    pndiode.nlp.exodiff
    pndiode.loca.dd.yaml
    pndiode.dd-result.loca.gold.exo
    pndiode.loca.exodiff
    currents-loca.dat.gold

  # run the NLP problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.nlp.output

  # Clean up any existing current file so this run starts with a new file
  TEST_2 CMND rm ARGS -f currents-loca.dat
    PASS_ANY

  # test the correctness of the nlp exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.nlp-result.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         pndiode.nlp-result.exo
         pndiode.nlp-result.gold.exo
    NUM_MPI_PROCS 1

  # run the LOCA gummel sweep
  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.loca.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.loca.dd.output

  # test the output current values
  TEST_6 CMND python3
    ARGS
      ${CMAKE_SOURCE_DIR}/scripts/utils/compare-text-value-files.py
      currents-loca.dat
      currents-loca.dat.gold
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the loca exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.dd-result.loca.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.loca.exodiff
         pndiode.dd-result.loca.exo
         pndiode.dd-result.loca.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
)
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_${testName} PROPERTY LABELS
    pndiode.loca.gummel.natural nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.loca.nlp.natural_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.inp
    pndiode.nlp.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.nlp.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.nlp.inp.xml pndiode.nlp.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.loca.nlp.natural_input PROPERTY LABELS
    pndiode.loca.nlp.natural_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.loca.gummel.natural_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.loca.dd.inp
    pndiode.loca.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.loca.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.loca.dd.inp.xml pndiode.loca.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.loca.gummel.natural_input PROPERTY LABELS
    pndiode.loca.gummel.natural_input inputVerification)
ENDIF()
