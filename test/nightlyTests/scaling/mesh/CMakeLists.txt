INCLUDE(TribitsAddAdvancedTest)

####################################################
# Check coordinate scaling
SET(testName charon_mp.mesh-scaling)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mesh-scale.sg.nlp.yaml
    mesh-scale.exodiff
    mesh-scale.gen.4.0
    mesh-scale.gen.4.1
    mesh-scale.gen.4.2
    mesh-scale.gen.4.3
    mesh-scale.nlp.gold.exo

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=mesh-scale.sg.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "mesh-scale.nlp.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f mesh-scale.exodiff mesh-scale.nlp.exo mesh-scale.nlp.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.mesh-scaling PROPERTY LABELS mesh-scaling
    nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_mesh-scaling_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mesh-scale.sg.nlp.inp
    mesh-scale.sg.nlp.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i mesh-scale.sg.nlp.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS mesh-scale.sg.nlp.inp.xml mesh-scale.sg.nlp.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mesh-scaling_input PROPERTY LABELS
    mesh-scaling_input inputVerification)
ENDIF()
