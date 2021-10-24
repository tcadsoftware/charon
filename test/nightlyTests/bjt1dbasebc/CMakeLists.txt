INCLUDE(TribitsAddAdvancedTest)

####################################################
# NPN 1D BJT test
SET(testName charon_mp_2n2222.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.yaml
    2n2222.dd.equ.exo.4.0
    2n2222.dd.equ.exo.4.1
    2n2222.dd.equ.exo.4.2
    2n2222.dd.equ.exo.4.3
    2n2222.dd.gold.exo
    2n2222.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=2n2222.dd.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "2n2222.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 2n2222.exodiff 2n2222.dd.exo 2n2222.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd PROPERTY LABELS 2n2222.dd.base
    nightly)
ENDIF()

####################################################
# NPN 1D BJT test
SET(testName charon_mp_2n2222.sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.sg.yaml
    2n2222.dd.equ.exo.4.0
    2n2222.dd.equ.exo.4.1
    2n2222.dd.equ.exo.4.2
    2n2222.dd.equ.exo.4.3
    2n2222.sg.gold.exo
    2n2222.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=2n2222.sg.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "2n2222.sg.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 2n2222.exodiff 2n2222.sg.exo 2n2222.sg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.sg PROPERTY LABELS 2n2222.sg
    nightly)
ENDIF()

####################################################
# PNP 1D BJT test
SET(testName charon_mp_2c2907.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2c2907.dd.yaml
    2c2907.dd.equ.exo.4.0
    2c2907.dd.equ.exo.4.1
    2c2907.dd.equ.exo.4.2
    2c2907.dd.equ.exo.4.3
    2c2907.dd.gold.exo
    2c2907.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=2c2907.dd.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "2c2907.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 2c2907.exodiff 2c2907.dd.exo 2c2907.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2c2907.dd PROPERTY LABELS 2c2907.dd
    nightly)
ENDIF()

####################################################
# PNP 1D BJT test
SET(testName charon_mp_2c2907.sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2c2907.sg.yaml
    2c2907.dd.equ.exo.4.0
    2c2907.dd.equ.exo.4.1
    2c2907.dd.equ.exo.4.2
    2c2907.dd.equ.exo.4.3
    2c2907.sg.gold.exo
    2c2907.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=2c2907.sg.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "2c2907.sg.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 2c2907.exodiff 2c2907.sg.exo 2c2907.sg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2c2907.sg PROPERTY LABELS 2c2907.sg
    nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2c2907.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2c2907.dd.inp
    2c2907.dd.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2c2907.dd.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2c2907.dd.inp.xml 2c2907.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2c2907.dd_input PROPERTY LABELS
    2c2907.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2c2907.sg_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2c2907.sg.inp
    2c2907.sg.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2c2907.sg.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2c2907.sg.inp.xml 2c2907.sg.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2c2907.sg_input PROPERTY LABELS
    2c2907.sg_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.inp
    2n2222.dd.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.inp.xml 2n2222.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd_input PROPERTY LABELS
    2n2222.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.sg_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.sg.inp
    2n2222.sg.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.sg.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.sg.inp.xml 2n2222.sg.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.sg_input PROPERTY LABELS
    2n2222.sg_input inputVerification)
ENDIF()
