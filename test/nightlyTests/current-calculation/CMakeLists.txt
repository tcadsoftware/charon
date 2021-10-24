INCLUDE(TribitsAddAdvancedTest)

####################################################
# SUPG-FEM test: 1st order current calculation
SET(testName charon_mp_pndiode.dd.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.current.yaml
    current-calculation.exodiff
    result_pndiode.dd.current.gold.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.current.yaml --current --current-lo
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.current.output

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_current.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         current-calculation.exodiff
         result_pndiode_dd_current.exo
         result_pndiode.dd.current.gold.exo
    NUM_MPI_PROCS 1

  # now test that the current is correct
  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.current.output
         --anodeCurrent   -1.26928602e-03
         --cathodeCurrent 1.26928978e-03
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.current PROPERTY LABELS
    pndiode.dd.current nightly)
ENDIF()

####################################################
# SUPG-FEM test: 2nd order current calculation
SET(testName charon_mp_pndiode.dd.hocurrent)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.hocurrent.yaml
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.hocurrent.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.hocurrent.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.hocurrent.output
         --anodeCurrent   -2.72738591e-02
         --cathodeCurrent 2.72738591e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.hocurrent PROPERTY LABELS
    pndiode.dd.hocurrent nightly)
ENDIF()

####################################################
# EFFPG-FEM test: 1st order current calculation
SET(testName charon_mp_pndiode.effpg.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.effpg.current.yaml
    current-calculation.exodiff
    result_pndiode.effpg.current.gold.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.effpg.current.yaml --current --current-lo
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.effpg.current.output

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_effpg_current.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f current-calculation.exodiff
         result_pndiode_effpg_current.exo
         result_pndiode.effpg.current.gold.exo
    NUM_MPI_PROCS 1

  # now test that the current is correct
  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.effpg.current.output
         --anodeCurrent   -1.26928602e-03
         --cathodeCurrent 1.26928978e-03
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.effpg.current PROPERTY LABELS pndiode.effpg.current nightly)
ENDIF()

####################################################
# EFFPG-FEM test: 2nd order current calculation
SET(testName charon_mp_pndiode.effpg.hocurrent)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.effpg.hocurrent.yaml
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.effpg.hocurrent.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.effpg.hocurrent.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.effpg.hocurrent.output
         --anodeCurrent   -2.72738591e-02
         --cathodeCurrent 2.72738591e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.effpg.hocurrent PROPERTY LABELS pndiode.effpg.hocurrent nightly)
ENDIF()

####################################################
# Scharfetter-Gummel test: 1st order current calc.
SET(testName charon_mp_pndiode.sg.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.current.yaml
    current-calculation.exodiff
    result_pndiode.sg.current.gold.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.current.yaml --current --current-lo
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.current.output

  # test the consistency of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_sg_current.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f current-calculation.exodiff
         result_pndiode_sg_current.exo
         result_pndiode.sg.current.gold.exo
    NUM_MPI_PROCS 1

  # now test that the current is correct
  TEST_4 CMND python3 NOEXEPREFIX NOEXESUFFIX
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.current.output
         --anodeCurrent   -1.26977150e-03
         --cathodeCurrent 1.26977531e-03
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.current PROPERTY LABELS pndiode.sg.current nightly)
ENDIF()

####################################################
# Scharfetter-Gummel test: 2nd order current calc.
SET(testName charon_mp_pndiode.sg.hocurrent)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.hocurrent.yaml
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.hocurrent.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.hocurrent.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.hocurrent.output
         --anodeCurrent   -2.72785150e-02
         --cathodeCurrent 2.72785150e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.hocurrent PROPERTY LABELS pndiode.sg.hocurrent nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.current_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.current.inp
    pndiode.dd.current.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.current.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.current.inp.xml pndiode.dd.current.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.current_input PROPERTY LABELS
    pndiode.dd.current_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.hocurrent_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.hocurrent.inp
    pndiode.dd.hocurrent.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.hocurrent.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.hocurrent.inp.xml pndiode.dd.hocurrent.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.hocurrent_input PROPERTY LABELS
    pndiode.dd.hocurrent inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.effpg.current_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.effpg.current.inp
    pndiode.effpg.current.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.effpg.current.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.effpg.current.inp.xml pndiode.effpg.current.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.effpg.current_input PROPERTY LABELS
    pndiode.effpg.current inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.effpg.hocurrent_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.effpg.hocurrent.inp
    pndiode.effpg.hocurrent.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.effpg.hocurrent.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.effpg.hocurrent.inp.xml pndiode.effpg.hocurrent.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.effpg.hocurrent_input PROPERTY LABELS
    pndiode.effpg.hocurrent inputVerification)
ENDIF()







####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.current_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.current.inp
    pndiode.sg.current.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.current.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.current.inp.xml pndiode.sg.current.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.current_input PROPERTY LABELS
    pndiode.sg.current inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.hocurrent_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.hocurrent.inp
    pndiode.sg.hocurrent.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.hocurrent.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.hocurrent.inp.xml pndiode.sg.hocurrent.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.hocurrent_input PROPERTY LABELS
    pndiode.sg.hocurrent inputVerification)
ENDIF()
