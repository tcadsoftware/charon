INCLUDE(TribitsAddAdvancedTest)

####################################################
# SUPG-FEM test: 2nd order current calculation
SET(testName charon_mp_pndiode.dd.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.current.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.current.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.current.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.current.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.current.output
         --anodeCurrent   -2.72738591e-02
         --cathodeCurrent 2.72738591e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.current PROPERTY LABELS
    pndiode.dd.current nightly)
ENDIF()

####################################################
# EFFPG-FEM test: 2nd order current calculation
SET(testName charon_mp_pndiode.effpg.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.effpg.current.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.effpg.current.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.effpg.current.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.effpg.current.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.effpg.current.output
         --anodeCurrent   -2.72738591e-02
         --cathodeCurrent 2.72738591e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.effpg.current PROPERTY LABELS pndiode.effpg.current nightly)
ENDIF()

####################################################
# Scharfetter-Gummel test: 2nd order current calc.
SET(testName charon_mp_pndiode.sg.current)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.current.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.current.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.current.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.current.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.current.output
         --anodeCurrent   -2.72785150e-02
         --cathodeCurrent 2.72785150e-02
         --tolerance      1e-12
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.current PROPERTY LABELS pndiode.sg.current nightly)
ENDIF()
