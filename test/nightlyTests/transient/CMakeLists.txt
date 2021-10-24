INCLUDE(TribitsAddAdvancedTest)

SET(testName charon_mp_pndiode.dd.lin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.lin.fixstep.yaml
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_dd_lin.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.lin.fixstep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.lin.fixstep PROPERTY LABELS
    pndiode.dd.lin.fixstep nightly)
ENDIF()

SET(testName charon_mp_pndiode.sg.lin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.fixstep.yaml
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_sg_lin.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.lin.fixstep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.fixstep PROPERTY LABELS
    pndiode.sg.lin.fixstep nightly)
ENDIF()

SET(testName charon_mp_pndiode.sg.lin.varstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.varstep.yaml
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_sg_lin_varstep.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.lin.varstep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.varstep PROPERTY LABELS
    pndiode.sg.lin.varstep nightly)
ENDIF()

SET(testName charon_mp_pndiode.sg.sin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.sin.fixstep.yaml
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3

  TEST_1 CMND rm ARGS -f current_time_sg_sin.csv
    PASS_ANY

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.sin.fixstep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.sin.fixstep PROPERTY LABELS
    pndiode.sg.sin.fixstep  nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.lin.fixstep_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.lin.fixstep.inp
    pndiode.dd.lin.fixstep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.lin.fixstep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.lin.fixstep.inp.xml pndiode.dd.lin.fixstep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.lin.fixstep_input PROPERTY LABELS
    pndiode.dd.lin.fixstep_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.lin.fixstep_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.fixstep.inp
    pndiode.sg.lin.fixstep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.lin.fixstep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.lin.fixstep.inp.xml pndiode.sg.lin.fixstep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.fixstep_input PROPERTY LABELS
    pndiode.sg.lin.fixstep_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.lin.varstep_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.varstep.inp
    pndiode.sg.lin.varstep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.lin.varstep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.lin.varstep.inp.xml pndiode.sg.lin.varstep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.varstep_input PROPERTY LABELS
    pndiode.sg.lin.varstep_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.sin.fixstep_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.sin.fixstep.inp
    pndiode.sg.sin.fixstep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.sin.fixstep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.sin.fixstep.inp.xml pndiode.sg.sin.fixstep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.sin.fixstep_input PROPERTY LABELS
    pndiode.sg.sin.fixstep_input inputVerification)
ENDIF()
