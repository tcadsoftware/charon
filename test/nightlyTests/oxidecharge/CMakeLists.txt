INCLUDE(TribitsAddAdvancedTest)

##########################################################
# These two examples test the bulk oxide charge capability
##########################################################

SET(testName charon_mp_oxide.lapl.bq-)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    oxide.jou
    oxide.exo.4.0
    oxide.exo.4.1
    oxide.exo.4.2
    oxide.exo.4.3
    oxide.lapl.inp
    oxide.lapl.bq-.inp
    oxide.lapl.bq-.gold.exo
    oxide.lapl.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.lapl.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.lapl.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.lapl.bq-.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.lapl.bq-.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_5  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "oxide.lapl.bq-.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # Compare the exodus
  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         oxide.lapl.exodiff
         oxide.lapl.bq-.exo
         oxide.lapl.bq-.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_oxide.lapl.bq- PROPERTY LABELS mp_oxide.lapl.bq-
    nightly)
ENDIF()


##########################################################
SET(testName charon_mp_oxide.sg.lapl.bq)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    oxide.jou
    oxide.exo.4.0
    oxide.exo.4.1
    oxide.exo.4.2
    oxide.exo.4.3
    oxide.lapl.inp
    oxide.sg.lapl.bq.inp
    oxide.sg.lapl.bq.gold.exo
    oxide.lapl.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.lapl.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.lapl.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.sg.lapl.bq.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.sg.lapl.bq.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_5  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "oxide.sg.lapl.bq.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # Compare the exodus
  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         oxide.lapl.exodiff
         oxide.sg.lapl.bq.exo
         oxide.sg.lapl.bq.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_oxide.sg.lapl.bq PROPERTY LABELS mp_oxide.sg.lapl.bq
    nightly)
ENDIF()


##########################################################
SET(testName charon_mp_oxide.lapl.bq-.loca)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    oxide.lapl.exo.4.0
    oxide.lapl.exo.4.1
    oxide.lapl.exo.4.2
    oxide.lapl.exo.4.3
    oxide.lapl.bq-.loca.inp
    oxide.lapl.bq-.loca.gold.exo
    oxide.lapl.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.lapl.bq-.loca.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.lapl.bq-.loca.inp.xml 
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_3  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "oxide.lapl.bq-.loca.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # Compare the exodus
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         oxide.lapl.exodiff
         oxide.lapl.bq-.loca.exo
         oxide.lapl.bq-.loca.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_oxide.lapl.bq-.loca PROPERTY LABELS mp_oxide.lapl.bq-.loca
    nightly)
ENDIF()


##########################################################
SET(testName charon_mp_oxide.sg.lapl.bq.loca)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    oxide.lapl.exo.4.0
    oxide.lapl.exo.4.1
    oxide.lapl.exo.4.2
    oxide.lapl.exo.4.3
    oxide.sg.lapl.bq.loca.inp
    oxide.sg.lapl.bq.loca.gold.exo
    oxide.lapl.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i oxide.sg.lapl.bq.loca.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=oxide.sg.lapl.bq.loca.inp.xml 
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_3  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "oxide.sg.lapl.bq.loca.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # Compare the exodus
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         oxide.lapl.exodiff
         oxide.sg.lapl.bq.loca.exo
         oxide.sg.lapl.bq.loca.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_oxide.sg.lapl.bq.loca PROPERTY LABELS mp_oxide.sg.lapl.bq.loca
    nightly)
ENDIF()

