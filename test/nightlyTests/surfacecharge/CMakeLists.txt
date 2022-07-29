INCLUDE(TribitsAddAdvancedTest)

####################################################
SET(testName charon_mp_surfacecharge_dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    surfacecharge.dd.hfmob.inp
    surfacecharge.dd.gold.exo
    surfacecharge.dd.hfmob.exodiff
    surfacecharge.nlp.exo.4.0
    surfacecharge.nlp.exo.4.1
    surfacecharge.nlp.exo.4.2
    surfacecharge.nlp.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i surfacecharge.dd.hfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=surfacecharge.dd.hfmob.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "surfacecharge.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         surfacecharge.dd.hfmob.exodiff
         surfacecharge.dd.exo
         surfacecharge.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_surfacecharge_dd PROPERTY LABELS surfacecharge_dd nightly)
ENDIF()


####################################################
SET(testName charon_mp_surfacecharge_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    surfacecharge.sg.hfmob.inp
    surfacecharge.sg.gold.exo
    surfacecharge.sg.hfmob.exodiff
    surfacecharge.nlp.exo.4.0
    surfacecharge.nlp.exo.4.1
    surfacecharge.nlp.exo.4.2
    surfacecharge.nlp.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i surfacecharge.sg.hfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=surfacecharge.sg.hfmob.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "surfacecharge.sg.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         surfacecharge.sg.hfmob.exodiff
         surfacecharge.sg.exo
         surfacecharge.sg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_surfacecharge_sg PROPERTY LABELS surfacecharge_sg nightly)
ENDIF()


####################################################
SET(testName charon_mp_interfacecharge_dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.ifqm.nlp.inp
    nmosfet.ifqm.dd.equ.inp
    nmosfet.exo.4.0
    nmosfet.exo.4.1
    nmosfet.exo.4.2
    nmosfet.exo.4.3
    nmosfet.ifqm.dd.equ.gold.exo
    nmosfet.ifqm.dd.equ.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.ifqm.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.ifqm.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.ifqm.dd.equ.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.ifqm.dd.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.ifqm.dd.equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.ifqm.dd.equ.exodiff
         nmosfet.ifqm.dd.equ.exo
         nmosfet.ifqm.dd.equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_interfacecharge_dd PROPERTY LABELS interfacecharge_dd nightly)
ENDIF()


####################################################
SET(testName charon_mp_interfacecharge_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.ifqp.nlp.inp
    nmosfet.ifqp.sg.equ.inp
    nmosfet.exo.4.0
    nmosfet.exo.4.1
    nmosfet.exo.4.2
    nmosfet.exo.4.3
    nmosfet.ifqp.sg.equ.gold.exo
    nmosfet.ifqp.sg.equ.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.ifqp.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.ifqp.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.ifqp.sg.equ.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.ifqp.sg.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.ifqp.sg.equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.ifqp.sg.equ.exodiff
         nmosfet.ifqp.sg.equ.exo
         nmosfet.ifqp.sg.equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_interfacecharge_sg PROPERTY LABELS interfacecharge_sg nightly)
ENDIF()


####################################################
SET(testName charon_mp_surfacetrap_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.nlp.inp
    nmosfet.sg.equ.inp
    nmosfet.sg.it.equ.inp
    nmosfet.sg.it.loca.inp
    nmosfet.jou
    nmosfet.exo.4.0
    nmosfet.exo.4.1
    nmosfet.exo.4.2
    nmosfet.exo.4.3
    nmosfet.sg.it.loca.gold.exo
    nmosfet.sg.it.loca.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.sg.equ.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.sg.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.sg.it.equ.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.sg.it.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_7 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.sg.it.loca.inp

  TEST_8 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.sg.it.loca.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.sg.it.loca.exo.4."
    NUM_MPI_PROCS 1

  TEST_10 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.sg.it.loca.exodiff
         nmosfet.sg.it.loca.exo
         nmosfet.sg.it.loca.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_surfacetrap_sg PROPERTY LABELS surfacetrap_sg nightly debugexclude)
ENDIF()


#########################################################################################
SET(testName charon_mp_surfacerecomb_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.sg.sr.equ.inp
    nmosfet.sg.equ.exo.4.0
    nmosfet.sg.equ.exo.4.1
    nmosfet.sg.equ.exo.4.2
    nmosfet.sg.equ.exo.4.3
    nmosfet.sg.sr.equ.gold.exo
    nmosfet.sg.it.loca.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.sg.sr.equ.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.sg.sr.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.sg.sr.equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.sg.it.loca.exodiff
         nmosfet.sg.sr.equ.exo
         nmosfet.sg.sr.equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_surfacerecomb_sg PROPERTY LABELS surfacerecomb_sg nightly)
ENDIF()


#########################################################################################
SET(testName charon_mp_surfacetrap_dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.dd.it.a1e11.equ.inp
    nmosfet.sg.equ.exo.4.0
    nmosfet.sg.equ.exo.4.1
    nmosfet.sg.equ.exo.4.2
    nmosfet.sg.equ.exo.4.3
    nmosfet.dd.it.a1e11.equ.gold.exo
    nmosfet.sg.it.loca.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i nmosfet.dd.it.a1e11.equ.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.it.a1e11.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.it.a1e11.equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.sg.it.loca.exodiff
         nmosfet.dd.it.a1e11.equ.exo
         nmosfet.dd.it.a1e11.equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_surfacetrap_dd PROPERTY LABELS surfacetrap_dd nightly debugexclude)
ENDIF()
