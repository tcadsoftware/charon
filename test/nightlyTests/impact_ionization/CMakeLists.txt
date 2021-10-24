INCLUDE(TribitsAddAdvancedTest)

####################################################
# Crowell-Sze avalanche 
SET(testName charon_mp_CrowellSzeAva)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si_diode.brk.crowellsze.dd.xml
    si_diode.brk.pre.crowellsze.dd.xml
    si_diode.brk.pre.crowellsze.dd.exo.4.0
    si_diode.brk.pre.crowellsze.dd.exo.4.1
    si_diode.brk.pre.crowellsze.dd.exo.4.2
    si_diode.brk.pre.crowellsze.dd.exo.4.3
    opto_generation.dat
    si-diode01.jou
    si-diode.nlp.xml
    si-diode.dd.V0.0.xml
    si_diode.brk.crowellsze.dd.gold.exo
    si_diode.brk.crowellsze.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si_diode.brk.crowellsze.dd.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si_diode.brk.crowellsze.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si_diode.brk.crowellsze.exodiff
         si_diode.brk.crowellsze.dd.exo
         si_diode.brk.crowellsze.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_CrowellSzeAva PROPERTY LABELS
    CrowellSzeAva nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_CrowellSzeAva_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si_diode.brk.crowellsze.dd.inp
    si_diode.brk.crowellsze.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i si_diode.brk.crowellsze.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS si_diode.brk.crowellsze.dd.inp.xml si_diode.brk.crowellsze.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_CrowellSzeAva_input PROPERTY LABELS
    CrowellSzeAva inputVerification)
ENDIF()
