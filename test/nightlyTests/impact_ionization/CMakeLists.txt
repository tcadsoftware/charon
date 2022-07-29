INCLUDE(TribitsAddAdvancedTest)

####################################################
# Crowell-Sze avalanche 
SET(testName charon_mp_CrowellSzeAva)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si_diode.brk.crowellsze.dd.inp
    si_diode.brk.pre.crowellsze.dd.exo.4.0
    si_diode.brk.pre.crowellsze.dd.exo.4.1
    si_diode.brk.pre.crowellsze.dd.exo.4.2
    si_diode.brk.pre.crowellsze.dd.exo.4.3
    opto_generation.dat
    si_diode.brk.crowellsze.dd.gold.exo
    si_diode.brk.crowellsze.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i si_diode.brk.crowellsze.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si_diode.brk.crowellsze.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si_diode.brk.crowellsze.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -steps
         last
         -f
         si_diode.brk.crowellsze.exodiff
         si_diode.brk.crowellsze.dd.exo
         si_diode.brk.crowellsze.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_CrowellSzeAva PROPERTY LABELS
    CrowellSzeAva nightly debugexclude)
ENDIF()
