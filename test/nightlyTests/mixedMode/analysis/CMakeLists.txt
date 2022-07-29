
INCLUDE(TribitsAddAdvancedTest)

####################################################
# Run manufactured and mixed mode problems, then compare
SET(testName charon_mp_mixedMode.analysis)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    cmosInverter_rLoad/cmos_inverter_load.cir
    cmosInverter_rLoad/cmos_inverter_rLoad.dd.gold.exo
    cmosInverter_rLoad/cmos_inverter_rLoad.dd.inp

    pullup_pulldown_res/nmos_pullup.cir
    pullup_pulldown_res/nmos_rPullup.dd.gold.exo
    pullup_pulldown_res/nmos_rPullup.dd.inp
    pullup_pulldown_res/pmos_pulldown.cir
    pullup_pulldown_res/pmos_rPulldown.dd.gold.exo
    pullup_pulldown_res/pmos_rPulldown.dd.inp

    resistor.dd.equ.exo
    resistor.dd.equ.exo.4.0
    resistor.dd.equ.exo.4.1
    resistor.dd.equ.exo.4.2
    resistor.dd.equ.exo.4.3
    tol.exodiff

### begin cmos inverter rload test ###

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i cmos_inverter_rLoad.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=cmos_inverter_rLoad.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE cmos_inverter_rLoad.dd.inp.log

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "cmos_inverter_rLoad.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         cmos_inverter_rLoad.dd.exo
         cmos_inverter_rLoad.dd.gold.exo
    NUM_MPI_PROCS 1

### end cmos inverter rload test ###

### begin nmos pullup test ###

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmos_rPullup.dd.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmos_rPullup.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE nmos_rPullup.dd.inp.log

  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmos_rPullup.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         nmos_rPullup.dd.exo
         nmos_rPullup.dd.gold.exo
    NUM_MPI_PROCS 1

### end nmos pullup test ###

### begin pmos pulldown test ###

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pmos_rPulldown.dd.inp

  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pmos_rPulldown.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pmos_rPulldown.dd.inp.log

  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pmos_rPulldown.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         tol.exodiff
         pmos_rPulldown.dd.exo
         pmos_rPulldown.dd.gold.exo
    NUM_MPI_PROCS 1

### end pmos pulldown test ###

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mixedMode.analysis PROPERTY LABELS mixedMode.analysis
    mixedmode)
ENDIF()
