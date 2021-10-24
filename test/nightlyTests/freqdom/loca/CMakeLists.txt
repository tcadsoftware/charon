INCLUDE(TribitsAddAdvancedTest)

####################################################
#  PN 1D longbase diode HB DD for CV
####################################################

SET(testName charon_mp.hb-dd.pndiode.longbase.cv)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode_dd_equ.exo.4.0
    pndiode_dd_equ.exo.4.1
    pndiode_dd_equ.exo.4.2
    pndiode_dd_equ.exo.4.3
    laux_data_longbase_s3a.txt
    pndiode-longbase.hb-dd.sweep.yaml
    gold_pndiode-longbase_hb-dd_sweep.exo
    currents-loca.longbase.dat
    compare-hb-current.py

  # run the hb rectifying diode test problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode-longbase.hb-dd.sweep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode-longbase.hb-dd.sweep.out

  # epu auto the result to prepare for exodiff
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode-longbase_hb-dd_sweep.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results to the gold
  TEST_3 CMND python3 NOEXEPREFIX
    NOEXESUFFIX
    ARGS
      compare-hb-current.py
      currents-loca.longbase.dat
      currents-loca.dat
    PASS_REGULAR_EXPRESSION "PASS"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.pndiode.longbase.cv  PROPERTY LABELS
    hb-dd.pndiode.longbase.cv nightly)
ENDIF()




####################################################
#  PN 1D shortbase diode HB DD for CV
####################################################

SET(testName charon_mp.hb-dd.pndiode.shortbase.cv)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode_dd_equ.exo.4.0
    pndiode_dd_equ.exo.4.1
    pndiode_dd_equ.exo.4.2
    pndiode_dd_equ.exo.4.3
    laux_data_shortbase_s3a.txt
    pndiode-shortbase.hb-dd.sweep.yaml
    gold_pndiode-shortbase_hb-dd_sweep.exo
    currents-loca.shortbase.dat
    compare-hb-current.py

  # run the hb rectifying diode test problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode-shortbase.hb-dd.sweep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode-shortbase.hb-dd.sweep.out

  # epu auto the result to prepare for exodiff
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode-shortbase_hb-dd_sweep.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results to the gold
  TEST_3 CMND python3 NOEXEPREFIX
    NOEXESUFFIX
    ARGS
      compare-hb-current.py
      currents-loca.shortbase.dat
      currents-loca.dat
    PASS_REGULAR_EXPRESSION "PASS"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.pndiode.shortbase.cv  PROPERTY LABELS
    hb-dd.pndiode.shortbase.cv nightly)
ENDIF()
