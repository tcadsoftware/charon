INCLUDE(TribitsAddAdvancedTest)

####################################################
#  PN 1D longbase diode HB DD for CV
####################################################

SET(testName charon_mp.hb-dd.pndiode.longbase.cv)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode_dd_equ.exo.4.0
    pndiode_dd_equ.exo.4.1
    pndiode_dd_equ.exo.4.2
    pndiode_dd_equ.exo.4.3
    laux_data_longbase_s3a.txt
    pndiode.hb-dd.loca.longbase.inp
    gold_pndiode-longbase_hb-dd_sweep.exo
    currents-loca.longbase.dat
    compare-hb-current.py


  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.hb-dd.loca.longbase.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-dd.loca.longbase.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-dd.loca.longbase.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode-longbase_hb-dd_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 CMND python3 NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode_dd_equ.exo.4.0
    pndiode_dd_equ.exo.4.1
    pndiode_dd_equ.exo.4.2
    pndiode_dd_equ.exo.4.3
    laux_data_shortbase_s3a.txt
    pndiode.hb-dd.loca.shortbase.inp
    gold_pndiode-shortbase_hb-dd_sweep.exo
    currents-loca.shortbase.dat
    compare-hb-current.py


  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.hb-dd.loca.shortbase.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-dd.loca.shortbase.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-dd.loca.shortbase.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode-shortbase_hb-dd_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 CMND python3 NOEXEPREFIX
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
