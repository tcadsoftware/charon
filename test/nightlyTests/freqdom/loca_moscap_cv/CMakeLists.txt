INCLUDE(TribitsAddAdvancedTest)

####################################################
#  moscap 1D HB small signal 1e0Hz oxide sweep for CV
####################################################

SET(testName charon_mp.hb.loca-moscap-cv.ox)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    moscap.dd.equ.exo.4.0
    moscap.dd.equ.exo.4.1
    moscap.dd.equ.exo.4.2
    moscap.dd.equ.exo.4.3
    moscap.hb.locaOxide.1e0Hz.inp
    moscap.hb.locaOxide.1e0Hz.gold.exo
    currents-loca.locaOxide.1e0Hz.gold.dat
    compare-hb-current.py

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i moscap.hb.locaOxide.1e0Hz.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=moscap.hb.locaOxide.1e0Hz.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE moscap.hb.locaOxide.1e0Hz.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "moscap.hb.locaOxide.1e0Hz.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 CMND python3 NOEXEPREFIX
    NOEXESUFFIX
    ARGS
      compare-hb-current.py
      currents-loca.locaOxide.1e0Hz.gold.dat
      currents-loca.dat
    PASS_REGULAR_EXPRESSION "PASS"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb.loca-moscap-cv.ox  PROPERTY LABELS
    hb.loca-moscap-cv.ox nightly)
ENDIF()




####################################################
#  moscap 1D HB small signal 1e6Hz substrate sweep for CV
####################################################

SET(testName charon_mp.hb.loca-moscap-cv.sub)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    moscap.dd.equ.exo.4.0
    moscap.dd.equ.exo.4.1
    moscap.dd.equ.exo.4.2
    moscap.dd.equ.exo.4.3
    moscap.hb.locaSubstrate.1e6Hz.inp
    moscap.hb.locaSubstrate.1e6Hz.gold.exo
    currents-loca.locaSubstrate.1e6Hz.gold.dat
    compare-hb-current.py

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i moscap.hb.locaSubstrate.1e6Hz.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=moscap.hb.locaSubstrate.1e6Hz.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE moscap.hb.locaSubstrate.1e6Hz.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "moscap.hb.locaSubstrate.1e6Hz.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 CMND python3 NOEXEPREFIX
    NOEXESUFFIX
    ARGS
      compare-hb-current.py
      currents-loca.locaSubstrate.1e6Hz.gold.dat
      currents-loca.dat
    PASS_REGULAR_EXPRESSION "PASS"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb.loca-moscap-cv.sub  PROPERTY LABELS
    hb.loca-moscap-cv.sub nightly)
ENDIF()
