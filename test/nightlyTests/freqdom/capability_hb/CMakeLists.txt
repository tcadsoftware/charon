INCLUDE(TribitsAddAdvancedTest)

####################################################
#  Test collection of mobility models for HB
####################################################

SET(testName charon_mp.hb-capability.mobility-hfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    accIncmplIoniz.txt
    donIncmplIoniz.txt
    pndiode.dd.ava.hb.inp
    pndiode.dd.hfmob.hb.inp
    pndiode.dd.ioniz.hb.inp
    pndiode.dd.lfmob.hb.inp
    pndiode.dd.recomb.hb.inp
    result_pndiode_equ.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    test1_dop1d.txt
    test2_dop1d.txt
    test3_dop2d.txt
    test4_dop2d.txt
    test_dop.txt

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.hfmob.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.hfmob.hb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.hfmob.hb.inp.xml.out

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-capability.mobility-hfmob  PROPERTY LABELS
    hb-capability.mobility-hfmob nightly)
ENDIF()

SET(testName charon_mp.hb-capability.mobility-lfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    accIncmplIoniz.txt
    donIncmplIoniz.txt
    pndiode.dd.ava.hb.inp
    pndiode.dd.hfmob.hb.inp
    pndiode.dd.ioniz.hb.inp
    pndiode.dd.lfmob.hb.inp
    pndiode.dd.recomb.hb.inp
    result_pndiode_equ.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    test1_dop1d.txt
    test2_dop1d.txt
    test3_dop2d.txt
    test4_dop2d.txt
    test_dop.txt

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.lfmob.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.lfmob.hb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.lfmob.hb.inp.xml.out

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-capability.mobility-lfmob  PROPERTY LABELS
    hb-capability.mobility-lfmob nightly)
ENDIF()

SET(testName charon_mp.hb-capability.mobility-recomb)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    accIncmplIoniz.txt
    donIncmplIoniz.txt
    pndiode.dd.ava.hb.inp
    pndiode.dd.hfmob.hb.inp
    pndiode.dd.ioniz.hb.inp
    pndiode.dd.lfmob.hb.inp
    pndiode.dd.recomb.hb.inp
    result_pndiode_equ.exo
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    test1_dop1d.txt
    test2_dop1d.txt
    test3_dop2d.txt
    test4_dop2d.txt
    test_dop.txt

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.recomb.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.recomb.hb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.recomb.hb.inp.xml.out

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-capability.mobility-recomb  PROPERTY LABELS
    hb-capability.mobility-recomb nightly)
ENDIF()
