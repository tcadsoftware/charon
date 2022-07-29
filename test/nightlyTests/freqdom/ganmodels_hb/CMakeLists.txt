INCLUDE(TribitsAddAdvancedTest)

####################################################
#  Test collection of gan models for HB
####################################################

####################################################
# PN 1D diode FEM GaN mobility tests
SET(testName charon_mp_hb-ganPN.dd_albrecht)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_albr_dd.hb.inp

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_albr_dd.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_albr_dd.hb.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hb-ganPN.dd_albrecht PROPERTY LABELS
    hb-ganPN.dd_albrecht nightly)
ENDIF()

SET(testName charon_mp_hb-ganPN.dd_farahmand_lf)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_lf_dd.hb.inp

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_fara_lf_dd.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_lf_dd.hb.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hb-ganPN.dd_farahmand_lf  PROPERTY LABELS
    hb-ganPN.dd_farahmand_lf nightly)
ENDIF()

SET(testName charon_mp_hb-ganPN.dd_farahmand_hf)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_hf_dd.hb.inp

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_fara_hf_dd.hb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_hf_dd.hb.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hb-ganPN.dd_farahmand_hf  PROPERTY LABELS
    hb-ganPN.dd_farahmand_hf nightly)
ENDIF()

####################################################
# GaN AlGaN interpolate bandgap, permittivity and DOS
#SET(testName charon_mp_hb-GaN_AlGaN)
#TRIBITS_ADD_ADVANCED_TEST(
#  ${testName}
#  CATEGORIES NIGHTLY
#  OVERALL_WORKING_DIRECTORY TEST_NAME
#
#  TEST_0 COPY_FILES_TO_TEST_DIR
#    ganalgan.gen.4.0
#    ganalgan.gen.4.1
#    ganalgan.gen.4.2
#    ganalgan.gen.4.3
#    ganalgan.nlp.hb.inp
#
#  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
#    ARGS -i ganalgan.nlp.hb.inp
#
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=ganalgan.nlp.hb.inp.xml
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  ADDED_TEST_NAME_OUT ${testName}_CREATED
#  )
#IF (${testName}_CREATED)
#  SET_PROPERTY(TEST Charon_charon_mp_hb-GaN_AlGaN  PROPERTY LABELS
#    hb-GaN_AlGaN nightly)
#ENDIF()
