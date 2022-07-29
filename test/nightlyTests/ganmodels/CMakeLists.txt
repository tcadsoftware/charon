INCLUDE(TribitsAddAdvancedTest)

####################################################
# PN 1D diode FEM GaN mobility tests
SET(testName charon_mp_ganPN.dd_albrecht)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_albr_dd.inp
    gan_pndiode_albr_dd_gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_albr_dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_albr_dd.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_albr_dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         gan.exodiff
         gan_pndiode_albr_dd.exo
         gan_pndiode_albr_dd_gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_albrecht PROPERTY LABELS
    ganPN.dd_albrecht nightly)
ENDIF()

SET(testName charon_mp_ganPN.dd_farahmand_lf)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_lf_dd.inp
    gan_pndiode_fara_lf_dd_gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_fara_lf_dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_lf_dd.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_fara_lf_dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         gan.exodiff
         gan_pndiode_fara_lf_dd.exo
         gan_pndiode_fara_lf_dd_gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_farahmand_lf  PROPERTY LABELS
    ganPN.dd_farahmand_lf nightly)
ENDIF()
SET(testName charon_mp_ganPN.dd_farahmand_hf)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_hf_dd.inp
    gan_pndiode_fara_hf_dd_gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan_pndiode_fara_hf_dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_hf_dd.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_fara_hf_dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         gan.exodiff
         gan_pndiode_fara_hf_dd.exo
         gan_pndiode_fara_hf_dd_gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_farahmand_hf  PROPERTY LABELS
    ganPN.dd_farahmand_hf nightly)
ENDIF()


####################################################
# PIN GaN diode test: breakdown current profile
SET(testName charon_mp_pin_GaN.reverse)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan-diode.dd.reverse-sweep0.exo.4.0
    gan-diode.dd.reverse-sweep0.exo.4.1
    gan-diode.dd.reverse-sweep0.exo.4.2
    gan-diode.dd.reverse-sweep0.exo.4.3
    gan-pin-diode.dd.exo.4.0
    gan-pin-diode.dd.exo.4.1
    gan-pin-diode.dd.exo.4.2
    gan-pin-diode.dd.exo.4.3
    gan-pin-diode.dd.reverse-sweep.inp
    gan-pin-diode.dd.reverse-sweep0.yaml
    opto_generation.dat
    gan_pin.exodiff
    gold_gan-diode.dd.reverse-sweep.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i gan-pin-diode.dd.reverse-sweep.inp

  # Clean up any existing current file so this run starts with a new file
  TEST_2 CMND rm ARGS -f currents-loca.dat
    PASS_ANY

  # run the reverse LOCA sweep
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan-pin-diode.dd.reverse-sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pin_GaN.reverse.output

  # epu auto the result to prepare for exodiff
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan-diode.dd.reverse-sweep.exo.4."
    NUM_MPI_PROCS 1

  # compare results
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         gan_pin.exodiff
         -steps
         last
         gan-diode.dd.reverse-sweep.exo
         gold_gan-diode.dd.reverse-sweep.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pin_GaN.reverse PROPERTY LABELS
     pin_GaN.reverse nightly debugexclude)
ENDIF()

####################################################
# GaN AlGaN interpolate bandgap, permittivity and DOS
SET(testName charon_mp_GaN_AlGaN)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    ganalgan.gen.4.0
    ganalgan.gen.4.1
    ganalgan.gen.4.2
    ganalgan.gen.4.3
    ganalgan.nlp.inp.yaml
    ganalgan.exodiff
    gold_ganalgan.nlp.exo

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=ganalgan.nlp.inp.yaml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "ganalgan.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         ganalgan.exodiff
         ganalgan.nlp.exo
         gold_ganalgan.nlp.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_GaN_AlGaN  PROPERTY LABELS
    GaN_AlGaN nightly)
ENDIF()
