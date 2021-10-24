INCLUDE(TribitsAddAdvancedTest)

####################################################
# PN 1D diode FEM GaN mobility tests
SET(testName charon_mp_ganPN.dd_albrecht)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TIMEOUT 900

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_albr_dd.yaml
    gan_pndiode_albr_dd_gold.exo

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_albr_dd.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_albr_dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_lf_dd.yaml
    gan_pndiode_fara_lf_dd_gold.exo

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_lf_dd.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_fara_lf_dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan.exodiff
    gan_pndiode_dd.gen.4.0
    gan_pndiode_dd.gen.4.1
    gan_pndiode_dd.gen.4.2
    gan_pndiode_dd.gen.4.3
    gan_pndiode_fara_hf_dd.yaml
    gan_pndiode_fara_hf_dd_gold.exo

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan_pndiode_fara_hf_dd.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan_pndiode_fara_hf_dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TIMEOUT 900

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan-diode.dd.reverse-sweep0.exo.4.0
    gan-diode.dd.reverse-sweep0.exo.4.1
    gan-diode.dd.reverse-sweep0.exo.4.2
    gan-diode.dd.reverse-sweep0.exo.4.3
    gan-pin-diode.dd.exo.4.0
    gan-pin-diode.dd.exo.4.1
    gan-pin-diode.dd.exo.4.2
    gan-pin-diode.dd.exo.4.3
    gan-pin-diode.dd.reverse-sweep.yaml
    gan-pin-diode.dd.reverse-sweep0.yaml
    opto_generation.dat
    gan_pin.exodiff
    gold_gan-diode.dd.reverse-sweep.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY

  # run the reverse LOCA sweep
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=gan-pin-diode.dd.reverse-sweep.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pin_GaN.reverse.output

  # epu auto the result to prepare for exodiff
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "gan-diode.dd.reverse-sweep.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         gan_pin.exodiff
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
# Interpreter test
SET(testName charon_mp_ganPN.dd_albrecht_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_albr_dd.inp
    gan_pndiode_albr_dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i gan_pndiode_albr_dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS gan_pndiode_albr_dd.inp.xml gan_pndiode_albr_dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_albrecht_input PROPERTY LABELS
    ganPN.dd_albrecht inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_ganPN.dd_farahmand_lf_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_fara_lf_dd.inp
    gan_pndiode_fara_lf_dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i gan_pndiode_fara_lf_dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS gan_pndiode_fara_lf_dd.inp.xml gan_pndiode_fara_lf_dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_farahmand_lf_input PROPERTY LABELS
    ganPN.dd_farahmand_lf inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_ganPN.dd_farahmand_hf_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan_pndiode_fara_hf_dd.inp
    gan_pndiode_fara_hf_dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i gan_pndiode_fara_hf_dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS gan_pndiode_fara_hf_dd.inp.xml gan_pndiode_fara_hf_dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_ganPN.dd_farahmand_hf_input PROPERTY LABELS
    ganPN.dd_farahmand_hf inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pin_GaN.reverse_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    gan-pin-diode.dd.reverse-sweep.inp
    gan-pin-diode.dd.reverse-sweep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i gan-pin-diode.dd.reverse-sweep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS gan-pin-diode.dd.reverse-sweep.inp.xml gan-pin-diode.dd.reverse-sweep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pin_GaN.reverse_input PROPERTY LABELS
    pin_GaN.reverse inputVerification)
ENDIF()
