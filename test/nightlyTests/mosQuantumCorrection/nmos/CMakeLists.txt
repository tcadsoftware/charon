INCLUDE(TribitsAddAdvancedTest)

####################################################
# Test the FEM NLP equation
SET(testName charon_mp_mosQuantumCorrection_nmos)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mos.exo
    nmos.nlp.exo
    nmos.sg.gold.exo
    nmos.sg.inp
    nmos.sg.qc.gold.exo
    nmos.sg.qc.inp

  ########################################################
  # classical charge distribution calculation using SG
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i nmos.sg.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmos.sg.inp.xml --current
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         nmos.sg.exo
         nmos.sg.gold.exo
    NUM_MPI_PROCS 1
 
  ########################################################
  # quantum charge distribution calculation using SG and FEM
  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i nmos.sg.qc.inp

  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmos.sg.qc.inp.xml --current
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         nmos.sg.qc.exo
         nmos.sg.qc.gold.exo
    NUM_MPI_PROCS 1


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mosQuantumCorrection_nmos PROPERTY LABELS mosQuantumCorrection_nmos nightly debugexclude)
ENDIF()
