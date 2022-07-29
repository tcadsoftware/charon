INCLUDE(TribitsAddAdvancedTest)

####################################################
# Test the FEM NLP equation
SET(testName charon_mp_mosQuantumCorrection_pmos)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    mos.exo
    pmos.nlp.exo
    pmos.sg.gold.exo
    pmos.sg.inp
    pmos.sg.qc.gold.exo
    pmos.sg.qc.inp

  ########################################################
  # classical charge distribution calculation using SG
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pmos.sg.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pmos.sg.inp.xml --current
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         pmos.sg.exo
         pmos.sg.gold.exo
    NUM_MPI_PROCS 1
 
  ########################################################
  # quantum charge distribution calculation using SG and FEM
  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pmos.sg.qc.inp

  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pmos.sg.qc.inp.xml --current
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         pmos.sg.qc.exo
         pmos.sg.qc.gold.exo
    NUM_MPI_PROCS 1


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_mosQuantumCorrection_pmos PROPERTY LABELS mosQuantumCorrection_pmos nightly debugexclude)
ENDIF()
