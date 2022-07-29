INCLUDE(TribitsAddAdvancedTest)

####################################################
# Displacement current on contacts
SET(testName charon_mp_DispContCurent_nmoscap_SGCVFEM)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc_msh.jou  
    nmosc_msh.gen.4.0
    nmosc_msh.gen.4.1
    nmosc_msh.gen.4.2
    nmosc_msh.gen.4.3
    nmosc.nlp.SGCVFEM.inp
    nmosc.nlp.SGCVFEM.exodiff
    nmosc.nlp.SGCVFEM.gold.exo
    nmosc.dd.SGCVFEM.inp
    nmosc.dd.SGCVFEM.exodiff
    nmosc.dd.SGCVFEM.gold.exo
    nmosc.dd.trans.SGCVFEM.inp
    nmosc.dd.trans.SGCVFEM.exodiff
    nmosc.dd.trans.SGCVFEM.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosc.nlp.SGCVFEM.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.nlp.SGCVFEM.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.nlp.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.nlp.SGCVFEM.exodiff
         nmosc.nlp.SGCVFEM.exo
         nmosc.nlp.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosc.dd.SGCVFEM.inp
    
  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.SGCVFEM.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.SGCVFEM.exodiff
         nmosc.dd.SGCVFEM.exo
         nmosc.dd.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosc.dd.trans.SGCVFEM.inp

  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.trans.SGCVFEM.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.trans.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.trans.SGCVFEM.exodiff
         nmosc.dd.trans.SGCVFEM.exo
         nmosc.dd.trans.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DispContCurent_nmoscap_SGCVFEM PROPERTY LABELS
     DispContCurent_nmoscap_SGCVFEM nightly)
ENDIF()
