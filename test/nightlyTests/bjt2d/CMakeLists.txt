INCLUDE(TribitsAddAdvancedTest)

####################################################
# Test the FEM NLP equation
SET(testName charon_mp_bjt2d)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    bjt2d.nlp.yaml
    bjt2d.nlp.gold.exo
    bjt2d.gen
    bjt2d.gen.4.0
    bjt2d.gen.4.1
    bjt2d.gen.4.2
    bjt2d.gen.4.3
    bjt2d.dd.equ.yaml
    bjt2d.dd.equ.gold.exo
    bjt2d.dd.current.yaml
    bjt2d.dd.current.gold.exo
    bjt2d.dd.hocurrent.yaml
    bjt2d.dd.hocurrent.gold.exo
    bjt2d.sg.current.yaml
    bjt2d.sg.current.gold.exo
    bjt2d.sg.hocurrent.yaml
    bjt2d.sg.hocurrent.gold.exo

  ########################################################
  # Nonlinear Poisson simulation
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.nlp.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative -tolerance 1.0e-8 bjt2d.nlp.exo bjt2d.nlp.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Equilibrium FEM simulation
  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.dd.equ.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.dd.equ.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative -tolerance 1.0e-8 bjt2d.dd.equ.exo bjt2d.dd.equ.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # First-order current caclulation using FEM
  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.dd.current.yaml --current --current-lo
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.dd.current.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.dd.current.exo
         bjt2d.dd.current.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Second-order current caclulation using FEM
  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.dd.hocurrent.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.dd.hocurrent.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.dd.hocurrent.exo
         bjt2d.dd.hocurrent.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # First-order current caclulation using SG
  TEST_13 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.sg.current.yaml --current --current-lo
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the consistency of the exodus output
  TEST_14 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.sg.current.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_15 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.sg.current.exo
         bjt2d.sg.current.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Second-order current caclulation using SG
  TEST_16 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.sg.hocurrent.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the consistency of the exodus output
  TEST_17 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.sg.hocurrent.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_18 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.sg.hocurrent.exo
         bjt2d.sg.hocurrent.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_bjt2d PROPERTY LABELS bjt2d nightly debugexclude)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_bjt2d.dd.current_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    bjt2d.dd.current.inp
    bjt2d.dd.current.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i bjt2d.dd.current.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS bjt2d.dd.current.inp.xml bjt2d.dd.current.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_bjt2d.dd.current_input PROPERTY LABELS
    bjt2d.dd.current_input inputVerification)
ENDIF()
