INCLUDE(TribitsAddAdvancedTest)

####################################################
# Test the FEM NLP equation
SET(testName charon_mp_bjt2d)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    bjt2d.nlp.inp
    bjt2d.nlp.gold.exo
    bjt2d.gen
    bjt2d.gen.4.0
    bjt2d.gen.4.1
    bjt2d.gen.4.2
    bjt2d.gen.4.3
    bjt2d.dd.equ.inp
    bjt2d.dd.equ.gold.exo
    bjt2d.dd.inp
    bjt2d.dd.gold.exo
    bjt2d.sg.inp
    bjt2d.sg.gold.exo

  ########################################################
  # Nonlinear Poisson simulation
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i bjt2d.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative -tolerance 1.0e-8 bjt2d.nlp.exo bjt2d.nlp.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Equilibrium FEM simulation
  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i bjt2d.dd.equ.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.dd.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.dd.equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative -tolerance 1.0e-8 bjt2d.dd.equ.exo bjt2d.dd.equ.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Second-order current caclulation using FEM
  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i bjt2d.dd.inp

  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.dd.exo
         bjt2d.dd.gold.exo
    NUM_MPI_PROCS 1

  ########################################################
  # Second-order current caclulation using SG
  TEST_13 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i bjt2d.sg.inp

  TEST_14 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=bjt2d.sg.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the consistency of the exodus output
  TEST_15 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "bjt2d.sg.exo.4."
    NUM_MPI_PROCS 1

  TEST_16 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1e-5
         -Floor 1e-12
         bjt2d.sg.exo
         bjt2d.sg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_bjt2d PROPERTY LABELS bjt2d nightly debugexclude)
ENDIF()
