INCLUDE(TribitsAddAdvancedTest)

####################################################
#  resistor 1D HB DD large-signal 
####################################################

SET(testName charon_mp.hb-dd.resistor.ls)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor.hb-dd.ls.yaml
    gold_resistor.hb-dd.ls.exo
    tolerance.exodiff

  # run hb-dd large signal problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.hb-dd.ls.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # epu auto the result to prepare for exodiff
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor.hb-dd.ls.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f tolerance.exodiff
         resistor.hb-dd.ls.exo
         gold_resistor.hb-dd.ls.exo
    NUM_MPI_PROCS 1

  # check anode_silicon_CosH0.000000__Current == anode_silicon_SinH1.000000__Current 
  #TEST_4 EXEC 

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.resistor.ls  PROPERTY LABELS
    hb-dd.resistor.ls nightly)
ENDIF()



####################################################
#  resistor 1D HB DD small-signal 
####################################################

SET(testName charon_mp.hb-dd.resistor.ss)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor.hb-dd.ss.yaml
    gold_resistor.hb-dd.ss.exo
    tolerance.exodiff

  # run hb-dd large signal problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.hb-dd.ss.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # epu auto the result to prepare for exodiff
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor.hb-dd.ss.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f tolerance.exodiff
         resistor.hb-dd.ss.exo
         gold_resistor.hb-dd.ss.exo
    NUM_MPI_PROCS 1

  # check anode_silicon_CosH0.000000__Current == anode_silicon_SinH1.000000__Current 
  #TEST_4 EXEC 

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.resistor.ss  PROPERTY LABELS
    hb-dd.resistor.ss nightly)
ENDIF()
