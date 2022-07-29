INCLUDE(TribitsAddAdvancedTest)

####################################################
#  resistor 1D HB laplace 
####################################################

SET(testName charon_mp.hb-laplace.resistor.ls)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor.hb-laplace.ls.inp
    gold_resistor.hb-laplace.ls.exo
    tol-laplace.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor.hb-laplace.ls.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.hb-laplace.ls.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor.hb-laplace.ls.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor.hb-laplace.ls.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f tol-laplace.exodiff
         resistor.hb-laplace.ls.exo
         gold_resistor.hb-laplace.ls.exo
    NUM_MPI_PROCS 1

  # check anode_silicon_CosH0.000000__Current == anode_silicon_SinH1.000000__Current 
  #TEST_5 EXEC 

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-laplace.resistor.ls  PROPERTY LABELS
    hb-laplace.resistor.ls nightly)
ENDIF()



####################################################
#  resistor 1D HB DD 
####################################################

SET(testName charon_mp.hb-dd.resistor.ls)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor.hb-laplace.ls.inp
    resistor.hb-dd.ls.inp
    gold_resistor.hb-dd.ls.exo
    tol-dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor.hb-laplace.ls.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.hb-laplace.ls.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor.hb-laplace.ls.inp.xml.out

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor.hb-dd.ls.inp

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=resistor.hb-dd.ls.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor.hb-dd.ls.inp.xml.out

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor.hb-dd.ls.exo.4."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f tol-dd.exodiff
         resistor.hb-dd.ls.exo
         gold_resistor.hb-dd.ls.exo
    NUM_MPI_PROCS 1

  # check anode_silicon_CosH0.000000__Current == anode_silicon_SinH1.000000__Current 
  #TEST_5 EXEC 

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.resistor.ls  PROPERTY LABELS
    hb-dd.resistor.ls nightly)
ENDIF()



