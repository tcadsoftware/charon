INCLUDE(TribitsAddAdvancedTest)

####################################################
# Check coordinate scaling
SET(testName charon_mp.dof-scaling)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    dof-scale.sg.dd.inp
    dof-scale.exodiff
    mesh.exo.4.0
    mesh.exo.4.1
    mesh.exo.4.2
    mesh.exo.4.3
    dof-scale.sg.dd.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i dof-scale.sg.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=dof-scale.sg.dd.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "dof-scale.sg.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f dof-scale.exodiff dof-scale.sg.dd.exo dof-scale.sg.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.dof-scaling PROPERTY LABELS dof-scaling
    nightly)
ENDIF()
