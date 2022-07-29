

####################################################
# NPN 1D BJT Dakota test
SET(testName charon_mp_2n2222.sg_dakota)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  ENVIRONMENT APREPRO_PATH=${SEACAS_BINARY_DIR}/applications/aprepro
  ENVIRONMENT CHARON_EXECUTABLE_PATH=${PACKAGE_BINARY_DIR}

  TEST_0 COPY_FILES_TO_TEST_DIR
    template_dir.tar
    dakota.in
    lhs.gold.dat

  TEST_1 CMND tar ARGS -xvf template_dir.tar

  TEST_2 CMND dakota
    ARGS  -i dakota.in
    NUM_MPI_PROCS 4
    
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/compareDakotaOutput.py
        lhs.dat
	lhs.gold.dat
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "files are the same"

    ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.sg_dakota PROPERTY LABELS 2n2222.sg_dakota
    dakota)
ENDIF()
