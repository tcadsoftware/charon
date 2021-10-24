INCLUDE(TribitsAddAdvancedTest)

###############################################################################
# test the incomplete ionization model for DD set in MB regime
SET(testName charon_mp_resistor_MBFull.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.dd.yaml
    resistor_MBFull.dd.gold.exo
    resistor_MBFull.dd.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.dd.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.dd.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MBFull.dd.exodiff
         resistor_MBFull.dd.exo
         resistor_MBFull.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.dd PROPERTY LABELS
    resistor_MBFull.dd nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DD set in MB regime
SET(testName charon_mp_resistor_MB.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.dd.yaml
    resistor_MB.dd.gold.exo
    resistor_MB.dd.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.dd.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.dd.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MB.dd.exodiff
         resistor_MB.dd.exo
         resistor_MB.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.dd PROPERTY LABELS
    resistor_MB.dd nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DD set in FD regime
SET(testName charon_mp_resistor_FDFull.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.dd.yaml
    resistor_FDFull.dd.gold.exo
    resistor_FDFull.dd.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.dd.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.dd.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FDFull.dd.exodiff
         resistor_FDFull.dd.exo
         resistor_FDFull.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.dd PROPERTY LABELS
    resistor_FDFull.dd nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DD set in FD regime
SET(testName charon_mp_resistor_FD.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.dd.yaml
    resistor_FD.dd.gold.exo
    resistor_FD.dd.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.dd.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FD.dd.output
         --anodeCurrent   -5.42799948
         --cathodeCurrent 5.42799948
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.dd.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FD.dd.exodiff
         resistor_FD.dd.exo
         resistor_FD.dd.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.dd PROPERTY LABELS
    resistor_FD.dd nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DDLattice set in MB regime
SET(testName charon_mp_resistor_MBFull.ddlattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.ddlattice.yaml
    resistor_MBFull.ddlattice.gold.exo
    resistor_MBFull.ddlattice.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.ddlattice.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.ddlattice.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.ddlattice.output
         --anodeCurrent   -5.60722256
         --cathodeCurrent 5.60722256
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.ddlattice.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MBFull.ddlattice.exodiff
         resistor_MBFull.ddlattice.exo
         resistor_MBFull.ddlattice.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.ddlattice PROPERTY LABELS
    resistor_MBFull.ddlattice nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DDLattice set in MB regime
SET(testName charon_mp_resistor_MB.ddlattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.ddlattice.yaml
    resistor_MB.ddlattice.gold.exo
    resistor_MB.ddlattice.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.ddlattice.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.ddlattice.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename resistor_MB.ddlattice.output
         --anodeCurrent -5.43381304
         --cathodeCurrent 5.43381304
         --tolerance 1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.ddlattice.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MB.ddlattice.exodiff
         resistor_MB.ddlattice.exo
         resistor_MB.ddlattice.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.ddlattice PROPERTY LABELS
    resistor_MB.ddlattice nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DDLattice set in FD regime
SET(testName charon_mp_resistor_FDFull.ddlattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.ddlattice.yaml
    resistor_FDFull.ddlattice.gold.exo
    resistor_FDFull.ddlattice.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.ddlattice.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.ddlattice.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.ddlattice.output
         --anodeCurrent   -5.60722256
         --cathodeCurrent 5.60722256
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.ddlattice.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FDFull.ddlattice.exodiff
         resistor_FDFull.ddlattice.exo
         resistor_FDFull.ddlattice.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.ddlattice PROPERTY LABELS
    resistor_FDFull.ddlattice nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for DDLattice set in FD regime
SET(testName charon_mp_resistor_FD.ddlattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.ddlattice.yaml
    resistor_FD.ddlattice.gold.exo
    resistor_FD.ddlattice.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.ddlattice.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.ddlattice.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename resistor_FD.ddlattice.output
         --anodeCurrent -5.43380637
         --cathodeCurrent 5.43380637
         --tolerance 1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.ddlattice.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FD.ddlattice.exodiff
         resistor_FD.ddlattice.exo
         resistor_FD.ddlattice.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.ddlattice PROPERTY LABELS
    resistor_FD.ddlattice nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for SGCVFEMDD set in MB regime
SET(testName charon_mp_resistor_MBFull.SGCVFEMDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.SGCVFEMDD.yaml
    resistor_MBFull.SGCVFEMDD.gold.exo
    resistor_MBFull.SGCVFEMDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.SGCVFEMDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.SGCVFEMDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.SGCVFEMDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.SGCVFEMDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MBFull.SGCVFEMDD.exodiff
         resistor_MBFull.SGCVFEMDD.exo
         resistor_MBFull.SGCVFEMDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.SGCVFEMDD PROPERTY LABELS
    resistor_MBFull.SGCVFEMDD nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for SGCVFEMDD set in MB regime
SET(testName charon_mp_resistor_MB.SGCVFEMDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.SGCVFEMDD.yaml
    resistor_MB.SGCVFEMDD.gold.exo
    resistor_MB.SGCVFEMDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.SGCVFEMDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.SGCVFEMDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.SGCVFEMDD.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.SGCVFEMDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MB.SGCVFEMDD.exodiff
         resistor_MB.SGCVFEMDD.exo
         resistor_MB.SGCVFEMDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.SGCVFEMDD PROPERTY LABELS
    resistor_MB.SGCVFEMDD nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for SGCVFEMDD set in FD regime
SET(testName charon_mp_resistor_FDFull.SGCVFEMDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.SGCVFEMDD.yaml
    resistor_FDFull.SGCVFEMDD.gold.exo
    resistor_FDFull.SGCVFEMDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.SGCVFEMDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.SGCVFEMDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.SGCVFEMDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.SGCVFEMDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FDFull.SGCVFEMDD.exodiff
         resistor_FDFull.SGCVFEMDD.exo
         resistor_FDFull.SGCVFEMDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.SGCVFEMDD PROPERTY LABELS
    resistor_FDFull.SGCVFEMDD nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for SGCVFEMDD set in FD regime
SET(testName charon_mp_resistor_FD.SGCVFEMDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.SGCVFEMDD.yaml
    resistor_FD.SGCVFEMDD.gold.exo
    resistor_FD.SGCVFEMDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.SGCVFEMDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.SGCVFEMDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FD.SGCVFEMDD.output
         --anodeCurrent   -5.42799953
         --cathodeCurrent 5.42799946
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.SGCVFEMDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_FD.SGCVFEMDD.exodiff
         resistor_FD.SGCVFEMDD.exo
         resistor_FD.SGCVFEMDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.SGCVFEMDD PROPERTY LABELS
    resistor_FD.SGCVFEMDD nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for EFFPGDD set in MB regime
SET(testName charon_mp_resistor_MBFull.EFFPGDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.EFFPGDD.yaml
    resistor_MBFull.EFFPGDD.gold.exo
    resistor_MBFull.EFFPGDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.EFFPGDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.EFFPGDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.EFFPGDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.EFFPGDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MBFull.EFFPGDD.exodiff
         resistor_MBFull.EFFPGDD.exo
         resistor_MBFull.EFFPGDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.EFFPGDD PROPERTY LABELS
    resistor_MBFull.EFFPGDD nightly)
ENDIF()


###############################################################################
# test the incomplete ionization model for EFFPGDD set in MB regime
SET(testName charon_mp_resistor_MB.EFFPGDD)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.EFFPGDD.yaml
    resistor_MB.EFFPGDD.gold.exo
    resistor_MB.EFFPGDD.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.EFFPGDD.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.EFFPGDD.output

  # test that the current is correct
  TEST_2 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.EFFPGDD.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.EFFPGDD.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         resistor_MB.EFFPGDD.exodiff
         resistor_MB.EFFPGDD.exo
         resistor_MB.EFFPGDD.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.EFFPGDD PROPERTY LABELS
    resistor_MB.EFFPGDD nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FD.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FD.dd.inp
    resistor_FD.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FD.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FD.dd.inp.xml resistor_FD.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.dd_input PROPERTY LABELS
    resistor_FD.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FD.ddlattice_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FD.ddlattice.inp
    resistor_FD.ddlattice.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FD.ddlattice.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FD.ddlattice.inp.xml resistor_FD.ddlattice.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.ddlattice_input PROPERTY LABELS
    resistor_FD.ddlattice_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FDFull.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FDFull.dd.inp
    resistor_FDFull.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FDFull.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FDFull.dd.inp.xml resistor_FDFull.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.dd_input PROPERTY LABELS
    resistor_FDFull.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FDFull.ddlattice_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FDFull.ddlattice.inp
    resistor_FDFull.ddlattice.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FDFull.ddlattice.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FDFull.ddlattice.inp.xml resistor_FDFull.ddlattice.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.ddlattice_input PROPERTY LABELS
    resistor_FDFull.ddlattice_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FDFull.SGCVFEMDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FDFull.SGCVFEMDD.inp
    resistor_FDFull.SGCVFEMDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FDFull.SGCVFEMDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FDFull.SGCVFEMDD.inp.xml resistor_FDFull.SGCVFEMDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FDFull.SGCVFEMDD_input PROPERTY LABELS
    resistor_FDFull.SGCVFEMDD_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_FD.SGCVFEMDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_FD.SGCVFEMDD.inp
    resistor_FD.SGCVFEMDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_FD.SGCVFEMDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_FD.SGCVFEMDD.inp.xml resistor_FD.SGCVFEMDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_FD.SGCVFEMDD_input PROPERTY LABELS
    resistor_FD.SGCVFEMDD_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MB.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MB.dd.inp
    resistor_MB.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MB.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MB.dd.inp.xml resistor_MB.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.dd_input PROPERTY LABELS
    resistor_MB.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MB.ddlattice_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MB.ddlattice.inp
    resistor_MB.ddlattice.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MB.ddlattice.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MB.ddlattice.inp.xml resistor_MB.ddlattice.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.ddlattice_input PROPERTY LABELS
    resistor_MB.ddlattice_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MBFull.dd_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MBFull.dd.inp
    resistor_MBFull.dd.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MBFull.dd.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MBFull.dd.inp.xml resistor_MBFull.dd.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.dd_input PROPERTY LABELS
    resistor_MBFull.dd_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MBFull.ddlattice_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MBFull.ddlattice.inp
    resistor_MBFull.ddlattice.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MBFull.ddlattice.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MBFull.ddlattice.inp.xml resistor_MBFull.ddlattice.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.ddlattice_input PROPERTY LABELS
    resistor_MBFull.ddlattice_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MBFull.SGCVFEMDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MBFull.SGCVFEMDD.inp
    resistor_MBFull.SGCVFEMDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MBFull.SGCVFEMDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MBFull.SGCVFEMDD.inp.xml resistor_MBFull.SGCVFEMDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.SGCVFEMDD_input PROPERTY LABELS
    resistor_MBFull.SGCVFEMDD_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MB.SGCVFEMDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MB.SGCVFEMDD.inp
    resistor_MB.SGCVFEMDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MB.SGCVFEMDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MB.SGCVFEMDD.inp.xml resistor_MB.SGCVFEMDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.SGCVFEMDD_input PROPERTY LABELS
    resistor_MB.SGCVFEMDD_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MB.EFFPGDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MB.EFFPGDD.inp
    resistor_MB.EFFPGDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MB.EFFPGDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MB.EFFPGDD.inp.xml resistor_MB.EFFPGDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MB.EFFPGDD_input PROPERTY LABELS
    resistor_MB.EFFPGDD_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_resistor_MBFull.EFFPGDD_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor_MBFull.EFFPGDD.inp
    resistor_MBFull.EFFPGDD.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i resistor_MBFull.EFFPGDD.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS resistor_MBFull.EFFPGDD.inp.xml resistor_MBFull.EFFPGDD.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_resistor_MBFull.EFFPGDD_input PROPERTY LABELS
    resistor_MBFull.EFFPGDD_input inputVerification)
ENDIF()
