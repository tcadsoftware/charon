INCLUDE(TribitsAddAdvancedTest)

###############################################################################
# test the incomplete ionization model for DD set in MB regime
SET(testName charon_mp_resistor_MBFull.dd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.dd.inp
    resistor_MBFull.dd.gold.exo
    resistor_MBFull.dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MBFull.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.dd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.dd.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.dd.inp
    resistor_MB.dd.gold.exo
    resistor_MB.dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MB.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.dd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.dd.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.dd.inp
    resistor_FDFull.dd.gold.exo
    resistor_FDFull.dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FDFull.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.dd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.dd.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.dd.inp
    resistor_FD.dd.gold.exo
    resistor_FD.dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FD.dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.dd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FD.dd.output
         --anodeCurrent   -5.42799948
         --cathodeCurrent 5.42799948
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.ddlattice.inp
    resistor_MBFull.ddlattice.gold.exo
    resistor_MBFull.ddlattice.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MBFull.ddlattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.ddlattice.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.ddlattice.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.ddlattice.output
         --anodeCurrent   -5.60722256
         --cathodeCurrent 5.60722256
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.ddlattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.ddlattice.inp
    resistor_MB.ddlattice.gold.exo
    resistor_MB.ddlattice.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MB.ddlattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.ddlattice.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.ddlattice.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename resistor_MB.ddlattice.output
         --anodeCurrent -5.43381304
         --cathodeCurrent 5.43381304
         --tolerance 1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.ddlattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.ddlattice.inp
    resistor_FDFull.ddlattice.gold.exo
    resistor_FDFull.ddlattice.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FDFull.ddlattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.ddlattice.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.ddlattice.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.ddlattice.output
         --anodeCurrent   -5.60722256
         --cathodeCurrent 5.60722256
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.ddlattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.ddlattice.inp
    resistor_FD.ddlattice.gold.exo
    resistor_FD.ddlattice.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FD.ddlattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.ddlattice.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.ddlattice.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename resistor_FD.ddlattice.output
         --anodeCurrent -5.43380637
         --cathodeCurrent 5.43380637
         --tolerance 1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.ddlattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.SGCVFEMDD.inp
    resistor_MBFull.SGCVFEMDD.gold.exo
    resistor_MBFull.SGCVFEMDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MBFull.SGCVFEMDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.SGCVFEMDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.SGCVFEMDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.SGCVFEMDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.SGCVFEMDD.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.SGCVFEMDD.inp
    resistor_MB.SGCVFEMDD.gold.exo
    resistor_MB.SGCVFEMDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MB.SGCVFEMDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.SGCVFEMDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.SGCVFEMDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.SGCVFEMDD.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.SGCVFEMDD.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FDFull.SGCVFEMDD.inp
    resistor_FDFull.SGCVFEMDD.gold.exo
    resistor_FDFull.SGCVFEMDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FDFull.SGCVFEMDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FDFull.SGCVFEMDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FDFull.SGCVFEMDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FDFull.SGCVFEMDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FDFull.SGCVFEMDD.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_FD.SGCVFEMDD.inp
    resistor_FD.SGCVFEMDD.gold.exo
    resistor_FD.SGCVFEMDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_FD.SGCVFEMDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_FD.SGCVFEMDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_FD.SGCVFEMDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_FD.SGCVFEMDD.output
         --anodeCurrent   -5.42799953
         --cathodeCurrent 5.42799946
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_FD.SGCVFEMDD.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MBFull.EFFPGDD.inp
    resistor_MBFull.EFFPGDD.gold.exo
    resistor_MBFull.EFFPGDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MBFull.EFFPGDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MBFull.EFFPGDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MBFull.EFFPGDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MBFull.EFFPGDD.output
         --anodeCurrent   -5.6076177
         --cathodeCurrent 5.6076177
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MBFull.EFFPGDD.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    resistor.exo.4.0
    resistor.exo.4.1
    resistor.exo.4.2
    resistor.exo.4.3
    resistor_MB.EFFPGDD.inp
    resistor_MB.EFFPGDD.gold.exo
    resistor_MB.EFFPGDD.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i resistor_MB.EFFPGDD.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=resistor_MB.EFFPGDD.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE resistor_MB.EFFPGDD.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       resistor_MB.EFFPGDD.output
         --anodeCurrent   -5.42802007
         --cathodeCurrent 5.42802007
         --tolerance      1e-7
    PASS_REGULAR_EXPRESSION "Test Passed"

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "resistor_MB.EFFPGDD.exo.4."
    NUM_MPI_PROCS 1

    TEST_5
    EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
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
