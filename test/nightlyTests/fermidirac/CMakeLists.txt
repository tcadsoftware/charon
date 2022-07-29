INCLUDE(TribitsAddAdvancedTest)

###############################################################################
# test SUPG-FEM and Fermi Dirac and FD Formula = Schroeder (default)
SET(testName charon_mp_pndiode.dd.fd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.fd.inp
    pndiode.equ.exo.4.0
    pndiode.equ.exo.4.1
    pndiode.equ.exo.4.2
    pndiode.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.fd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.fd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.fd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.fd.output
         --anodeCurrent   -1.23406365e-04
         --cathodeCurrent 1.23406221e-04
         --tolerance      8.0e-07
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.fd  PROPERTY LABELS
    pndiode.dd.fd nightly)
ENDIF()

###############################################################################
# test CVFEM-SG and Fermi Dirac and FD Formula = Schroeder (default)
SET(testName charon_mp_pndiode.sg.fd)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.fd.inp
    pndiode.equ.exo.4.0
    pndiode.equ.exo.4.1
    pndiode.equ.exo.4.2
    pndiode.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.fd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.fd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.fd.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.fd.output
         --anodeCurrent   -1.22998203e-04
         --cathodeCurrent 1.22998122e-04
         --tolerance      2.0e-06
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.fd  PROPERTY LABELS
    pndiode.sg.fd nightly)
ENDIF()

###############################################################################
# test SUPG-FEM and Fermi Dirac and FD Formula = Diffusion
SET(testName charon_mp_pndiode.dd.fd.diff)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.fd.diff.inp
    pndiode.equ.exo.4.0
    pndiode.equ.exo.4.1
    pndiode.equ.exo.4.2
    pndiode.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.fd.diff.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.fd.diff.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.fd.diff.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.fd.diff.output
         --anodeCurrent   -1.91271675e-04
         --cathodeCurrent 1.91272646e-04
         --tolerance      3e-07
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.fd.diff  PROPERTY LABELS
    pndiode.dd.fd.diff nightly)
ENDIF()
