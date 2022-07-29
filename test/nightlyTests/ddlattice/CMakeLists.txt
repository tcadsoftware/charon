INCLUDE(TribitsAddAdvancedTest)

################################################################################
# Test the FEM NLP equation
SET(testName charon_mp_pndiode.nlp.ddlatice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.yaml
    pndiode.nlp.gold.exo
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.yaml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 1.0e-8
         -Floor 1.0e-13
         pndiode.nlp.exo
         pndiode.nlp.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.ddlatice PROPERTY LABELS  pndiode.nlp
    nightly)
ENDIF()

################################################################################
# Test the FEM Lattice equation
SET(testName charon_mp_pndiode.lattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.lattice.inp
    pndiode.lattice.gold.exo
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.lattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.lattice.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.lattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative -tolerance 1e-8 pndiode.lattice.exo pndiode.lattice.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.lattice PROPERTY LABELS
    pndiode.lattice nightly)
ENDIF()

################################################################################
# Test the SUPG-FEM DDLattice formulation with thermal Dirichlet BC
SET(testName charon_mp_pndiode.ddlattice)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.ddlattice.inp
    pndiode.ddlattice.gold.exo
    pndiode.ddlattice.equ.exo.4.0
    pndiode.ddlattice.equ.exo.4.1
    pndiode.ddlattice.equ.exo.4.2
    pndiode.ddlattice.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.ddlattice.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.ddlattice.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.ddlattice.output

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.ddlattice.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -relative
         -tolerance 2e-8
         pndiode.ddlattice.exo
         pndiode.ddlattice.gold.exo
    NUM_MPI_PROCS 1

  # now test that the current is correct
  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename pndiode.ddlattice.output
         --anodeCurrent -3.85558716e-02
         --cathodeCurrent 3.85558716e-02
         --tolerance 1e-6
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.ddlattice  PROPERTY LABELS
    pndiode.ddlattice  nightly)
ENDIF()

################################################################################
# Test the SUPG-FEM DDLattice formulation with constant thermal Neumann BC
SET(testName charon_mp_pndiode.ddlattice.cnbc)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.ddlattice.cnbc.inp
    pndiode.ddlattice.equ.exo.4.0
    pndiode.ddlattice.equ.exo.4.1
    pndiode.ddlattice.equ.exo.4.2
    pndiode.ddlattice.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.ddlattice.cnbc.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.ddlattice.cnbc.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.ddlattice.cnbc.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename pndiode.ddlattice.cnbc.output
         --anodeCurrent -3.48704528e-02
         --cathodeCurrent 3.48704528e-02
         --tolerance 1e-6
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.ddlattice.cnbc PROPERTY LABELS
    pndiode.ddlattice.cnbc  nightly)
ENDIF()

################################################################################
# Test the SUPG-FEM DDLattice formulation with Fermi-Dirac and Band Gap Narrowing
SET(testName charon_mp_pndiode.ddlattice.fd.bgn)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.ddlattice.fd.bgn.inp
    pndiode.fd.bgn.equ.exo
    pndiode.fd.bgn.equ.exo.4.0
    pndiode.fd.bgn.equ.exo.4.1
    pndiode.fd.bgn.equ.exo.4.2
    pndiode.fd.bgn.equ.exo.4.3

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.ddlattice.fd.bgn.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.ddlattice.fd.bgn.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.ddlattice.fd.bgn.output

  # test that the current is correct
  TEST_3 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename pndiode.ddlattice.fd.bgn.output
         --anodeCurrent -4.58731423e-02
         --cathodeCurrent 4.58731423e-02
         --tolerance 1e-6
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.ddlattice.fd.bgn  PROPERTY LABELS
    pndiode.ddlattice.fd.bgn  nightly)
ENDIF()
