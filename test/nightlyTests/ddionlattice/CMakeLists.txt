INCLUDE(TribitsAddAdvancedTest)

################################################################################
# Test the FEM-SUPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.supg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    memristor.nlp.inp
    memristor.dd.equ.inp
    memristor.ddionlatt.trans.supg.inp
    memristor.exo.4.0
    memristor.exo.4.1
    memristor.exo.4.2
    memristor.exo.4.3
    memristor.ddionlatt.exodiff
    memristor.ddionlatt.trans.supg.gold.exo

  # Remove old csv file
  TEST_1 CMND rm ARGS -f current_time_ddionlatt_supg.csv
    PASS_ANY

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.nlp.inp

  # Nonlinear Poisson simulation
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.nlp.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.dd.equ.inp

  # Equilibrium DD simulation
  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.dd.equ.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddionlatt.trans.supg.inp

  # Transient DDIonLattice simulation
  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.ddionlatt.trans.supg.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "memristor.ddionlatt.trans.supg.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         memristor.ddionlatt.exodiff
         memristor.ddionlatt.trans.supg.exo
         memristor.ddionlatt.trans.supg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.supg  PROPERTY LABELS tiox.ddionlatt.supg  nightly)
ENDIF()


################################################################################
# Test the FEM-EFFPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.effpg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    memristor.nlp.inp
    memristor.dd.equ.inp
    memristor.ddionlatt.trans.effpg.inp
    memristor.exo.4.0
    memristor.exo.4.1
    memristor.exo.4.2
    memristor.exo.4.3
    memristor.ddionlatt.exodiff
    memristor.ddionlatt.trans.effpg.gold.exo

  # Remove old csv file
  TEST_1 CMND rm ARGS -f current_time_ddionlatt_effpg.csv
    PASS_ANY

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.nlp.inp

  # Nonlinear Poisson simulation
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.nlp.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.dd.equ.inp

  # Equilibrium DD simulation
  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.dd.equ.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddionlatt.trans.effpg.inp

  # Transient DDIonLattice simulation
  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.ddionlatt.trans.effpg.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "memristor.ddionlatt.trans.effpg.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         memristor.ddionlatt.exodiff
         memristor.ddionlatt.trans.effpg.exo
         memristor.ddionlatt.trans.effpg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.effpg  PROPERTY LABELS
    tiox.ddionlatt.effpg  nightly debugexclude)
ENDIF()


################################################################################
# Test the FEM-SymEFFPG DDIonLattice formulation for a brick tiox device
SET(testName charon_mp_tiox.ddionlatt.symeffpg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    memristor.nlp.inp
    memristor.dd.equ.inp
    memristor.ddionlatt.trans.symeffpg.inp
    memristor.exo.4.0
    memristor.exo.4.1
    memristor.exo.4.2
    memristor.exo.4.3
    memristor.ddionlatt.exodiff
    memristor.ddionlatt.trans.symeffpg.gold.exo

  # Remove old csv file
  TEST_1 CMND rm ARGS -f current_time_ddionlatt_symeffpg.csv
    PASS_ANY

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.nlp.inp

  # Nonlinear Poisson simulation
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.nlp.inp.xml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.dd.equ.inp

  # Equilibrium DD simulation
  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.dd.equ.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i memristor.ddionlatt.trans.symeffpg.inp

  # Transient DDIonLattice simulation
  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=memristor.ddionlatt.trans.symeffpg.inp.xml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "memristor.ddionlatt.trans.symeffpg.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         memristor.ddionlatt.exodiff
         memristor.ddionlatt.trans.symeffpg.exo
         memristor.ddionlatt.trans.symeffpg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_tiox.ddionlatt.symeffpg  PROPERTY LABELS
    tiox.ddionlatt.symeffpg  nightly)
ENDIF()

