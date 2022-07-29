INCLUDE(TribitsAddAdvancedTest)

####################################################
# Heterointerface diode with charge, potential heterointerface
SET(testName charon_mp_hetero_diode_charge_hi)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    hetero_diode.jou
    hetero_diode.nlp.inp
    hetero_diode.dd.inp
    hetero_diode.dd.1.inp
    hetero_diode.dd.1.exo.4.0
    hetero_diode.dd.1.exo.4.1
    hetero_diode.dd.1.exo.4.2
    hetero_diode.dd.1.exo.4.3
    hetero_diode.dd.forward.inp
    hetero_diode.dd.forward.exodiff
    hetero_diode.dd.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=hetero_diode.dd.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=hetero_diode.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "hetero_diode.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         hetero_diode.dd.forward.exodiff
         hetero_diode.dd.forward.exo
	 hetero_diode.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hetero_diode_charge_hi PROPERTY LABELS
       hetero_diode_charge_hi nightly debugexclude)
ENDIF()

####################################################
# Heterointerface diode with charge, potential heterointerface, cvfem
SET(testName charon_mp_hetero_diode_charge_hi_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    hetero_diode.sg.dd.inp
    hetero_diode.sg.dd.1.inp
    hetero_diode.sg.dd.1.exo.4.0
    hetero_diode.sg.dd.1.exo.4.1
    hetero_diode.sg.dd.1.exo.4.2
    hetero_diode.sg.dd.1.exo.4.3
    hetero_diode.sg.dd.forward.inp
    hetero_diode.sg.dd.forward.exodiff
    hetero_diode.sg.dd.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=hetero_diode.sg.dd.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=hetero_diode.sg.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "hetero_diode.sg.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         hetero_diode.sg.dd.forward.exodiff
         hetero_diode.sg.dd.forward.exo
	 hetero_diode.sg.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hetero_diode_charge_hi_sg PROPERTY LABELS
       hetero_diode_charge_hi_sg nightly debugexclude)
ENDIF()

####################################################
# Heterointerface diode with surface charge
SET(testName charon_mp_hetero_diode_sc)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    hetero_diode_sc.nlp.inp
    hetero_diode_sc.dd.inp
    hetero_diode_sc.dd.1.inp
    hetero_diode_sc.dd.1.exo.4.0
    hetero_diode_sc.dd.1.exo.4.1
    hetero_diode_sc.dd.1.exo.4.2
    hetero_diode_sc.dd.1.exo.4.3
    hetero_diode_sc.dd.forward.inp
    hetero_diode_sc.dd.forward.exodiff
    hetero_diode_sc.dd.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=hetero_diode_sc.dd.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=hetero_diode_sc.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "hetero_diode_sc.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         hetero_diode_sc.dd.forward.exodiff
         hetero_diode_sc.dd.forward.exo
	 hetero_diode_sc.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hetero_diode_sc PROPERTY LABELS
       hetero_diode_sc nightly debugexclude)
ENDIF()

####################################################
# Heterointerface diode with surface charge, sg
SET(testName charon_mp_hetero_diode_sc_sg)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    hetero_diode_sc.sg.dd.inp
    hetero_diode_sc.sg.dd.1.inp
    hetero_diode_sc.sg.dd.1.exo.4.0
    hetero_diode_sc.sg.dd.1.exo.4.1
    hetero_diode_sc.sg.dd.1.exo.4.2
    hetero_diode_sc.sg.dd.1.exo.4.3
    hetero_diode_sc.sg.dd.forward.inp
    hetero_diode_sc.sg.dd.forward.exodiff
    hetero_diode_sc.sg.dd.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=hetero_diode_sc.sg.dd.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=hetero_diode_sc.sg.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "hetero_diode_sc.sg.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         hetero_diode_sc.sg.dd.forward.exodiff
         hetero_diode_sc.sg.dd.forward.exo
	 hetero_diode_sc.sg.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hetero_diode_sc_sg PROPERTY LABELS
       hetero_diode_sc_sg nightly debugexclude)
ENDIF()
