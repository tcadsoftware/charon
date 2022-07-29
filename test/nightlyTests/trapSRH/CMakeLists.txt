INCLUDE(TribitsAddAdvancedTest)
SET(TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS 30)

####################################################
# TrapSRH in bulk, level type
SET(testName charon_mp_trapSRH_level_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.jou
    si-diode.nlp.inp
    si-diode.dd.inp
    si-diode.dd.exo.4.0
    si-diode.dd.exo.4.1
    si-diode.dd.exo.4.2
    si-diode.dd.exo.4.3
    si-diode.dd.forward.inp
    si-diode.dd.forward.exodiff
    si-diode.dd.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.forward.exodiff
         si-diode.dd.forward.exo
	 si-diode.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_level_bulk PROPERTY LABELS
      trapSRH_level_bulk nightly debugexclude)
ENDIF()

####################################################
# TrapSRH in bulk, uniform type
SET(testName charon_mp_trapSRH_uni_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.uni.inp
    si-diode.dd.uni.exo.4.0
    si-diode.dd.uni.exo.4.1
    si-diode.dd.uni.exo.4.2
    si-diode.dd.uni.exo.4.3
    si-diode.dd.uni.forward.inp
    si-diode.dd.uni.forward.exodiff
    si-diode.dd.uni.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.uni.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.uni.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.uni.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.uni.forward.exodiff
         si-diode.dd.uni.forward.exo
	 si-diode.dd.uni.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_uni_bulk PROPERTY LABELS
      trapSRH_uni_bulk nightly debugexclude)
ENDIF()

####################################################
# TrapSRH in bulk, exponential type
SET(testName charon_mp_trapSRH_exp_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.exp.inp
    si-diode.dd.exp.exo.4.0
    si-diode.dd.exp.exo.4.1
    si-diode.dd.exp.exo.4.2
    si-diode.dd.exp.exo.4.3
    si-diode.dd.exp.forward.inp
    si-diode.dd.exp.forward.exodiff
    si-diode.dd.exp.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.exp.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.exp.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.exp.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.exp.forward.exodiff
         si-diode.dd.exp.forward.exo
	 si-diode.dd.exp.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_exp_bulk PROPERTY LABELS
      trapSRH_exp_bulk nightly debugexclude)
ENDIF()

####################################################
# TrapSRH in bulk, gaussian type
SET(testName charon_mp_trapSRH_gau_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.gau.inp
    si-diode.dd.gau.exo.4.0
    si-diode.dd.gau.exo.4.1
    si-diode.dd.gau.exo.4.2
    si-diode.dd.gau.exo.4.3
    si-diode.dd.gau.forward.inp
    si-diode.dd.gau.forward.exodiff
    si-diode.dd.gau.forward.gold.exo

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.gau.forward.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.gau.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.gau.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.gau.forward.exodiff
         si-diode.dd.gau.forward.exo
	 si-diode.dd.gau.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_gau_bulk PROPERTY LABELS
      trapSRH_gau_bulk nightly debugexclude)
ENDIF()

####################################################
# TrapSRH on interface, level type
SET(testName charon_mp_trapSRH_level_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet_msh.jou
    nmosfet.nlp.inp
    nmosfet.dd.inp
    nmosfet.dd.0.1.inp
    nmosfet.dd.0.1.exo.4.0
    nmosfet.dd.0.1.exo.4.1
    nmosfet.dd.0.1.exo.4.2
    nmosfet.dd.0.1.exo.4.3
    nmosfet.dd.gate_sweep.inp
    nmosfet.dd.gate_sweep.gold.exo
    nmosfet.dd.gate_sweep.exodiff 
    
  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.gate_sweep.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.gate_sweep.exodiff
         nmosfet.dd.gate_sweep.exo
	 nmosfet.dd.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_level_inter PROPERTY LABELS
      trapSRH_level_inter nightly debugexclude)
ENDIF()

####################################################
# TrapSRH on interface, uniform type
SET(testName charon_mp_trapSRH_uni_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.dd.uni.inp
    nmosfet.dd.uni.0.1.inp
    nmosfet.dd.uni.0.1.exo.4.0
    nmosfet.dd.uni.0.1.exo.4.1
    nmosfet.dd.uni.0.1.exo.4.2
    nmosfet.dd.uni.0.1.exo.4.3
    nmosfet.dd.uni.gate_sweep.inp
    nmosfet.dd.uni.gate_sweep.gold.exo
    nmosfet.dd.uni.gate_sweep.exodiff 
    
  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.uni.gate_sweep.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.uni.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.uni.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.uni.gate_sweep.exodiff
         nmosfet.dd.uni.gate_sweep.exo
	 nmosfet.dd.uni.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_uni_inter PROPERTY LABELS
      trapSRH_uni_inter nightly debugexclude)
ENDIF()

####################################################
# TrapSRH on interface, exponential type
SET(testName charon_mp_trapSRH_exp_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.dd.exp.inp
    nmosfet.dd.exp.0.1.inp
    nmosfet.dd.exp.0.1.exo.4.0
    nmosfet.dd.exp.0.1.exo.4.1
    nmosfet.dd.exp.0.1.exo.4.2
    nmosfet.dd.exp.0.1.exo.4.3
    nmosfet.dd.exp.gate_sweep.inp
    nmosfet.dd.exp.gate_sweep.gold.exo
    nmosfet.dd.exp.gate_sweep.exodiff 
    
  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.exp.gate_sweep.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.exp.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.exp.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.exp.gate_sweep.exodiff
         nmosfet.dd.exp.gate_sweep.exo
	 nmosfet.dd.exp.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_exp_inter PROPERTY LABELS
      trapSRH_exp_inter nightly debugexclude)
ENDIF()

####################################################
# TrapSRH on interface, gaussian type
SET(testName charon_mp_trapSRH_gau_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosfet.dd.gau.inp
    nmosfet.dd.gau.0.1.inp
    nmosfet.dd.gau.0.1.exo.4.0
    nmosfet.dd.gau.0.1.exo.4.1
    nmosfet.dd.gau.0.1.exo.4.2
    nmosfet.dd.gau.0.1.exo.4.3
    nmosfet.dd.gau.gate_sweep.inp
    nmosfet.dd.gau.gate_sweep.gold.exo
    nmosfet.dd.gau.gate_sweep.exodiff 
    
  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
      
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.gau.gate_sweep.inp
    
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.gau.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.gau.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.gau.gate_sweep.exodiff
         nmosfet.dd.gau.gate_sweep.exo
	 nmosfet.dd.gau.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_trapSRH_gau_inter PROPERTY LABELS
      trapSRH_gau_inter nightly debugexclude)
ENDIF()
