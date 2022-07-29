INCLUDE(TribitsAddAdvancedTest)
SET(TRIBITS_ADD_ADVANCED_TEST_MAX_NUM_TEST_BLOCKS 30)

####################################################
# Dynamic Traps in bulk, level type
SET(testName charon_mp_DynamicTraps_level_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.jou
    si-diode_msh.exo.4.0
    si-diode_msh.exo.4.1
    si-diode_msh.exo.4.2
    si-diode_msh.exo.4.3
    si-diode.nlp.inp
    si-diode.nlp.exodiff
    si-diode.nlp.gold.exo
    si-diode.dd.inp
    si-diode.dd.exodiff
    si-diode.dd.gold.exo
    si-diode.dd.forward.inp
    si-diode.dd.forward.exodiff
    si-diode.dd.forward.gold.exo
     
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX 
    NOEXESUFFIX
    ARGS -f
         si-diode.nlp.exodiff
	 si-diode.nlp.exo	
    	 si-diode.nlp.gold.exo	
    NUM_MPI_PROCS 1

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.inp
    
  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.exodiff
         si-diode.dd.exo
	 si-diode.dd.gold.exo
    NUM_MPI_PROCS 1

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.forward.inp
    
  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.forward.exodiff
         si-diode.dd.forward.exo
	 si-diode.dd.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_level_bulk PROPERTY LABELS
      DynamicTraps_level_bulk nightly debugexclude)
ENDIF()

####################################################
# Dynamic Traps in bulk, continuous type
SET(testName charon_mp_DynamicTraps_cont_uni_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.uni.ini.exo.4.0
    si-diode.dd.uni.ini.exo.4.1
    si-diode.dd.uni.ini.exo.4.2
    si-diode.dd.uni.ini.exo.4.3
    si-diode.dd.uni.forward.inp
    si-diode.dd.uni.forward.exodiff
    si-diode.dd.uni.forward.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.uni.forward.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.uni.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.uni.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.uni.forward.exodiff
         si-diode.dd.uni.forward.exo
	 si-diode.dd.uni.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_cont_uni_bulk PROPERTY LABELS
      DynamicTraps_cont_uni_bulk nightly debugexclude)
ENDIF()

####################################################
# Dynamic Traps in bulk, continuous type
SET(testName charon_mp_DynamicTraps_cont_exp_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.exp.ini.exo.4.0
    si-diode.dd.exp.ini.exo.4.1
    si-diode.dd.exp.ini.exo.4.2
    si-diode.dd.exp.ini.exo.4.3
    si-diode.dd.exp.forward.inp
    si-diode.dd.exp.forward.exodiff
    si-diode.dd.exp.forward.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.exp.forward.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.exp.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.exp.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.exp.forward.exodiff
         si-diode.dd.exp.forward.exo
	 si-diode.dd.exp.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_cont_exp_bulk PROPERTY LABELS
      DynamicTraps_cont_exp_bulk nightly debugexclude)
ENDIF()

####################################################
# Dynamic Traps in bulk, continuous type
SET(testName charon_mp_DynamicTraps_cont_gau_bulk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME
    
  TEST_0 COPY_FILES_TO_TEST_DIR
    si-diode.dd.gau.ini.exo.4.0
    si-diode.dd.gau.ini.exo.4.1
    si-diode.dd.gau.ini.exo.4.2
    si-diode.dd.gau.ini.exo.4.3
    si-diode.dd.gau.forward.inp
    si-diode.dd.gau.forward.exodiff
    si-diode.dd.gau.forward.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=si-diode.dd.gau.forward.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=si-diode.dd.gau.forward.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "si-diode.dd.gau.forward.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         si-diode.dd.gau.forward.exodiff
         si-diode.dd.gau.forward.exo
	 si-diode.dd.gau.forward.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_cont_gau_bulk PROPERTY LABELS
      DynamicTraps_cont_gau_bulk nightly debugexclude)
ENDIF()

####################################################
# Dynamic Traps on interface, level 
SET(testName charon_mp_DynamicTraps_lev_inter)
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
    nmosfet.dd.lev.gate_sweep.inp
    nmosfet.dd.lev.gate_sweep.exodiff
    nmosfet.dd.lev.gate_sweep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.lev.gate_sweep.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.lev.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.lev.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.lev.gate_sweep.exodiff
         nmosfet.dd.lev.gate_sweep.exo
         nmosfet.dd.lev.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_lev_inter PROPERTY LABELS
      DynamicTraps_lev_inter nightly debugexclude)
ENDIF()

####################################################
# Dynamic Traps on interface, uniform distribution
SET(testName charon_mp_DynamicTraps_uni_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME
    
  TEST_0 COPY_FILES_TO_TEST_DIR 
    nmosfet.dd.0.1.uni.exo.4.0
    nmosfet.dd.0.1.uni.exo.4.1
    nmosfet.dd.0.1.uni.exo.4.2
    nmosfet.dd.0.1.uni.exo.4.3
    nmosfet.dd.uni.gate_sweep.inp
    nmosfet.dd.uni.gate_sweep.exodiff
    nmosfet.dd.uni.gate_sweep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.uni.gate_sweep.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.uni.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.uni.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.uni.gate_sweep.exodiff
         nmosfet.dd.uni.gate_sweep.exo
         nmosfet.dd.uni.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_uni_inter PROPERTY LABELS
      DynamicTraps_uni_inter nightly debugexclude)
ENDIF()


####################################################
# Dynamic Traps on interface, exponential distribution
SET(testName charon_mp_DynamicTraps_exp_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME
    
  TEST_0 COPY_FILES_TO_TEST_DIR 
    nmosfet.dd.0.1.exp.exo.4.0
    nmosfet.dd.0.1.exp.exo.4.1
    nmosfet.dd.0.1.exp.exo.4.2
    nmosfet.dd.0.1.exp.exo.4.3
    nmosfet.dd.exp.gate_sweep.inp
    nmosfet.dd.exp.gate_sweep.exodiff
    nmosfet.dd.exp.gate_sweep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.exp.gate_sweep.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.exp.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.exp.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.exp.gate_sweep.exodiff
         nmosfet.dd.exp.gate_sweep.exo
         nmosfet.dd.exp.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_exp_inter PROPERTY LABELS
      DynamicTraps_exp_inter nightly debugexclude)
ENDIF()


####################################################
# Dynamic Traps on interface, gaussian distribution
SET(testName charon_mp_DynamicTraps_gau_inter)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME
    
  TEST_0 COPY_FILES_TO_TEST_DIR 
    nmosfet.dd.0.1.gau.exo.4.0
    nmosfet.dd.0.1.gau.exo.4.1
    nmosfet.dd.0.1.gau.exo.4.2
    nmosfet.dd.0.1.gau.exo.4.3
    nmosfet.dd.gau.gate_sweep.inp
    nmosfet.dd.gau.gate_sweep.exodiff
    nmosfet.dd.gau.gate_sweep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=nmosfet.dd.gau.gate_sweep.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosfet.dd.gau.gate_sweep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosfet.dd.gau.gate_sweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosfet.dd.gau.gate_sweep.exodiff
         nmosfet.dd.gau.gate_sweep.exo
         nmosfet.dd.gau.gate_sweep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_DynamicTraps_gau_inter PROPERTY LABELS
      DynamicTraps_gau_inter nightly debugexclude)
ENDIF()
