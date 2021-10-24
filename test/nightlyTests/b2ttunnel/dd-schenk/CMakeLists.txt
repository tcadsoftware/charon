INCLUDE(TribitsAddAdvancedTest)

#######################################################
# test the band-to-trap tunneling capability for 
# Hole Tunneling Model = Schenk ConstFDOS and FEM-SUPG
#######################################################
SET(testName charon_mp_npn-hbt.dd.trap.schenk)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME
  TIMEOUT 1200
  
  TEST_0 COPY_FILES_TO_TEST_DIR
    GaAs_LowField_EMob.txt
    GaAs_LowField_HMob.txt
    npn-hbt.dd.telt.0.8V.exo.4.0
    npn-hbt.dd.telt.0.8V.exo.4.1
    npn-hbt.dd.telt.0.8V.exo.4.2    
    npn-hbt.dd.telt.0.8V.exo.4.3
    npn-hbt.dd.trap.schenk.xml
    currents-loca.gold.dat
    npn-hbt.dd.trap.schenk.gold.exo
    npn-hbt.exodiff

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY
  
  # Run the Gummel LOCA sweep (5 voltage points for short time)
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=npn-hbt.dd.trap.schenk.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_3  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "npn-hbt.dd.trap.schenk.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1
    
  # Compare the exodus
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 
         npn-hbt.exodiff 
         npn-hbt.dd.trap.schenk.exo 
         npn-hbt.dd.trap.schenk.gold.exo
    NUM_MPI_PROCS 1
  
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_npn-hbt.dd.trap.schenk PROPERTY LABELS
     npn-hbt.dd.trap.schenk nightly)
ENDIF()  


####################################################
# Interpreter test
SET(testName charon_mp_npn-hbt.dd.trap.schenk_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    npn-hbt.dd.trap.schenk.inp
    npn-hbt.dd.trap.schenk.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i npn-hbt.dd.trap.schenk.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS npn-hbt.dd.trap.schenk.inp.xml npn-hbt.dd.trap.schenk.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_npn-hbt.dd.trap.schenk_input PROPERTY LABELS
    npn-hbt.dd.trap.schenk_input inputVerification)
ENDIF()
    
