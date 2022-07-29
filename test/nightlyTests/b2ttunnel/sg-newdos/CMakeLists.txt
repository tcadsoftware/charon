INCLUDE(TribitsAddAdvancedTest)

#######################################################
# test the band-to-trap tunneling capability for 
# Hole Tunneling Model = Schenk NewDOS and FEM-SUPG
#######################################################
SET(testName charon_mp_npn-hbt.sg.trap.short.newdos)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME
  
  TEST_0 COPY_FILES_TO_TEST_DIR
    GaAs_LowField_EMob.txt
    GaAs_LowField_HMob.txt
    npn-hbt.sg.telt.0.8V.exo.4.0
    npn-hbt.sg.telt.0.8V.exo.4.1
    npn-hbt.sg.telt.0.8V.exo.4.2    
    npn-hbt.sg.telt.0.8V.exo.4.3
    npn-hbt.sg.trap.newdos.inp
    currents-loca.gold.dat
    npn-hbt.sg.trap.newdos.gold.exo
    npn-hbt.exodiff

  # Clean up any existing current file so this run starts with a new file
  TEST_1 CMND rm ARGS -f currents-loca.dat
    PASS_ANY

  # Run the interpreter
  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i npn-hbt.sg.trap.newdos.inp

  # Run the Gummel LOCA sweep (5 voltage points for short time)
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=npn-hbt.sg.trap.newdos.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the parallel exodus
  TEST_4  EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX
    NOEXESUFFIX
    ARGS  "-auto" "npn-hbt.sg.trap.newdos.exo.4."
    NUM_MPI_PROCS 1
    
  # Compare the exodus
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f 
         npn-hbt.exodiff 
         npn-hbt.sg.trap.newdos.exo 
         npn-hbt.sg.trap.newdos.gold.exo
    NUM_MPI_PROCS 1
  
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_npn-hbt.sg.trap.short.newdos PROPERTY LABELS
     npn-hbt.sg.trap.short.newdos nightly)
ENDIF()  
