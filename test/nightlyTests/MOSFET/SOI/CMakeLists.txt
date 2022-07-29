INCLUDE(TribitsAddAdvancedTest)

####################################################
# MOSFET Mobility coverage
SET(testName charon_mp_Mobility_SOIMOSFET)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  ENVIRONMENT APREPRO_PATH=${SEACAS_BINARY_DIR}/applications/aprepro

  TEST_0 COPY_FILES_TO_TEST_DIR
  	  bodyDoping.inp
	  gatePolyDoping.inp
	  geometry.inp
	  insulatorProperties.inp
	  mobility.inp
	  scaling.inp
	  soimosfet.sg.out.gold.exo
	  soimosfet.sg.in.exo.4.0
	  soimosfet.sg.in.exo.4.1
	  soimosfet.sg.in.exo.4.2
	  soimosfet.sg.in.exo.4.3
	  soimosfet.sg.inp
          soimosfet.exodiff
	  solverSettings.inp
	  substrateDoping.inp
	  sweepOptions.inp
	  TpetraSolverSettings.inp

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=soimosfet.sg.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=soimosfet.sg.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "soimosfet.sg.out.exo.4. "
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         soimosfet.exodiff
         soimosfet.sg.out.exo
         soimosfet.sg.out.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Mobility_SOIMOSFET PROPERTY LABELS
     Mobility_SOIMOSFET nightly debugexclude)
ENDIF()
