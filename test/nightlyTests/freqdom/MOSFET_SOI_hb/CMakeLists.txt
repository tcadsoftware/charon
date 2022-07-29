INCLUDE(TribitsAddAdvancedTest)

####################################################
# MOSFET Mobility coverage
SET(testName charon_mp_hb-Mobility_SOIMOSFET)
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
	  soimosfet.sg.hb.inp
          soimosfet.exodiff
	  solverSettings.inp
	  substrateDoping.inp
	  sweepOptions.inp
	  TpetraSolverSettings.inp

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=soimosfet.sg.hb.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=soimosfet.sg.hb.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_hb-Mobility_SOIMOSFET PROPERTY LABELS
     hb-Mobility_SOIMOSFET nightly debugexclude)
ENDIF()
