INCLUDE(TribitsAddAdvancedTest)

####################################################
# Anode:    Constant Current
# Cathode:  Empty
SET(testName charon_mp_pndiode.dd.AC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AC.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.AC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.AC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.AC.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-11
         --anodeVoltage   4.46464522e-1
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.AC
    PROPERTY LABELS pndiode.dd.AC nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.AC
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Empty
# Cathode:  Constant Current
SET(testName charon_mp_pndiode.dd.CC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CC.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.CC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.CC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.CC.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-9
         --cathodeVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.CC
    PROPERTY LABELS pndiode.dd.CC nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.CC
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Resistor Contact
# Cathode:  Empty
SET(testName charon_mp_pndiode.dd.AR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AR.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.AR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.AR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.AR.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-11
         --anodeVoltage   4.46464522e-1
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.AR
    PROPERTY LABELS pndiode.dd.AR nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.AR
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Empty
# Cathode:  Resistor Contact
SET(testName charon_mp_pndiode.dd.CR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CR.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.CR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.CR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.CR.output
         --anodeCurrent   -4.0e-6
         --cathodeCurrent 4.0e-6
         --tolerance      5.0e-10
         --cathodeVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.CR
    PROPERTY LABELS pndiode.dd.CR nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.CR
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Constant Current
# Cathode:  Resistor Contact
SET(testName charon_mp_pndiode.dd.ACCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ACCR.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ACCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ACCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.ACCR.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-8
         --anodeVoltage   4.46464522e-1
         --cathodeVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ACCR
    PROPERTY LABELS pndiode.dd.ACCR nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ACCR
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Resistor Contact
# Cathode:  Constant Current
SET(testName charon_mp_pndiode.dd.ARCC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCC.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ARCC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ARCC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.ARCC.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-8
         --anodeVoltage   4.46464522e-1
         --cathodeVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ARCC
    PROPERTY LABELS pndiode.dd.ARCC nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ARCC
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Resistor Contact
# Cathode:  Resistor Contact
SET(testName charon_mp_pndiode.dd.ARCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCR.yaml
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ARCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ARCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.dd.ARCR.output
         --anodeCurrent   -4e-6
         --cathodeCurrent 4e-6
         --tolerance      1e-10
         --anodeVoltage   4.46464522e-1
         --cathodeVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ARCR
    PROPERTY LABELS pndiode.dd.ARCR nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.dd.ARCR
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

###############################################################################
# Test the Resistor Contact BC capability for a pndiode with tri mesh and SG
# scheme.
SET(testName charon_mp_pndiode.sg.AR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AR.yaml
    pndiode.sg.AR.gold.exo
    pndiode.sg.AR.exodiff
    pndiode.sg.exo
    pndiode.sg.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.sg.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.sg.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.sg.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.sg.AR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.AR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.AR.output
         --anodeCurrent   -9.82385162e-04
         --cathodeCurrent 9.82385162e-04
         --tolerance      1e-10
         --anodeVoltage   5.90176148e-01
    PASS_REGULAR_EXPRESSION "Test Passed"

  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result.pndiode.sg.AR.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.sg.AR.exodiff
         result.pndiode.sg.AR.exo
         pndiode.sg.AR.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.sg.AR
    PROPERTY LABELS pndiode.sg.AR nightly)
ENDIF()

###############################################################################
# Test the Constant Current BC capability for a pndiode with tri mesh and SG
# scheme.
SET(testName charon_mp_pndiode.sg.AC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AC.yaml
    pndiode.sg.exo
    pndiode.sg.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.sg.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         pndiode.sg.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.sg.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp
    NOEXEPREFIX
    DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.sg.AC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.AC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.AC.output
         --anodeCurrent   -3.36926779e-5
         --cathodeCurrent 3.36926779e-5
         --tolerance      1e-10
         --anodeVoltage   0.5
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.sg.AC
    PROPERTY LABELS pndiode.sg.AC nightly)
ENDIF()

###############################################################################
# Test a case with both a Constant Current constraint and a BJT1D Base Contact
# boundary condition on the base.
SET(testName charon_mp_2n2222.dd.BJT1D.BC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BJT1D.BC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BJT1D.BC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BJT1D.BC.output
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --emitterCurrent   6.16854554e-3
         --tolerance        1e-10
         --baseVoltage      1.91262625e-10
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BJT1D.BC
    PROPERTY LABELS 2n2222.dd.BJT1D.BC nightly)
ENDIF()

###############################################################################
# Test a case with both a Constant Current constraint and a BJT1D Base Contact
# boundary condition on the base.
SET(testName charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.VoltageSweep.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp
    NOEXEPREFIX
    DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BJT1D.BC.VoltageSweep.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BJT1D.BC.VoltageSweep.output

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep
    PROPERTY LABELS 2n2222.dd.BJT1D.BC.VoltageSweep nightly)
ENDIF()

###############################################################################
# Emitter:    Constant Current
# Base:       Empty
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.EC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.EC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.EC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.EC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.EC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.EC
    PROPERTY LABELS 2n2222.dd.EC nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Constant Current
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.BC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --baseVoltage      0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BC
    PROPERTY LABELS 2n2222.dd.BC nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Empty
# Collector:  Constant Current
SET(testName charon_mp_2n2222.dd.CC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.CC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.CC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.CC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        6e-9
         --collectorVoltage 6.15489283e-09
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.CC
    PROPERTY LABELS 2n2222.dd.CC nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Empty
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.ER)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ER.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ER.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ER.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ER.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ER
    PROPERTY LABELS 2n2222.dd.ER nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Resistor Contact
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.BR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --baseVoltage      0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BR
    PROPERTY LABELS 2n2222.dd.BR nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Empty
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.CR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.CR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.CR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.CR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.CR
    PROPERTY LABELS 2n2222.dd.CR nightly)
ENDIF()

###############################################################################
# Emitter:    Constant Current
# Base:       Resistor Contact
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.ECBR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECBR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECBR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ECBR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
         --baseVoltage      0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ECBR
    PROPERTY LABELS 2n2222.dd.ECBR nightly)
ENDIF()

###############################################################################
# Emitter:    Constant Current
# Base:       Empty
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.ECCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ECCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ECCR
    PROPERTY LABELS 2n2222.dd.ECCR nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Constant Current
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.ERBC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-8
         --emitterVoltage   -0.7
         --baseVoltage      0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERBC
    PROPERTY LABELS 2n2222.dd.ERBC nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Constant Current
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.BCCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BCCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BCCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BCCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BCCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        1e-8
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BCCR
    PROPERTY LABELS 2n2222.dd.BCCR nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Resistor Contact
# Collector:  Constant Current
SET(testName charon_mp_2n2222.dd.BRCC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BRCC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BRCC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BRCC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        6e-9
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BRCC
    PROPERTY LABELS 2n2222.dd.BRCC nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Empty
# Collector:  Constant Current
SET(testName charon_mp_2n2222.dd.ERCC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERCC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERCC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERCC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        7e-9
         --emitterVoltage   -0.7
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERCC
    PROPERTY LABELS 2n2222.dd.ERCC nightly)
ENDIF()

###############################################################################
# Emitter:    Empty
# Base:       Resistor Contact
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.BRCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BRCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BRCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BRCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BRCR
    PROPERTY LABELS 2n2222.dd.BRCR nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Empty
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.ERCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERCR
    PROPERTY LABELS 2n2222.dd.ERCR nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Resistor Contact
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.ERBR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
         --baseVoltage      0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERBR
    PROPERTY LABELS 2n2222.dd.ERBR nightly)
ENDIF()

###############################################################################
# Emitter:    Constant Current
# Base:       Resistor Contact
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.ECBRCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBRCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECBRCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECBRCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ECBRCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        3e-9
         --emitterVoltage   -0.7
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ECBRCR
    PROPERTY LABELS 2n2222.dd.ECBRCR nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Constant Current
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.ERBCCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBCCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBCCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBCCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBCCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        1e-4
         --emitterVoltage   -0.7
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERBCCR
    PROPERTY LABELS 2n2222.dd.ERBCCR nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Resistor Contact
# Collector:  Constant Current
SET(testName charon_mp_2n2222.dd.ERBRCC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCC.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBRCC.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBRCC.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBRCC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        1e-8
         --emitterVoltage   -0.7
         --baseVoltage      0
         --collectorVoltage 0
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERBRCC
    PROPERTY LABELS 2n2222.dd.ERBRCC nightly)
ENDIF()

###############################################################################
# Emitter:    Resistor Contact
# Base:       Resistor Contact
# Collector:  Resistor Contact
SET(testName charon_mp_2n2222.dd.ERBRCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCR.yaml
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=${MPI_EXEC_MAX_NUMPROCS}
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBRCR.yaml
         --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBRCR.output

  TEST_4 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBRCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -0.0061681509
         --tolerance        3e-9
         --emitterVoltage   -0.699999831
         --baseVoltage      1.64131489e-07
         --collectorVoltage 4.13088498e-06
    PASS_REGULAR_EXPRESSION "Test Passed"

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.ERBRCR
    PROPERTY LABELS 2n2222.dd.ERBRCR nightly)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.AC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AC.inp
    pndiode.sg.AC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i pndiode.sg.AC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.AC.inp.xml pndiode.sg.AC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.AC_input PROPERTY LABELS
    pndiode.sg.AC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.sg.AR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AR.inp
    pndiode.sg.AR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i pndiode.sg.AR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.sg.AR.inp.xml pndiode.sg.AR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.AR_input PROPERTY LABELS
    pndiode.sg.AR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.EC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.EC.inp
    2n2222.dd.EC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.EC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.EC.inp.xml 2n2222.dd.EC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.EC_input PROPERTY LABELS
    2n2222.dd.EC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.CC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CC.inp
    2n2222.dd.CC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.CC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.CC.inp.xml 2n2222.dd.CC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.CC_input PROPERTY LABELS
    2n2222.dd.CC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BC.inp
    2n2222.dd.BC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.BC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BC.inp.xml 2n2222.dd.BC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BC_input PROPERTY LABELS
    2n2222.dd.BC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ER_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ER.inp
    2n2222.dd.ER.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ER.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ER.inp.xml 2n2222.dd.ER.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ER_input PROPERTY LABELS
    2n2222.dd.ER_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BR.inp
    2n2222.dd.BR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.BR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BR.inp.xml 2n2222.dd.BR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BR_input PROPERTY LABELS
    2n2222.dd.BR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.CR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CR.inp
    2n2222.dd.CR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.CR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.CR.inp.xml 2n2222.dd.CR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.CR_input PROPERTY LABELS
    2n2222.dd.CR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ECBR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBR.inp
    2n2222.dd.ECBR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ECBR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ECBR.inp.xml 2n2222.dd.ECBR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ECBR_input PROPERTY LABELS
    2n2222.dd.ECBR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ECCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECCR.inp
    2n2222.dd.ECCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ECCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ECCR.inp.xml 2n2222.dd.ECCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ECCR_input PROPERTY LABELS
    2n2222.dd.ECCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERBC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBC.inp
    2n2222.dd.ERBC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERBC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERBC.inp.xml 2n2222.dd.ERBC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERBC_input PROPERTY LABELS
    2n2222.dd.ERBC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BCCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BCCR.inp
    2n2222.dd.BCCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.BCCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BCCR.inp.xml 2n2222.dd.BCCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BCCR_input PROPERTY LABELS
    2n2222.dd.BCCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BRCC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCC.inp
    2n2222.dd.BRCC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.BRCC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BRCC.inp.xml 2n2222.dd.BRCC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BRCC_input PROPERTY LABELS
    2n2222.dd.BRCC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERCC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCC.inp
    2n2222.dd.ERCC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERCC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERCC.inp.xml 2n2222.dd.ERCC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERCC_input PROPERTY LABELS
    2n2222.dd.ERCC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BRCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCR.inp
    2n2222.dd.BRCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.BRCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BRCR.inp.xml 2n2222.dd.BRCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BRCR_input PROPERTY LABELS
    2n2222.dd.BRCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCR.inp
    2n2222.dd.ERCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERCR.inp.xml 2n2222.dd.ERCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERCR_input PROPERTY LABELS
    2n2222.dd.ERCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERBR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBR.inp
    2n2222.dd.ERBR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERBR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERBR.inp.xml 2n2222.dd.ERBR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERBR_input PROPERTY LABELS
    2n2222.dd.ERBR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ECBRCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBRCR.inp
    2n2222.dd.ECBRCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ECBRCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ECBRCR.inp.xml 2n2222.dd.ECBRCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ECBRCR_input PROPERTY LABELS
    2n2222.dd.ECBRCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERBCCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBCCR.inp
    2n2222.dd.ERBCCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERBCCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERBCCR.inp.xml 2n2222.dd.ERBCCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERBCCR_input PROPERTY LABELS
    2n2222.dd.ERBCCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERBRCC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCC.inp
    2n2222.dd.ERBRCC.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERBRCC.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERBRCC.inp.xml 2n2222.dd.ERBRCC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERBRCC_input PROPERTY LABELS
    2n2222.dd.ERBRCC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.ERBRCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCR.inp
    2n2222.dd.ERBRCR.inp.gold.xml

  TEST_1 CMND ../../../../src/charonInterpreter 
    ARGS -i 2n2222.dd.ERBRCR.inp

  TEST_2 CMND ../../../../src/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.ERBRCR.inp.xml 2n2222.dd.ERBRCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.ERBRCR_input PROPERTY LABELS
    2n2222.dd.ERBRCR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.ARCC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCC.inp
    pndiode.dd.ARCC.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.ARCC.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.ARCC.inp.xml pndiode.dd.ARCC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.ARCC_input PROPERTY LABELS
    pndiode.dd.ARCC inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.ARCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCR.inp
    pndiode.dd.ARCR.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.ARCR.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.ARCR.inp.xml pndiode.dd.ARCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.ARCR_input PROPERTY LABELS
    pndiode.dd.ARCR inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.AC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AC.inp
    pndiode.dd.AC.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.AC.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.AC.inp.xml pndiode.dd.AC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.AC_input PROPERTY LABELS
    pndiode.dd.AC inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.CC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CC.inp
    pndiode.dd.CC.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.CC.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.CC.inp.xml pndiode.dd.CC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.CC_input PROPERTY LABELS
    pndiode.dd.CC_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.AR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AR.inp
    pndiode.dd.AR.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.AR.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.AR.inp.xml pndiode.dd.AR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.AR_input PROPERTY LABELS
    pndiode.dd.AR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.CR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CR.inp
    pndiode.dd.CR.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.CR.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.CR.inp.xml pndiode.dd.CR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.CR_input PROPERTY LABELS
    pndiode.dd.CR_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.dd.ACCR_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ACCR.inp
    pndiode.dd.ACCR.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.ACCR.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.dd.ACCR.inp.xml pndiode.dd.ACCR.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.ACCR_input PROPERTY LABELS
    pndiode.dd.ACCR_input inputVerification)
ENDIF()


##########################################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BJT1D.BC_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.inp
    2n2222.dd.BJT1D.BC.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BJT1D.BC.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BJT1D.BC.inp.xml 2n2222.dd.BJT1D.BC.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BJT1D.BC_input PROPERTY LABELS
    2n2222.dd.BJT1D.BC_input inputVerification)
ENDIF()


##########################################################################
# Interpreter test
SET(testName charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.VoltageSweep.inp
    2n2222.dd.BJT1D.BC.VoltageSweep.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BJT1D.BC.VoltageSweep.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS 2n2222.dd.BJT1D.BC.VoltageSweep.inp.xml 2n2222.dd.BJT1D.BC.VoltageSweep.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep_input PROPERTY LABELS
    2n2222.dd.BJT1D.BC.VoltageSweep_input inputVerification)
ENDIF()



