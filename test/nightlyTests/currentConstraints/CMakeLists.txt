INCLUDE(TribitsAddAdvancedTest)

####################################################
# Anode:    Constant Current
# Cathode:  Empty
SET(testName charon_mp_pndiode.dd.AC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AC.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.AC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.AC.inp.xml
      --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.AC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CC.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.CC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.CC.inp.xml
      --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.CC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.AR.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.AR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.AR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.AR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.CR.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.CR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.CR.inp.xml
      --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.CR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ACCR.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.ACCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ACCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ACCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCC.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.ARCC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ARCC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ARCC.output

  TEST_5 CMND python3
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
# Cathode:  Constant Current
# LOCA Householder Constraint Algorithm
SET(testName charon_mp_pndiode.loca.dd.ARCC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.loca.dd.ARCC.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.loca.dd.ARCC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.loca.dd.ARCC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.loca.dd.ARCC.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.loca.dd.ARCC.output
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
    TEST Charon_charon_mp_pndiode.loca.dd.ARCC
    PROPERTY LABELS pndiode.loca.dd.ARCC nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.loca.dd.ARCC
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

####################################################
# Anode:    Resistor Contact
# Cathode:  Resistor Contact
SET(testName charon_mp_pndiode.dd.ARCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ARCR.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.ARCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.dd.ARCR.inp.xml
      --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.dd.ARCR.output

  TEST_5 CMND python3
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

####################################################
# Anode:    Resistor Contact
# Cathode:  Resistor Contact
# LOCA Householder Constraint Algorithm
SET(testName charon_mp_pndiode.loca.dd.ARCR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.loca.dd.ARCR.inp
    pndiode.dd.exo
    pndiode.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.loca.dd.ARCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.dd.nem
         -m mesh=4
         pndiode.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.loca.dd.ARCR.inp.xml
      --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.loca.dd.ARCR.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.loca.dd.ARCR.output
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
    TEST Charon_charon_mp_pndiode.loca.dd.ARCR
    PROPERTY LABELS pndiode.loca.dd.ARCR nightly)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.loca.dd.ARCR
    PROPERTY ENVIRONMENT "ASAN_OPTIONS=detect_leaks=false")
ENDIF()

###############################################################################
# Test the Resistor Contact BC capability for a pndiode with tri mesh and SG
# scheme.
SET(testName charon_mp_pndiode.sg.AR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AR.inp
    pndiode.sg.AR.gold.exo
    pndiode.sg.AR.exodiff
    pndiode.sg.exo
    pndiode.sg.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.AR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.sg.nem
         -m mesh=4
         pndiode.sg.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.sg.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.sg.AR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.AR.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.sg.AR.output
         --anodeCurrent   -9.82385162e-04
         --cathodeCurrent 9.82385162e-04
         --tolerance      1e-10
         --anodeVoltage   5.90176148e-01
    PASS_REGULAR_EXPRESSION "Test Passed"

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result.pndiode.sg.AR.exo.4."
    NUM_MPI_PROCS 1

  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff
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
# Test the Resistor Contact BC capability for a pndiode with tri mesh and SG
# scheme with Householder LOCA algorithm for constraints.
SET(testName charon_mp_pndiode.loca.sg.AR)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.loca.sg.AR.inp
    pndiode.loca.sg.AR.gold.exo
    pndiode.loca.sg.AR.exodiff
    pndiode.sg.exo
    pndiode.sg.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.loca.sg.AR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.sg.nem
         -m mesh=4
         pndiode.sg.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.sg.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.loca.sg.AR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.loca.sg.AR.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification.py
         --filename       pndiode.loca.sg.AR.output
         --anodeCurrent   -9.82385162e-04
         --cathodeCurrent 9.82385162e-04
         --tolerance      1e-10
         --anodeVoltage   5.90176148e-01
    PASS_REGULAR_EXPRESSION "Test Passed"

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result.pndiode.loca.sg.AR.exo.4."
    NUM_MPI_PROCS 1

  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.loca.sg.AR.exodiff
         result.pndiode.loca.sg.AR.exo
         pndiode.loca.sg.AR.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_pndiode.loca.sg.AR
    PROPERTY LABELS pndiode.loca.sg.AR nightly)
ENDIF()

###############################################################################
# Test the Constant Current BC capability for a pndiode with tri mesh and SG
# scheme.
SET(testName charon_mp_pndiode.sg.AC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.AC.inp
    pndiode.sg.exo
    pndiode.sg.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.AC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o pndiode.sg.nem
         -m mesh=4
         pndiode.sg.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS pndiode.sg.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp
    NOEXEPREFIX
    DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=pndiode.sg.AC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.sg.AC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BJT1D.BC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BJT1D.BC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BJT1D.BC.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.BJT1D.BC.output
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --emitterCurrent   6.16854554e-3
         --tolerance        2.0e-10
         --baseVoltage      0.0
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.VoltageSweep.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BJT1D.BC.VoltageSweep.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp
    NOEXEPREFIX
    DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BJT1D.BC.VoltageSweep.inp.xml
         --current
    NUM_MPI_PROCS 4
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
# Test a case with both a Constant Current constraint and a BJT1D Base Contact
# boundary condition on the base.
# Uses the LOCA Householder Constraint algorithm
SET(testName charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep.householder)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BJT1D.BC.VoltageSweep.householder.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BJT1D.BC.VoltageSweep.householder.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp
    NOEXEPREFIX
    DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BJT1D.BC.VoltageSweep.householder.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BJT1D.BC.VoltageSweep.householder.output

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(
    TEST Charon_charon_mp_2n2222.dd.BJT1D.BC.VoltageSweep.householder
    PROPERTY LABELS 2n2222.dd.BJT1D.BC.VoltageSweep.householder nightly)
ENDIF()

###############################################################################
# Emitter:    Constant Current
# Base:       Empty
# Collector:  Empty
SET(testName charon_mp_2n2222.dd.EC)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.EC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.EC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.EC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.EC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.CC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.CC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.CC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ER.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ER.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ER.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ER.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.CR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.CR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.CR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.CR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ECBR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1
    
  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECBR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECBR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ECCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERBC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BCCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BCCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BCCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BCCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BRCC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BRCC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BRCC.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERCC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERCC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERCC.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERCC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        1.0e-7
         --emitterVoltage   -0.7
         --collectorVoltage 0.0
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.BRCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.BRCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.BRCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.BRCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERBR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ECBRCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ECBRCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ECBRCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ECBRCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBCCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERBCCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBCCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBCCR.output

  TEST_5 CMND python3
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCC.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERBRCC.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBRCC.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBRCC.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBRCC.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -6.16815485e-3
         --tolerance        1.0e-7
         --emitterVoltage   -0.7
         --baseVoltage      0.0
         --collectorVoltage 0.0
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
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    2n2222.dd.ERBRCR.inp
    2n2222.dd.exo
    2n2222.dd.pex

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i 2n2222.dd.ERBRCR.inp

  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/nem_slice/nem_slice
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS -e
         -l inertial
         -c
         -o 2n2222.dd.nem
         -m mesh=4
         2n2222.dd.exo
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/nem_spread/nem_spread
    NOEXEPREFIX
    NOEXESUFFIX
    ARGS 2n2222.dd.pex
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS --i=2n2222.dd.ERBRCR.inp.xml
         --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE 2n2222.dd.ERBRCR.output

  TEST_5 CMND python3
    ARGS ${CMAKE_SOURCE_DIR}/scripts/utils/current_verification_bjt.py
         --filename         2n2222.dd.ERBRCR.output
         --emitterCurrent   6.16854554e-3
         --baseCurrent      -3.90689884e-7
         --collectorCurrent -0.0061681509
         --tolerance        5.0e-6
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
