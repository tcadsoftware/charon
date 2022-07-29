INCLUDE(TribitsAddAdvancedTest)

##########################################################
# Test the Band2Band_Local_Tunneling Hurkx model
##########################################################

SET(testName charon_mp_b2btunnel_Hurkx)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    Esaki.diode.jou Esaki.diode.nlp.inp Esaki.diode.equ.mb.inp Esaki.diode.equ.fd.inp
    Esaki.diode.loca.dd.inp Esaki.diode.loca.sg.inp 
    Esaki.diode.exo.4.0 Esaki.diode.exo.4.1 Esaki.diode.exo.4.2 Esaki.diode.exo.4.3
    Esaki.diode.loca.dd.gold.exo Esaki.diode.loca.sg.gold.exo
    Esaki.diode.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.nlp.inp

  # Nonlinear Poisson simulation
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.equ.mb.inp

  # Equilibrium MB simulation
  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.equ.mb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.equ.fd.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.equ.fd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_7 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.loca.dd.inp

  # Loca simulation (FEM-SUPG)
  TEST_8 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.loca.dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.loca.sg.inp

  # Loca simulation (CVFEM-SG)
  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.loca.sg.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS "-auto" "Esaki.diode.loca.dd.exo.4." 
    NUM_MPI_PROCS 1

  # Combine the exodus files for comparison
  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "Esaki.diode.loca.sg.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_13 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         Esaki.diode.exodiff
         Esaki.diode.loca.dd.exo
         Esaki.diode.loca.dd.gold.exo
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_14 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         Esaki.diode.exodiff
         Esaki.diode.loca.sg.exo
         Esaki.diode.loca.sg.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
)

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_b2btunnel_Hurkx PROPERTY LABELS b2btunnel_Hurkx nightly debugexclude)
ENDIF()


##########################################################
# Test the Band2Band_Local_Tunneling Kane model
##########################################################

SET(testName charon_mp_b2btunnel_Kane)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  # Copy files
  TEST_0 COPY_FILES_TO_TEST_DIR
    Esaki.diode.jou Esaki.diode.nlp.inp Esaki.diode.equ.mb.inp Esaki.diode.equ.fd.inp
    Esaki.diode.kane.dd.useF0.inp Esaki.diode.kane.sg.useF0.inp
    Esaki.diode.exo.4.0 Esaki.diode.exo.4.1 Esaki.diode.exo.4.2 Esaki.diode.exo.4.3
    Esaki.diode.kane.dd.useF0.gold.exo Esaki.diode.kane.sg.useF0.gold.exo
    Esaki.diode.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.nlp.inp

  # Nonlinear Poisson simulation
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.equ.mb.inp

  # Equilibrium MB simulation
  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.equ.mb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.equ.fd.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.equ.fd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_7 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.kane.dd.useF0.inp

  # Loca simulation (FEM-SUPG)
  TEST_8 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.kane.dd.useF0.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_9 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i Esaki.diode.kane.sg.useF0.inp

  # Loca simulation (CVFEM-SG)
  TEST_10 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=Esaki.diode.kane.sg.useF0.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # Combine the exodus files for comparison
  TEST_11 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS "-auto" "Esaki.diode.kane.dd.useF0.exo.4."
    NUM_MPI_PROCS 1

  # Combine the exodus files for comparison
  TEST_12 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "Esaki.diode.kane.sg.useF0.exo.4."
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_13 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         Esaki.diode.exodiff
         Esaki.diode.kane.dd.useF0.exo
         Esaki.diode.kane.dd.useF0.gold.exo
    NUM_MPI_PROCS 1

  # Compare simulation results
  TEST_14 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         Esaki.diode.exodiff
         Esaki.diode.kane.sg.useF0.exo
         Esaki.diode.kane.sg.useF0.gold.exo
   NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
)

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_b2btunnel_Kane PROPERTY LABELS b2btunnel_Kane nightly debugexclude)
ENDIF()




