INCLUDE(TribitsAddAdvancedTest)

##########################################################
# These tests simply execute various capabilities to make
# sure they run. Add comparison, otherwise, no indication 
# of test failure.
##########################################################

##################################################################
# Test FEM NLP on a pmoscap at 300 K
##################################################################
SET(testName charon_mp_pmoscap.nlp)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pmoscap.nlp.inp
    pmoscap.exo.4.0
    pmoscap.exo.4.1
    pmoscap.exo.4.2
    pmoscap.exo.4.3
    pmoscap.nlp.gold.exo
    pndiode.nlp.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pmoscap.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pmoscap.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pmoscap.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         pmoscap.nlp.exo
         pmoscap.nlp.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pmoscap.nlp PROPERTY LABELS pmoscap.nlp nightly)
ENDIF()


####################################################################
# Test SGCVFEM NLP on a pmoscap at 4 K using Fermi-Dirac statistics
####################################################################
SET(testName charon_mp_pmoscap.nlp.4K)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pmoscap.nlp.4K.inp
    pmoscap.exo.4.0
    pmoscap.exo.4.1
    pmoscap.exo.4.2
    pmoscap.exo.4.3
    pmoscap.nlp.4K.gold.exo
    pndiode.nlp.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pmoscap.nlp.4K.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pmoscap.nlp.4K.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pmoscap.nlp.4K.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         pmoscap.nlp.4K.exo
         pmoscap.nlp.4K.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pmoscap.nlp.4K PROPERTY LABELS pmoscap.nlp.4K nightly)
ENDIF()


##################################################################
# Test FEM NLP on a pn diode with uniform doping
##################################################################
SET(testName charon_mp_pndiode.nlp)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.exodiff
    potential_pndiode.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.nlp.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "potential_pndiode.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         potential_pndiode.exo
         potential_pndiode.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp PROPERTY LABELS mp_pndiode.nlp nightly)
ENDIF()


##################################################################
# Test FEM NLP on a pn diode with linear doping
##################################################################
SET(testName charon_mp_pndiode.nlp.linear)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.linear.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.exodiff
    potential_pndiode_linear.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.nlp.linear.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.linear.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "potential_pndiode_linear.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         potential_pndiode_linear.exo
         potential_pndiode_linear.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.linear PROPERTY LABELS 
    pndiode.nlp.linear nightly)
ENDIF()


##################################################################
# Test FEM NLP on a pn diode with Gauss doping
##################################################################
SET(testName charon_mp_pndiode.nlp.gauss)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.gauss.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.exodiff
    potential_pndiode_gauss.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.nlp.gauss.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.gauss.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "potential_pndiode_gauss.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         potential_pndiode_gauss.exo
         potential_pndiode_gauss.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.gauss PROPERTY LABELS
    pndiode.nlp.gauss nightly)
ENDIF()


##################################################################
# Test FEM NLP on a pn diode with Erfc doping
##################################################################
SET(testName charon_mp_pndiode.nlp.erfc)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.erfc.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.exodiff
    potential_pndiode_erfc.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.nlp.erfc.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.erfc.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "potential_pndiode_erfc.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         potential_pndiode_erfc.exo
         potential_pndiode_erfc.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.erfc PROPERTY LABELS
    pndiode.nlp.erfc nightly)
ENDIF()


##################################################################
# Test FEM DD on a pn diode with constant SRH life time
##################################################################
SET(testName charon_mp_pndiode.dd.equ)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.equ.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_dd_equ.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.equ.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_equ.exo
         result_pndiode_dd_equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.equ PROPERTY LABELS
    pndiode.dd.equ nightly)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.equ PROPERTY DEPENDS
    Charon_charon_mp_pndiode.nlp)
ENDIF()


##################################################################
# Test FEM DD on a pn diode with Auger and conc-dependent SRH
##################################################################
SET(testName charon_mp_pndiode.dd.recomb)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.recomb.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_dd_recomb.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.recomb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.recomb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_recomb.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_recomb.exo
         result_pndiode_dd_recomb.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.recomb PROPERTY LABELS
    pndiode.dd.recomb nightly)
ENDIF()


##############################################################################
# Test FEM DD on a pn diode with Philips-Thomas and Masetti mobility models
##############################################################################
SET(testName charon_mp_pndiode.dd.lfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.lfmob.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_dd_lfmob.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.lfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.lfmob.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_lfmob.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_lfmob.exo
         result_pndiode_dd_lfmob.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.lfmob PROPERTY LABELS
    pndiode.dd.lfmob nightly)
ENDIF()


##############################################################################
# Test FEM DD on a pn diode with incomplete ionization model
##############################################################################
SET(testName charon_mp_pndiode.dd.ioniz)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ioniz.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    accIncmplIoniz.txt
    donIncmplIoniz.txt
    pndiode.nlp.exodiff
    result_pndiode_dd_ioniz.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.ioniz.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.ioniz.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_ioniz.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_ioniz.exo
         result_pndiode_dd_ioniz.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.ioniz PROPERTY LABELS
    pndiode.dd.ioniz nightly)
ENDIF()


##############################################################################
# Test FEM DD on a pn diode with Lucent and Philips-Thomas mobility models
##############################################################################
SET(testName charon_mp_pndiode.dd.hfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.hfmob.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_dd_hfmob.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.hfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.hfmob.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_hfmob.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_hfmob.exo
         result_pndiode_dd_hfmob.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.hfmob PROPERTY LABELS
    pndiode.dd.hfmob nightly)
ENDIF()


##############################################################################
# Test FEM DD on a pn diode with Van Overstraeten avalanche model
##############################################################################
SET(testName charon_mp_pndiode.dd.ava)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.ava.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_dd_ava.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.dd.ava.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.ava.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_dd_ava.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_dd_ava.exo
         result_pndiode_dd_ava.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.ava PROPERTY LABELS
    pndiode.dd.ava nightly)
ENDIF()


##############################################################################
# Test SGCVFEM DD on a pn diode at equilibrium
##############################################################################
SET(testName charon_mp_pndiode.sg.equ)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.equ.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_sg_equ.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.equ.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.equ.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_sg_equ.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_sg_equ.exo
         result_pndiode_sg_equ.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.equ PROPERTY LABELS
    pndiode.sg.equ nightly)
#  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.equ PROPERTY DEPENDS
#    Charon_charon_mp_pndiode.nlp)
ENDIF()


##############################################################################
# Test SGCVFEM DD on a pn diode with Auger and conc-dependent SRH models
##############################################################################
SET(testName charon_mp_pndiode.sg.recomb)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.recomb.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_sg_recomb.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.recomb.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.recomb.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_sg_recomb.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_sg_recomb.exo
         result_pndiode_sg_recomb.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.recomb PROPERTY LABELS
    pndiode.sg.recomb nightly)
ENDIF()


################################################################################
# Test SGCVFEM DD on a pn diode with Philips-Thomas and Masetti mobility models 
################################################################################
SET(testName charon_mp_pndiode.sg.lfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lfmob.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_sg_lfmob.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.lfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.lfmob.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_sg_lfmob.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_sg_lfmob.exo
         result_pndiode_sg_lfmob.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lfmob PROPERTY LABELS
    pndiode.sg.lfmob nightly)
ENDIF()


################################################################################
# Test SGCVFEM DD on a pn diode with Lucent and Philips-Thomas mobility models 
################################################################################
SET(testName charon_mp_pndiode.sg.hfmob)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.hfmob.inp
    result_pndiode_equ.exo.4.0
    result_pndiode_equ.exo.4.1
    result_pndiode_equ.exo.4.2
    result_pndiode_equ.exo.4.3
    pndiode.nlp.exodiff
    result_pndiode_sg_hfmob.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.sg.hfmob.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.hfmob.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_sg_hfmob.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         result_pndiode_sg_hfmob.exo
         result_pndiode_sg_hfmob.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.hfmob PROPERTY LABELS
    pndiode.sg.hfmob nightly)
ENDIF()


##################################################################
# test Gauss decay for Uniform and File1D doping
##################################################################

SET(testName charon_mp_pndiode.nlp.gauss-decay)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.gauss-decay.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.gauss-decay.gold.exo
    pndiode.nlp.exodiff
    test_dop.txt

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i  pndiode.nlp.gauss-decay.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.gauss-decay.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.nlp.gauss-decay.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         pndiode.nlp.gauss-decay.exo
         pndiode.nlp.gauss-decay.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.gauss-decay PROPERTY LABELS pndiode.nlp.gauss-decay nightly)
ENDIF()


##################################################################
# test to read multiple 1d and 2d doping files
##################################################################

SET(testName charon_mp_pndiode.nlp.multi-files)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.nlp.multi-files.inp
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.multi-files.gold.exo
    pndiode.nlp.exodiff
    test1_dop1d.txt
    test2_dop1d.txt
    test3_dop2d.txt
    test4_dop2d.txt

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS  -i  pndiode.nlp.multi-files.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.multi-files.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.nlp.multi-files.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.nlp.exodiff
         pndiode.nlp.multi-files.exo
         pndiode.nlp.multi-files.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.nlp.multi-files PROPERTY LABELS pndiode.nlp.multi-files nightly)
ENDIF()


