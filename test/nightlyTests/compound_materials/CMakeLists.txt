INCLUDE(TribitsAddAdvancedTest)

####################################################
# Simple Al(x)Ga(1-x)N resistor example
SET(testName charon_mp_AlGaN_resist)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    algan_resist_msh.jou 
    algan_resist_msh.exo.4.0
    algan_resist_msh.exo.4.1
    algan_resist_msh.exo.4.2
    algan_resist_msh.exo.4.3
    algan_resist.nlp.inp
    algan_resist.nlp.exodiff
    algan_resist.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=algan_resist.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=algan_resist.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "algan_resist.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         algan_resist.nlp.exodiff
         algan_resist.nlp.exo
         algan_resist.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_AlGaN_resist PROPERTY LABELS
     AlGaN_resist nightly)
ENDIF()

####################################################
# Simple Si(1-x)Ge(x) resistor example
SET(testName charon_mp_SiGe_resist)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    sige_resist_msh.jou 
    sige_resist_msh.exo.4.0
    sige_resist_msh.exo.4.1
    sige_resist_msh.exo.4.2
    sige_resist_msh.exo.4.3
    sige_resist.nlp.inp
    sige_resist.nlp.exodiff
    sige_resist.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=sige_resist.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=sige_resist.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "sige_resist.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         sige_resist.nlp.exodiff
         sige_resist.nlp.exo
         sige_resist.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_SiGe_resist PROPERTY LABELS
     SiGe_resist nightly)
ENDIF()


####################################################
# Simple Si(1-x)Ge(x) resistor example with gaussian mf
SET(testName charon_mp_SiGe_resist_gauss)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    sige_resist_gauss_msh.jou 
    sige_resist_gauss_msh.exo.4.0
    sige_resist_gauss_msh.exo.4.1
    sige_resist_gauss_msh.exo.4.2
    sige_resist_gauss_msh.exo.4.3
    sige_resist_gauss.nlp.inp
    sige_resist_gauss.nlp.exodiff
    sige_resist_gauss.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=sige_resist_gauss.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=sige_resist_gauss.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "sige_resist_gauss.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         sige_resist_gauss.nlp.exodiff
         sige_resist_gauss.nlp.exo
         sige_resist_gauss.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_SiGe_resist_gauss PROPERTY LABELS
     SiGe_resist_gauss nightly)
ENDIF()

####################################################
# Simple Si(1-x)Ge(x) resistor example with erfc mf
SET(testName charon_mp_SiGe_resist_erfc)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    sige_resist_erfc_msh.jou 
    sige_resist_erfc_msh.exo.4.0
    sige_resist_erfc_msh.exo.4.1
    sige_resist_erfc_msh.exo.4.2
    sige_resist_erfc_msh.exo.4.3
    sige_resist_erfc.nlp.inp
    sige_resist_erfc.nlp.exodiff
    sige_resist_erfc.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=sige_resist_erfc.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=sige_resist_erfc.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "sige_resist_erfc.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         sige_resist_erfc.nlp.exodiff
         sige_resist_erfc.nlp.exo
         sige_resist_erfc.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_SiGe_resist_erfc PROPERTY LABELS
     SiGe_resist_erfc nightly)
ENDIF()

####################################################
# Simple Si(1-x)Ge(x) resistor example with mgauss mf
SET(testName charon_mp_SiGe_resist_mgauss)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    sige_resist_mgauss_msh.jou 
    sige_resist_mgauss_msh.exo.4.0
    sige_resist_mgauss_msh.exo.4.1
    sige_resist_mgauss_msh.exo.4.2
    sige_resist_mgauss_msh.exo.4.3
    sige_resist_mgauss.nlp.inp
    sige_resist_mgauss.nlp.exodiff
    sige_resist_mgauss.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=sige_resist_mgauss.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=sige_resist_mgauss.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "sige_resist_mgauss.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         sige_resist_mgauss.nlp.exodiff
         sige_resist_mgauss.nlp.exo
         sige_resist_mgauss.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_SiGe_resist_mgauss PROPERTY LABELS
     SiGe_resist_mgauss nightly)
ENDIF()

####################################################
# Simple Si(1-x)Ge(x) resistor example with halo mf
SET(testName charon_mp_SiGe_resist_halo)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    sige_resist_halo_msh.jou 
    sige_resist_halo_msh.exo.4.0
    sige_resist_halo_msh.exo.4.1
    sige_resist_halo_msh.exo.4.2
    sige_resist_halo_msh.exo.4.3
    sige_resist_halo.nlp.inp
    sige_resist_halo.nlp.exodiff
    sige_resist_halo.nlp.gold.exo
    
  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i=sige_resist_halo.nlp.inp
    
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=sige_resist_halo.nlp.inp.xml
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    
  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "sige_resist_halo.nlp.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         sige_resist_halo.nlp.exodiff
         sige_resist_halo.nlp.exo
         sige_resist_halo.nlp.gold.exo
    NUM_MPI_PROCS 1
   
  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_SiGe_resist_halo PROPERTY LABELS
     SiGe_resist_halo nightly)
ENDIF()




