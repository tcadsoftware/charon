INCLUDE(TribitsAddAdvancedTest)

####################################################
# Kimpton TID test
SET(testName charon_mp_Kimpton_nmoscap_FEM)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc_msh.jou
    nmosc_msh.gen.4.0
    nmosc_msh.gen.4.1
    nmosc_msh.gen.4.2
    nmosc_msh.gen.4.3
    nmosc.nlp.FEM.xml
    nmosc.nlp.FEM.gold.exo
    nmosc.nlp.FEM.exodiff
    nmosc.dd.FEM.xml
    nmosc.dd.FEM.gold.exo
    nmosc.dd.FEM.exodiff
    nmosc.dd.sweep.FEM.xml
    nmosc.dd.sweep.FEM.gold.exo
    nmosc.dd.sweep.FEM.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.nlp.FEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.nlp.FEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.nlp.FEM.exodiff
         nmosc.nlp.FEM.exo
         nmosc.nlp.FEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.FEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.FEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.FEM.exodiff
         nmosc.dd.FEM.exo
         nmosc.dd.FEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.sweep.FEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.sweep.FEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.sweep.FEM.exodiff
         nmosc.dd.sweep.FEM.exo
         nmosc.dd.sweep.FEM.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_FEM PROPERTY LABELS
     Kimpton_nmoscap_FEM nightly)
ENDIF()


####################################################
# Kimpton TID test
SET(testName charon_mp_Kimpton_nmoscap_SGCVFEM)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc_msh.jou
    nmosc_msh.gen.4.0
    nmosc_msh.gen.4.1
    nmosc_msh.gen.4.2
    nmosc_msh.gen.4.3
    nmosc.nlp.SGCVFEM.xml
    nmosc.nlp.SGCVFEM.gold.exo
    nmosc.nlp.SGCVFEM.exodiff
    nmosc.dd.SGCVFEM.xml
    nmosc.dd.SGCVFEM.gold.exo
    nmosc.dd.SGCVFEM.exodiff
    nmosc.dd.sweep.SGCVFEM.xml
    nmosc.dd.sweep.SGCVFEM.gold.exo
    nmosc.dd.sweep.SGCVFEM.exodiff

  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.nlp.SGCVFEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_2 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.nlp.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.nlp.SGCVFEM.exodiff
         nmosc.nlp.SGCVFEM.exo
         nmosc.nlp.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_4 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.SGCVFEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

    # test the correctness of the exodus output
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_6 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.SGCVFEM.exodiff
         nmosc.dd.SGCVFEM.exo
         nmosc.dd.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1

  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.sweep.SGCVFEM.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.sweep.SGCVFEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.sweep.SGCVFEM.exodiff
         nmosc.dd.sweep.SGCVFEM.exo
         nmosc.dd.sweep.SGCVFEM.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_SGCVFEM PROPERTY LABELS
     Kimpton_nmoscap_SGCVFEM nightly)
ENDIF()


####################################################
# Interpreter test
SET(testName charon_mp_Kimpton_nmoscap_FEM_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc.nlp.FEM.inp
    nmosc.nlp.FEM.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.nlp.FEM.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.nlp.FEM.inp.xml nmosc.nlp.FEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  TEST_3 COPY_FILES_TO_TEST_DIR
    nmosc.dd.FEM.inp
    nmosc.dd.FEM.inp.gold.xml

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.FEM.inp

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.dd.FEM.inp.xml nmosc.dd.FEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  TEST_6 COPY_FILES_TO_TEST_DIR
    nmosc.dd.sweep.FEM.inp
    nmosc.dd.sweep.FEM.inp.gold.xml

  TEST_7 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.sweep.FEM.inp

  TEST_8 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.dd.sweep.FEM.inp.xml nmosc.dd.sweep.FEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_FEM_input PROPERTY LABELS
     Kimpton_nmoscap_FEM inputVerification)
ENDIF()


####################################################
# Interpreter test
SET(testName charon_mp_Kimpton_nmoscap_SGCVFEM_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc.nlp.SGCVFEM.inp
    nmosc.nlp.SGCVFEM.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.nlp.SGCVFEM.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.nlp.SGCVFEM.inp.xml nmosc.nlp.SGCVFEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  TEST_3 COPY_FILES_TO_TEST_DIR
    nmosc.dd.SGCVFEM.inp
    nmosc.dd.SGCVFEM.inp.gold.xml

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.SGCVFEM.inp

  TEST_5 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.dd.SGCVFEM.inp.xml nmosc.dd.SGCVFEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  TEST_6 COPY_FILES_TO_TEST_DIR
    nmosc.dd.sweep.SGCVFEM.inp
    nmosc.dd.sweep.SGCVFEM.inp.gold.xml

  TEST_7 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.sweep.SGCVFEM.inp

  TEST_8 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS nmosc.dd.sweep.SGCVFEM.inp.xml nmosc.dd.sweep.SGCVFEM.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_SGCVFEM_input PROPERTY LABELS
     Kimpton_nmoscap_SGCVFEM inputVerification)
ENDIF()



####################################################
# Kimpton TID test
SET(testName charon_mp_Kimpton_nmoscap_FEM_TDSweep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc_msh.jou
    nmosc_msh.gen.4.0
    nmosc_msh.gen.4.1
    nmosc_msh.gen.4.2
    nmosc_msh.gen.4.3
    nmosc.nlp.FEM.inp
    nmosc.nlp.FEM.gold.exo
    nmosc.nlp.FEM.exodiff
    nmosc.dd.TDSweep.inp
    nmosc.dd.TDSweep.gold.exo
    nmosc.dd.TDSweep.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.nlp.FEM.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.nlp.FEM.inp.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.nlp.FEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.nlp.FEM.exodiff
         nmosc.nlp.FEM.exo
         nmosc.nlp.FEM.gold.exo
    NUM_MPI_PROCS 1
 
  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.TDSweep.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.TDSweep.inp.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.TDSweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.TDSweep.exodiff
         nmosc.dd.TDSweep.exo
         nmosc.dd.TDSweep.gold.exo
    NUM_MPI_PROCS 1


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_FEM_TDSweep PROPERTY LABELS
     Kimpton_nmoscap_FEM_TDSweep nightly)
ENDIF()



####################################################
# Kimpton TID test
SET(testName charon_mp_Kimpton_nmoscap_FEM_ICCSweep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    nmosc_msh.jou
    nmosc_msh.gen.4.0
    nmosc_msh.gen.4.1
    nmosc_msh.gen.4.2
    nmosc_msh.gen.4.3
    nmosc.nlp.FEM.inp
    nmosc.nlp.FEM.gold.exo
    nmosc.nlp.FEM.exodiff
    nmosc.dd.ICCSweep.inp
    nmosc.dd.ICCSweep.gold.exo
    nmosc.dd.ICCSweep.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.nlp.FEM.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.nlp.FEM.inp.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.nlp.FEM.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.nlp.FEM.exodiff
         nmosc.nlp.FEM.exo
         nmosc.nlp.FEM.gold.exo
    NUM_MPI_PROCS 1
 
  TEST_5 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i nmosc.dd.ICCSweep.inp

  TEST_6 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=nmosc.dd.ICCSweep.inp.xml 
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # test the correctness of the exodus output
  TEST_7 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "nmosc.dd.ICCSweep.exo.4."
    NUM_MPI_PROCS 1

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         nmosc.dd.ICCSweep.exodiff
         nmosc.dd.ICCSweep.exo
         nmosc.dd.ICCSweep.gold.exo
    NUM_MPI_PROCS 1


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_Kimpton_nmoscap_FEM_ICCSweep PROPERTY LABELS
     Kimpton_nmoscap_FEM_ICCSweep nightly)
ENDIF()
