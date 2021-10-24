INCLUDE(TribitsAddAdvancedTest)

####################################################
#  PN 1D diode HB DD 1 tone rectifying test 
####################################################
#
#SET(testName charon_mp.hb-dd.pndiode.ac1-rectifying)
#TRIBITS_ADD_ADVANCED_TEST(
#  ${testName}
#  OVERALL_WORKING_DIRECTORY TEST_NAME
#
#  TEST_0 COPY_FILES_TO_TEST_DIR
#    pndiode.exo.4.0
#    pndiode.exo.4.1
#    pndiode.exo.4.2
#    pndiode.exo.4.3
#    pndiode.nlp.yaml
#    pndiode.dd.yaml
#    pndiode.hb-dd.ac1-rectifying.yaml
#    gold_pndiode_hb-dd_ac1-rectifying.exo
#    pndiode.hb-dd.ac1-rectifying.exodiff
#
#  # run nlp problem, resulting with input exodus for dd problem
#  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.nlp.yaml
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb rectifying diode test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac1-rectifying.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac1-rectifying.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac1-rectifying.exo.${MPI_EXEC_MAX_NUMPROCS}."
#    NUM_MPI_PROCS 1
#
#  # compare results
#  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
#    NOEXESUFFIX
#    ARGS -f
#         pndiode.hb-dd.ac1-rectifying.exodiff
#         result_pndiode_hb-dd_ac1-rectifying.exo
#         gold_pndiode_hb-dd_ac1-rectifying.exo
#    NUM_MPI_PROCS 1
#
#  ADDED_TEST_NAME_OUT ${testName}_CREATED
#  )
#
#IF (${testName}_CREATED)
#  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.pndiode.ac1-rectifying  PROPERTY LABELS
#    hb-dd.pndiode.ac1-rectifying nightly)
#ENDIF()



####################################################
#  PN 1D diode HB DD 1 tone linear LS test
####################################################
#
#SET(testName charon_mp.hb-dd.pndiode.ac1-linear-LS)
#TRIBITS_ADD_ADVANCED_TEST(
#  ${testName}
#  OVERALL_WORKING_DIRECTORY TEST_NAME
#
#  TEST_0 COPY_FILES_TO_TEST_DIR
#    pndiode.exo.4.0
#    pndiode.exo.4.1
#    pndiode.exo.4.2
#    pndiode.exo.4.3
#    pndiode.nlp.yaml
#    pndiode.dd.yaml
#    pndiode.hb-dd.ac1-linear-LS.yaml
#    gold_pndiode_hb-dd_ac1-linear-LS.exo
#    pndiode.hb-dd.ac1-linear-LS.exodiff
#
#  # run nlp problem, resulting with input exodus for dd problem
#  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.nlp.yaml
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb rectifying diode test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac1-linear-LS.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac1-linear-LS.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac1-linear-LS.exo.${MPI_EXEC_MAX_NUMPROCS}."
#    NUM_MPI_PROCS 1
#
#  # compare results
#  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
#    NOEXESUFFIX
#    ARGS -f
#         pndiode.hb-dd.ac1-linear-LS.exodiff
#         result_pndiode_hb-dd_ac1-linear-LS.exo
#         gold_pndiode_hb-dd_ac1-linear-LS.exo
#    NUM_MPI_PROCS 1
#
#  ADDED_TEST_NAME_OUT ${testName}_CREATED
#  )
#
#IF (${testName}_CREATED)
#  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.pndiode.ac1-linear-LS  PROPERTY LABELS
#    hb-dd.pndiode.ac1-linear-LS nightly)
#ENDIF()



####################################################
#  PN 1D diode HB DD 2 tone envelope test
####################################################
#SET(testName charon_mp.hb-dd.pndiode.ac2-envelope)
#TRIBITS_ADD_ADVANCED_TEST(
#  ${testName}
#  OVERALL_WORKING_DIRECTORY TEST_NAME
#
#  TEST_0 COPY_FILES_TO_TEST_DIR
#    pndiode.exo.4.0
#    pndiode.exo.4.1
#    pndiode.exo.4.2
#    pndiode.exo.4.3
#    pndiode.nlp.yaml
#    pndiode.dd.equ.yaml
#    pndiode.hb-dd.ac2-envelope.yaml
#    gold_pndiode_hb-dd_ac2-envelope.exo
#    pndiode.hb-dd.ac2-envelope.exodiff
#
#  # run nlp problem, resulting with input exodus for dd problem
#  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.nlp.yaml
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.equ.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb diode envelope test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac2-envelope.yaml --current
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac2-envelope.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac2-envelope.exo.${MPI_EXEC_MAX_NUMPROCS}."
#    NUM_MPI_PROCS 1
#
#  # compare results
#  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
#    NOEXESUFFIX
#    ARGS -f
#         pndiode.hb-dd.ac2-envelope.exodiff
#         result_pndiode_hb-dd_ac2-envelope.exo
#         gold_pndiode_hb-dd_ac2-envelope.exo
#    NUM_MPI_PROCS 1
#
#  ADDED_TEST_NAME_OUT ${testName}_CREATED
#  )
#
#IF (${testName}_CREATED)
#  SET_PROPERTY(TEST Charon_charon_mp.hb-dd.pndiode.ac2-envelope  PROPERTY LABELS
#    hb-dd.pndiode.ac2-envelope nightly)
#ENDIF()




####################################################
#  PN 1D diode HB SGCVFEM Drift Diffusion 1 tone linear LS test
####################################################

SET(testName charon_mp.hb-sgcvfem-dd.pndiode.ac1-linear-SS)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.yaml
    pndiode.dd.yaml
    pndiode.hb-sgcvfem-dd.ac1-linear-SS.yaml
    gold_pndiode_hb-sgcvfem-dd_ac1-linear-SS.exo
    pndiode.hb-sgcvfem-dd.ac1-linear-SS.exodiff

  # run nlp problem, resulting with input exodus for dd problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # run dd problem, resulting with input exodus for hb-sgcvfem-dd problem
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # run the hb-sgcvfem-dd rectifying diode test problem
  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-sgcvfem-dd.ac1-linear-SS.yaml --current
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-sgcvfem-dd.ac1-linear-SS.current.output

  # epu auto the result to prepare for exodiff
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_hb-sgcvfem-dd_ac1-linear-SS.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.hb-sgcvfem-dd.ac1-linear-SS.exodiff
         result_pndiode_hb-sgcvfem-dd_ac1-linear-SS.exo
         gold_pndiode_hb-sgcvfem-dd_ac1-linear-SS.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-sgcvfem-dd.pndiode.ac1-linear-SS  PROPERTY LABELS
    hb-sgcvfem-dd.pndiode.ac1-linear-SS nightly)
ENDIF()




####################################################
#  PN 1D diode HB SGCVFEM Laplace 1 tone linear LS test
####################################################

SET(testName charon_mp.hb-sgcvfem-laplace.pndiode.ac1-linear-SS)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.yaml
    pndiode.hb-sgcvfem-laplace.ac1-linear-SS.yaml
    gold_pndiode_hb-sgcvfem-laplace_ac1-linear-SS.exo
    pndiode.hb-sgcvfem-laplace.ac1-linear-SS.exodiff

  # run nlp problem, resulting with input exodus for hb-sgcvfem-laplace problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # run the hb-sgcvfem-laplace pn diode test problem
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-sgcvfem-laplace.ac1-linear-SS.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-sgcvfem-laplace.ac1-linear-SS.output

  # epu auto the result to prepare for exodiff
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_hb-sgcvfem-laplace_ac1-linear-SS.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.hb-sgcvfem-laplace.ac1-linear-SS.exodiff
         result_pndiode_hb-sgcvfem-laplace_ac1-linear-SS.exo
         gold_pndiode_hb-sgcvfem-laplace_ac1-linear-SS.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-sgcvfem-laplace.pndiode.ac1-linear-SS  PROPERTY LABELS
    hb-sgcvfem-laplace.pndiode.ac1-linear-SS nightly)
ENDIF()



####################################################
#  PN 1D diode HB Laplace 1 tone linear LS test
####################################################

SET(testName charon_mp.hb-laplace.pndiode.ac1-linear-SS)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.exo.4.0
    pndiode.exo.4.1
    pndiode.exo.4.2
    pndiode.exo.4.3
    pndiode.nlp.yaml
    pndiode.hb-laplace.ac1-linear-SS.yaml
    gold_pndiode_hb-laplace_ac1-linear-SS.exo
    pndiode.hb-laplace.ac1-linear-SS.exodiff

  # run nlp problem, resulting with input exodus for hb-laplace problem
  TEST_1 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.nlp.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."

  # run the hb-laplace pn diode test problem
  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-laplace.ac1-linear-SS.yaml
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-laplace.ac1-linear-SS.output

  # epu auto the result to prepare for exodiff
  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_pndiode_hb-laplace_ac1-linear-SS.exo.${MPI_EXEC_MAX_NUMPROCS}."
    NUM_MPI_PROCS 1

  # compare results
  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.hb-laplace.ac1-linear-SS.exodiff
         result_pndiode_hb-laplace_ac1-linear-SS.exo
         gold_pndiode_hb-laplace_ac1-linear-SS.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-laplace.pndiode.ac1-linear-SS  PROPERTY LABELS
    hb-laplace.pndiode.ac1-linear-SS nightly)
ENDIF()




####################################################
# Interpreter test
SET(testName charon_mp_pndiode.hb-sgcvfem-dd.ac1-linear-SS_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.hb-sgcvfem-dd.ac1-linear-SS.inp
    pndiode.hb-sgcvfem-dd.ac1-linear-SS.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.hb-sgcvfem-dd.ac1-linear-SS.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.hb-sgcvfem-dd.ac1-linear-SS.inp.xml pndiode.hb-sgcvfem-dd.ac1-linear-SS.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.hb-sgcvfem-dd.ac1-linear-SS_input PROPERTY LABELS
    pndiode.hb-sgcvfem-dd.ac1-linear-SS_input inputVerification)
ENDIF()



####################################################
# Interpreter test
SET(testName charon_mp_pndiode.hb-sgcvfem-laplace.ac1-linear-SS_input)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.hb-sgcvfem-laplace.ac1-linear-SS.inp
    pndiode.hb-sgcvfem-laplace.ac1-linear-SS.inp.gold.xml

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.hb-sgcvfem-laplace.ac1-linear-SS.inp

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS pndiode.hb-sgcvfem-laplace.ac1-linear-SS.inp.xml pndiode.hb-sgcvfem-laplace.ac1-linear-SS.inp.gold.xml
    PASS_REGULAR_EXPRESSION "Files are the same."


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.hb-sgcvfem-laplace.ac1-linear-SS_input PROPERTY LABELS
    pndiode.hb-sgcvfem-laplace.ac1-linear-SS_input inputVerification)
ENDIF()
