INCLUDE(TribitsAddAdvancedTest)

####################################################
#  PN 1D diode HB SGCVFEM Laplace 1 tone test
####################################################

SET(testName charon_mp.hb-sgcvfem-laplace.pndiode)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.hb-sgcvfem-laplace.inp
    gold_pndiode_hb-sgcvfem-laplace.exo
    pndiode.hb-sgcvfem-laplace.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.hb-sgcvfem-laplace.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-sgcvfem-laplace.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-sgcvfem-laplace.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode_hb-sgcvfem-laplace.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f pndiode.hb-sgcvfem-laplace.exodiff
         pndiode_hb-sgcvfem-laplace.exo
         gold_pndiode_hb-sgcvfem-laplace.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-sgcvfem-laplace.pndiode  PROPERTY LABELS
    hb-sgcvfem-laplace.pndiode nightly)
ENDIF()



####################################################
#  PN 1D diode HB SGCVFEM Drift Diffusion 1 tone test
####################################################

SET(testName charon_mp.hb-sgcvfem-dd.pndiode)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    result_pndiode_dd.exo.4.0
    result_pndiode_dd.exo.4.1
    result_pndiode_dd.exo.4.2
    result_pndiode_dd.exo.4.3
    pndiode.hb-sgcvfem-dd.inp
    gold_pndiode_hb-sgcvfem-dd.exo
    pndiode.hb-sgcvfem-dd.exodiff

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter
    ARGS -i pndiode.hb-sgcvfem-dd.inp

  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.hb-sgcvfem-dd.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."
    OUTPUT_FILE pndiode.hb-sgcvfem-dd.inp.xml.out

  TEST_3 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode_hb-sgcvfem-dd.exo.4."
    NUM_MPI_PROCS 1

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS
         -f pndiode.hb-sgcvfem-dd.exodiff
         pndiode_hb-sgcvfem-dd.exo
         gold_pndiode_hb-sgcvfem-dd.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )

IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp.hb-sgcvfem-dd.pndiode  PROPERTY LABELS
    hb-sgcvfem-dd.pndiode nightly)
ENDIF()



####################################################
#  PN 1D diode HB DD 1 tone rectifying test 
####################################################
#
#SET(testName charon_mp.hb-dd.pndiode.ac1-rectifying)
#TRIBITS_ADD_ADVANCED_TEST(
#  ${testName}
#  CATEGORIES NIGHTLY
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
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb rectifying diode test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac1-rectifying.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac1-rectifying.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac1-rectifying.exo.4."
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
#  CATEGORIES NIGHTLY
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
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb rectifying diode test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac1-linear-LS.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac1-linear-LS.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac1-linear-LS.exo.4."
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
#  CATEGORIES NIGHTLY
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
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run dd problem, resulting with input exodus for hb-dd problem
#  TEST_2 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.dd.equ.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#
#  # run the hb diode envelope test problem
#  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
#    ARGS  --i=pndiode.hb-dd.ac2-envelope.yaml --current
#    NUM_MPI_PROCS 4
#    PASS_REGULAR_EXPRESSION "Charon run completed."
#    OUTPUT_FILE pndiode.hb-dd.ac2-envelope.current.output
#
#  # epu auto the result to prepare for exodiff
#  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
#    ARGS  "-auto" "result_pndiode_hb-dd_ac2-envelope.exo.4."
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




