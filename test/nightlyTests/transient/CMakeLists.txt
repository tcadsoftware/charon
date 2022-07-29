INCLUDE(TribitsAddAdvancedTest)

SET(testName charon_mp_pndiode.dd.lin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.dd.lin.fixstep.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    result_dd_lin_fixstep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.dd.lin.fixstep.inp

  TEST_2 CMND rm ARGS -f current_time_dd_lin.csv
    PASS_ANY

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.dd.lin.fixstep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_dd_lin_fixstep.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         result_dd_lin_fixstep.exo
         result_dd_lin_fixstep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.dd.lin.fixstep PROPERTY LABELS
    pndiode.dd.lin.fixstep nightly)
ENDIF()

SET(testName charon_mp_pndiode.sg.lin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.fixstep.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    result_sg_lin_fixstep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.lin.fixstep.inp

  TEST_2 CMND rm ARGS -f current_time_sg_lin.csv
    PASS_ANY

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.lin.fixstep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_sg_lin_fixstep.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         result_sg_lin_fixstep.exo
         result_sg_lin_fixstep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.fixstep PROPERTY LABELS
    pndiode.sg.lin.fixstep nightly)
ENDIF()

SET(testName charon_mp_pndiode.sg.lin.varstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.lin.varstep.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    result_sg_lin_varstep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.lin.varstep.inp

  TEST_2 CMND rm ARGS -f current_time_sg_lin_varstep.csv
    PASS_ANY

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.lin.varstep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_sg_lin_varstep.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         result_sg_lin_varstep.exo
         result_sg_lin_varstep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.lin.varstep PROPERTY LABELS
    pndiode.sg.lin.varstep nightly)
ENDIF()

##---------------------------------------------------------
##---------------------------------------------------------
SET(testName charon_mp_pndiode.sg.sin.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.sin.fixstep.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    result_sg_sin_fixstep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.sin.fixstep.inp

  TEST_2 CMND rm ARGS -f current_time_sg_sin.csv
    PASS_ANY

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.sin.fixstep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_sg_sin_fixstep.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         result_sg_sin_fixstep.exo
         result_sg_sin_fixstep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.sin.fixstep PROPERTY LABELS
    pndiode.sg.sin.fixstep  nightly)
ENDIF()

##---------------------------------------------------------
##---------------------------------------------------------
SET(testName charon_mp_pndiode.sg.tra.fixstep)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.tra.fixstep.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    result_sg_tra_fixstep.gold.exo

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.tra.fixstep.inp

  TEST_2 CMND rm ARGS -f current_time_sg_tra.csv
    PASS_ANY

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.tra.fixstep.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "result_sg_tra_fixstep.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_5 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         result_sg_tra_fixstep.exo
         result_sg_tra_fixstep.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.tra.fixstep PROPERTY LABELS
    pndiode.sg.tra.fixstep  nightly debugexclude)
ENDIF()

##---------------------------------------------------------
##---------------------------------------------------------
SET(testName charon_mp_pndiode.sg.trans.opt)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  CATEGORIES NIGHTLY
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    pndiode.sg.equ.inp
    pndiode.sg.-1V.inp
    pndiode.sg.trans.opt.inp
    potential_pndiode.exo.4.0
    potential_pndiode.exo.4.1
    potential_pndiode.exo.4.2
    potential_pndiode.exo.4.3
    pndiode.exodiff
    pndiode.sg.trans.opt.gold.exo
    timeDepG.txt 
    space1DFactor.txt
    space2DFactor.txt

  TEST_1 CMND rm ARGS -f current_time_trans_opt.csv
    PASS_ANY

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.equ.inp

  TEST_3 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.equ.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_4 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.-1V.inp

  TEST_5 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.-1V.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_6 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i pndiode.sg.trans.opt.inp

  TEST_7 EXEC charon_mp NOEXEPREFIX DIRECTORY ${PACKAGE_BINARY_DIR}
    ARGS  --i=pndiode.sg.trans.opt.inp.xml --current
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Charon run completed."

  TEST_8 EXEC ${SEACAS_BINARY_DIR}/applications/epu/epu NOEXEPREFIX NOEXESUFFIX
    ARGS  "-auto" "pndiode.sg.trans.opt.exo.4."
    NUM_MPI_PROCS 1
    
  TEST_9 EXEC ${SEACAS_BINARY_DIR}/applications/exodiff/exodiff NOEXEPREFIX
    NOEXESUFFIX
    ARGS -f
         pndiode.exodiff
         pndiode.sg.trans.opt.exo
         pndiode.sg.trans.opt.gold.exo
    NUM_MPI_PROCS 1

  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_pndiode.sg.trans.opt PROPERTY LABELS
    pndiode.sg.trans.opt  nightly)
ENDIF()
