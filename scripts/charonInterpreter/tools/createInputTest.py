#! /usr/bin/env python3
from __future__ import print_function

import sys

boilerPlate = """
####################################################
# Interpreter test
SET(testName charon_mp_<testName>)
TRIBITS_ADD_ADVANCED_TEST(
  ${testName}
  OVERALL_WORKING_DIRECTORY TEST_NAME

  TEST_0 COPY_FILES_TO_TEST_DIR
    <inputFileName>
    <inputFileGold>

  TEST_1 CMND ${PACKAGE_BINARY_DIR}/charonInterpreter 
    ARGS -i <inputFileName>

  TEST_2 CMND ${PACKAGE_BINARY_DIR}/interpreter/charonInterpreter/tools/compareParameterLists.py
    ARGS <inputFileName>.xml <inputFileGold>
    PASS_REGULAR_EXPRESSION \"Files are the same.\"


  ADDED_TEST_NAME_OUT ${testName}_CREATED
  )
IF (${testName}_CREATED)
  SET_PROPERTY(TEST Charon_charon_mp_<testName> PROPERTY LABELS
    <testName> inputVerification)
ENDIF()
"""

if len(sys.argv) < 4:
    print("usage: createInputTest.py  testName inputFileName inputFileGold")
    sys.exit(1)

testName = sys.argv[1]
inputFileName = sys.argv[2]
inputFileGold = sys.argv[3]

boilerPlate = boilerPlate.replace("<testName>",testName)
boilerPlate = boilerPlate.replace("<inputFileName>",inputFileName)
boilerPlate = boilerPlate.replace("<inputFileGold>",inputFileGold)

print(boilerPlate)

testFile = open("CMakeLists.txt","a+")
testFile.write("\n\n"+boilerPlate)
