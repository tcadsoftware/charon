# To be safe, define your minimum CMake version.  This may be newer than the
# minimum version required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0 FATAL_ERROR)

# Set the project name as a variable at the very beginning before including
# anything else.  This happens in a separate file so CTest scripts can use it.
INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)

IF(tcad-charon_ENABLE_GCOVR AND CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  INCLUDE(${CMAKE_SOURCE_DIR}/CodeCoverage.cmake)
ENDIF()

IF(NOT IS_DIRECTORY src/)
       execute_process (COMMAND mkdir src)
ENDIF()
IF(NOT IS_DIRECTORY src/interpreter)
       execute_process (COMMAND mkdir src/interpreter)
ENDIF()
IF(IS_DIRECTORY src/interpreter/parsers)
        execute_process (COMMAND rm -fr src/interpreter/parsers)
        execute_process (COMMAND rm -fr src/interpreter/modifiers)
ENDIF()
execute_process (COMMAND rm -fr src/interpreter )
execute_process (COMMAND rsync -av  --exclude generateParserDocumentation --exclude parsers --exclude modifiers --progress ${CMAKE_SOURCE_DIR}/scripts/charonInterpreter  src/interpreter --delete)
execute_process (WORKING_DIRECTORY src/interpreter/charonInterpreter/parseGenerator COMMAND python3 generateInterpreter.py --verbosity=20)
execute_process (COMMAND find ./src/interpreter \( -name "*pyc" \) -exec rm {} \;)
execute_process (COMMAND python3 --version)
execute_process (COMMAND mv src/interpreter/charonInterpreter/parseGenerator parseGenerator)
execute_process (WORKING_DIRECTORY src/interpreter/charonInterpreter/ COMMAND python3 -m compileall ./)
execute_process (WORKING_DIRECTORY src/ COMMAND chmod -R g+rX ./)
execute_process (WORKING_DIRECTORY src/ COMMAND ln -s interpreter/charonInterpreter/charonInterpreter.py charonInterpreter )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -s interpreter/charonInterpreter/tools/compareParameterLists.py compareParameterLists )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -s interpreter/charonInterpreter/tools/createInputTest.py createInputTest )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -s interpreter/charonInterpreter/tools/dagify.py dagify )
execute_process (WORKING_DIRECTORY src/ COMMAND ln -s interpreter/charonInterpreter/tools/xmlToLCM.py xmlToLCM )
execute_process (COMMAND mv parseGenerator src/interpreter/charonInterpreter/)


# CMake requires that you declare the CMake project in the top-level file and
# not in an include file.
PROJECT(${PROJECT_NAME} NONE)

# Pull in the TriBITS system and execute.
IF("$ENV{TRIBITS_BASE_DIR}" STREQUAL "")
  MESSAGE(FATAL_ERROR
    "You must set TRIBITS_BASE_DIR to point to your local tribits directory.")
ENDIF()
SET(${PROJECT_NAME}_TRIBITS_DIR "$ENV{TRIBITS_BASE_DIR}/tribits" CACHE STRING
  "TriBITS base directory.")
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

# Turn these on for DEBUG builds
IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  OPTION(tcad-charon_ENABLE_DEBUG "Enable Debugging" ON)
  OPTION(tcad-charon_ENABLE_DEBUG_SYMBOLS "Enable Debug symbols" ON)
ELSE()
  OPTION(tcad-charon_ENABLE_DEBUG "Enable Debugging" OFF)
  OPTION(tcad-charon_ENABLE_DEBUG_SYMBOLS "Enable Debug symbols" OFF)
ENDIF()

IF(tcad-charon_ENABLE_GCOVR AND CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  SET(COVERAGE_GCOVR_EXCLUDES ".*/Trilinos/.*")
  SETUP_TARGET_FOR_COVERAGE_GCOVR_HTML(
    NAME charon_coverage
    EXECUTABLE ctest -LE debugexclude -j 21
    DEPENDENCIES charon_mp.exe
    )
ENDIF()

# Do all of the processing for this TriBITS project.
TRIBITS_PROJECT()

