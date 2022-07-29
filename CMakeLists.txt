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

IF (EXISTS ${CMAKE_SOURCE_DIR}/Interpreter.cmake)
   INCLUDE (${CMAKE_SOURCE_DIR}/Interpreter.cmake)
ENDIF()

IF (EXISTS ${CMAKE_SOURCE_DIR}/pyMesh.cmake)
   INCLUDE (${CMAKE_SOURCE_DIR}/pyMesh.cmake)
ENDIF()

IF (EXISTS ${CMAKE_SOURCE_DIR}/DakotaDrivers.cmake)
   INCLUDE (${CMAKE_SOURCE_DIR}/DakotaDrivers.cmake)
ENDIF()

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

# Don't tell make to ignore errors (add "-i" to invocation). Doing so
# will cause ctest to proceed with running tests even when the build
# fails. I have NO idea why this isn't set by default. Perhaps
# something in Tribts/Trilinos is setting it to WARN, which is what it
# looks like it's being set to when traced via the debugger.
CMAKE_POLICY(SET CMP0061 NEW)
