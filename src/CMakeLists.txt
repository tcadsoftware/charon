INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

####################################
# A) Define the package
TRIBITS_PACKAGE(Charon)

####################################
# B) Set up package-specific options
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  CHARON_DEBUG
  "Build with debug code enabled."
  OFF )

# If this is a git repo we're working from the generate some
# additional version info from SHAs
IF (EXISTS ${CMAKE_SOURCE_DIR}/.git)
  SET(Charon_VERSION_STRING "R\"===(------------------------------------------------------------\n")

  TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}
    Charon_TcadCharon_VERSION_STRING)
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}tcad-charon:\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------\n${Charon_TcadCharon_VERSION_STRING}\n")

  TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/src
    Charon_Src_VERSION_STRING)
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}src:\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}----\n${Charon_Src_VERSION_STRING}\n")

  TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/test/nightlyTests
    Charon_NightlyTests_VERSION_STRING)
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}nightlyTests:\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}-------------\n${Charon_NightlyTests_VERSION_STRING}\n")

  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/nightlyTestsOUO)
    TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/test/nightlyTestsOUO
      Charon_NightlyTestsOUO_VERSION_STRING)
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}nightlyTestsOUO:\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}----------------\n${Charon_NightlyTestsOUO_VERSION_STRING}\n")
  ENDIF()

  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/heavyTests)
    TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/test/heavyTests
      Charon_HeavyTests_VERSION_STRING)
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}heavyTests:\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}-----------\n${Charon_HeavyTests_VERSION_STRING}\n")
  ENDIF()

  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/heavyTestsOUO)
    TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/test/heavyTestsOUO
      Charon_HeavyTestsOUO_VERSION_STRING)
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}heavyTestsOUO:\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}--------------\n${Charon_HeavyTestsOUO_VERSION_STRING}\n")
  ENDIF()

  IF (EXISTS ${CMAKE_SOURCE_DIR}/docs)
    TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/docs
      Charon_Docs_VERSION_STRING)
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}docs:\n")
    SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}-----\n${Charon_Docs_VERSION_STRING}\n")
  ENDIF()

  TRIBITS_GENERATE_SINGLE_REPO_VERSION_STRING(${CMAKE_SOURCE_DIR}/Trilinos
    Charon_Trilinos_VERSION_STRING)
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}Trilinos:\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}---------\n${Charon_Trilinos_VERSION_STRING}\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING}------------------------------------------------------------\n")
  SET(Charon_VERSION_STRING "${Charon_VERSION_STRING})===\"")
ELSE()
  SET(Charon_VERSION_STRING "R\"===(------------------------------------------------------------)===\"\n")
ENDIF()
FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/charonVersionString.txt"
  "${Charon_VERSION_STRING}")

# Explicit Template Instantiation
# force eti flag to mirror panzer's instead of trilinos default
# the tribits macro has been copied out and the default changed
##TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
  HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
  "Enable the use of explicit template instantiation (defaults to Panzer value)."
  ${Panzer_ENABLE_EXPLICIT_INSTANTIATION}
  )
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DECL_EXT "_decl.hpp")
SET(${PARENT_PACKAGE_NAME}_TEMPLATE_DEF_EXT "_impl.hpp")

####################################
# C) Add the library
INCLUDE(TribitsLibraryMacros)

#
# C.1) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# C.2) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

#################################################################
# add the option to couple to Xyce for cluster defect problems
#################################################################

##I dont need this  SET(ENABLE_XYCE_CLUSTER      "couple Charon to Xyce for cluster defect problems" )

IF (ENABLE_XYCE_CLUSTER)

  MESSAGE(STATUS "Attempting to link to Xyce for cluster calculations.")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/xyceCluster)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/xyceCluster)
  APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/xyceCluster/*.hpp)
  APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/xyceCluster/*.cpp)

ENDIF()
#################################################################
# add the option to couple to Xyce for cluster defect problems-END
#################################################################


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/bc_strategies)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/equation_sets)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/responses)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/solver)

# auto generate headers, add to dependencies and include directories
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/bc_strategies)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/equation_sets)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/evaluators)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/responses)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR}/solver)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/bc_strategies/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/equation_sets/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.hpp)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/solver/*.hpp)

APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/bc_strategies/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/equation_sets/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/evaluators/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/responses/*.cpp)
APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solver/*.cpp)

STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/Charon_Main.cpp;" "" SOURCES "${SOURCES}")
STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/xmlToYaml.cpp;" "" SOURCES "${SOURCES}")

#
# C.3) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  charon
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

####################################
# D) Add the executables and tests
INCLUDE(TribitsAddExecutable)

Tribits_ADD_EXECUTABLE(
  charon_mp
  SOURCES Charon_Main.cpp ${UNIT_TEST_DRIVER}
  NOEXEPREFIX
  INSTALLABLE
  )

Tribits_ADD_EXECUTABLE(
  xmlToYaml
  SOURCES xmlToYaml.cpp
  NOEXEPREFIX
  NOEXESUFFIX
  INSTALLABLE
  )

IF (${PACKAGE_NAME}_ENABLE_TESTS)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/test/unitTests
    ${CMAKE_BINARY_DIR}/test/unitTests)
  ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/test/nightlyTests
    ${CMAKE_BINARY_DIR}/test/nightlyTests)
  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/nightlyTestsOUO)
    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/test/nightlyTestsOUO
      ${CMAKE_BINARY_DIR}/test/nightlyTestsOUO)
  ENDIF()
  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/heavyTests)
    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/test/heavyTests
      ${CMAKE_BINARY_DIR}/test/heavyTests)
  ENDIF()
  IF (EXISTS ${CMAKE_SOURCE_DIR}/test/heavyTestsOUO)
    ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/test/heavyTestsOUO
      ${CMAKE_BINARY_DIR}/test/heavyTestsOUO)
  ENDIF()
ENDIF()

# Exclude files for source package.
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

####################################
# E) Do standard postprocessing
TRIBITS_PACKAGE_POSTPROCESS()
