#
# Usage: ctest -S ctest_regression.cmake
#   -DTYPE:STRING=<OPT|DBG> -- Compiles using optimization (default is debug build)
#   -DMODEL:STRING=<Experimental|Nightly|Continuous> (default is Nightly)
#   -DCOMPILER:STRING=<compiler version information> (Example: -DCOMPILER=OMPI_GCC_4.4.7)
#   -DDISTRIB:STRING=<computer OS information>  (Example: -DDISTRIB:STRING=Mac_OSX_11.4.2)
#   -DPROCESSORCOUNT:INT=<number of processors for parallel make> (Example: -DPROCESSORCOUNT:INT=10)
#   -DBSCRIPTARGS:STRING=<arguments passed directly to the python build script>
#   -DSITEDESC:STRING=<name of the CTEST_SITE to use in CDash> (default: set based on host name)
#   -DOUOTESTS:BOOL=<TRUE|FALSE> (default is TRUE)
#   -DTBLDDIR:STRING=<name of the build directory for testing> (default: TEST.OPT, TEST.DBG or TEST.COV)
#   -DINSTALLDIR:PATH=<path to base installation directory> (default: not needed unless installing build)
#   -DBUILDONLY:BOOL=<TRUE|FALSE> (default is FALSE)
#   -DCOVERAGE:BOOL=<TRUE|FALSE> (default is FALSE)
#   -DTESTTRACK:STRING=<Name of CDash track> (no default)
#   -DDONOTUPDATE:BOOL=<TRUE|FALSE> (default is FALSE)
#   -DDONOTCLONE:BOOL=<TRUE|FALSE> (default is FALSE)
#   -DCIINVOCATION:BOOL=<TRUE|FALSE> (default is FALSE)
#
# This is a ctest script for running Charon nightly tests
SET(CTEST_PROJECT_NAME "Charon")
SET(CTEST_DROP_METHOD "https")
SET(CTEST_DROP_SITE "charon-cdash.sandia.gov")
SET(CTEST_DROP_LOCATION "/submit.php?project=Charon")

if(NOT CIINVOCATION)
  SET(CITEST FALSE)
else()
  SET(CITEST TRUE)
endif()

if(NOT OUOTESTS)
  SET(OUOTESTS TRUE)
endif()

if (NOT COVERAGE)
  SET(COVERAGE FALSE)
endif()

if(NOT DONOTUPDATE)
  SET(DONOTUPDATE FALSE)
endif()
if(NOT DONOTCLONE)
  SET(DONOTCLONE FALSE)
endif()

SET(TCAD_CHARON_REPO "git@cee-gitlab.sandia.gov:Charon/tcad-charon")
SET(SRC_REPO "git@cee-gitlab.sandia.gov:Charon/src")
SET(NIGHTLY_TESTS_REPO "git@cee-gitlab.sandia.gov:Charon/nightlyTests")
if (OUOTESTS)
  SET(NIGHTLY_TESTS_OUO_REPO
    "git@cee-gitlab.sandia.gov:Charon/nightlyTestsOUO")
endif()
SET(TRILINOS_REPO "git@cee-gitlab.sandia.gov:Charon/Trilinos")

# Set default values for variables
if(NOT DEBUGLEVEL)
  SET(DEBUGLEVEL 0)
endif()

# Test for directories
if ("${BASETESTDIR}" STREQUAL "")
  MESSAGE(FATAL_ERROR
    "ERROR: Must provide a BASETESTDIR option for base test directory")
else()
  if (NOT EXISTS "${BASETESTDIR}")
    MESSAGE(FATAL_ERROR
      "ERROR: Specified base directory ${BASETESTDIR} doesn't exist")
  endif()
endif()

SET(TCAD_CHARON_SOURCE_DIRECTORY "${BASETESTDIR}/tcad-charon")
SET(SRC_SOURCE_DIRECTORY "${TCAD_CHARON_SOURCE_DIRECTORY}/src")
SET(NIGHTLY_TESTS_SOURCE_DIRECTORY
  "${TCAD_CHARON_SOURCE_DIRECTORY}/test/nightlyTests")
if(OUOTESTS)
  SET(NIGHTLY_TESTS_OUO_SOURCE_DIRECTORY
    "${TCAD_CHARON_SOURCE_DIRECTORY}/test/nightlyTestsOUO")
endif()
SET(TRILINOS_SOURCE_DIRECTORY "${TCAD_CHARON_SOURCE_DIRECTORY}/Trilinos")

SET(CTEST_SOURCE_DIRECTORY "${TCAD_CHARON_SOURCE_DIRECTORY}")

####################################################################
####################################################################
# Nothing below this needs to be modified for other Linux platforms
# with the possible exception of the host-specific settings
####################################################################
####################################################################
SET(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 524288)
SET(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 1048576)
if("${TBLDDIR}" STREQUAL "")
  if(${TYPE} MATCHES OPT)
    SET(CTEST_BINARY_DIRECTORY "${BASETESTDIR}/TEST.OPT")
  elseif(${TYPE} MATCHES COV)
    SET(CTEST_BINARY_DIRECTORY "${BASETESTDIR}/TEST.COV")
  else()
    SET(CTEST_BINARY_DIRECTORY "${BASETESTDIR}/TEST.DBG")
  endif()
else()
  SET(CTEST_BINARY_DIRECTORY "${BASETESTDIR}/${TBLDDIR}")
endif()

# Get some information about the platform to name this test run
FIND_PROGRAM(UNAME NAMES uname)
MACRO(getuname name flag)
  EXECUTE_PROCESS(COMMAND "${UNAME}" "${flag}"
    OUTPUT_VARIABLE "${name}"
    OUTPUT_STRIP_TRAILING_WHITESPACE)
ENDMACRO()

getuname(osname -s)
getuname(proc -m)

if(DEBUGLEVEL GREATER 2)
  MESSAGE("DEBUG: osname: ${osname}")
  MESSAGE("DEBUG: proc: ${proc}")
endif()

FIND_PROGRAM(WHOAMI NAMES whoami)
EXECUTE_PROCESS(COMMAND "${WHOAMI}"
  OUTPUT_VARIABLE username
  OUTPUT_STRIP_TRAILING_WHITESPACE)

FIND_PROGRAM(HNAME NAMES hostname)
EXECUTE_PROCESS(COMMAND "${HNAME}" "-f"
  OUTPUT_VARIABLE myhost
  OUTPUT_STRIP_TRAILING_WHITESPACE)

FIND_PROGRAM(CTEST_GIT_COMMAND NAMES git)

if(DEBUGLEVEL GREATER 2)
  MESSAGE("DEBUG: username: ${username}")
  MESSAGE("DEBUG: myhost: ${myhost}")
endif()

if("${SITEDESC}" STREQUAL "")
  SET(CTEST_SITE "${myhost}")
else()
  SET(CTEST_SITE "${SITEDESC}")
endif()

#########################################
# Any host specific settings go here
#########################################
if(NOT CITEST)
  if(${TYPE} MATCHES OPT)
    SET(CTEST_BUILD_NAME "${osname}-${proc}-${DISTRIB}-${COMPILER}-OPT")
    SET(BLDARG "-b r")
  elseif(${TYPE} MATCHES COV)
    SET(CTEST_BUILD_NAME "${osname}-${proc}-${DISTRIB}-${COMPILER}-COV")
    SET(BLDARG "-b d")
  else()
    SET(CTEST_BUILD_NAME "${osname}-${proc}-${DISTRIB}-${COMPILER}-DBG")
    SET(BLDARG "-b d")
  endif()
else()
  # In the case of a CI test just use what comes in via the COMPILER
  # argument
  SET(CTEST_BUILD_NAME "${COMPILER}")
  if(${TYPE} MATCHES OPT)
    SET(BLDARG "-b r")
  else()
    SET(BLDARG "-b d")
  endif()
endif()
MESSAGE("Build name: ${CTEST_BUILD_NAME}")

FIND_PROGRAM(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
FIND_PROGRAM(CTEST_COVERAGE_COMMAND NAMES gcov)

# This is supposed to invoke make with "-j <processor count>" but
# somewhere during the build process this is frustratingly reset. This
# makes the testing process take a large amount of time because of the
# serial compile during part of the build.
if (NOT PROCESSORCOUNT)
  SET(PROCESSORCOUNT 1)
endif()

SET(CTEST_BUILD_FLAGS "-j${PROCESSORCOUNT}")

# Don't do updates if this is running as the user "jenkins"
if ("${username}" STREQUAL "jenkins")
  SET(DONOTUPDATE TRUE)
  SET(DONOTCLONE TRUE)
endif()

# if this is specified then don't clone
if (NOT DONOTCLONE)
  # Clone the Charon repositories.
  if(NOT EXISTS "${TCAD_CHARON_SOURCE_DIRECTORY}")
    MESSAGE("Getting tcad-charon repository...")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" clone "${TCAD_CHARON_REPO}"
      "${TCAD_CHARON_SOURCE_DIRECTORY}" RESULT_VARIABLE res_tcad_charon_co)
    MESSAGE("Result of tcad-charon checkout: ${res_tcad_charon_co}")
  endif()
  if(NOT EXISTS "${SRC_SOURCE_DIRECTORY}")
    MESSAGE("Getting src repository...")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" clone "${SRC_REPO}"
      "${SRC_SOURCE_DIRECTORY}" RESULT_VARIABLE res_src_co)
    MESSAGE("Result of src checkout: ${res_src_co}")
  endif()
  if(NOT EXISTS "${NIGHTLY_TESTS_SOURCE_DIRECTORY}")
    MESSAGE("Getting nightlyTests repository...")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" clone "${NIGHTLY_TESTS_REPO}"
      "${NIGHTLY_TESTS_SOURCE_DIRECTORY}" RESULT_VARIABLE res_nightly_tests_co)
    MESSAGE("Result of nightlyTests checkout: ${res_nightly_tests_co}")
  endif()
  if (OUOTESTS)
    if(NOT EXISTS "${NIGHTLY_TESTS_OUO_SOURCE_DIRECTORY}")
      MESSAGE("Getting nightlyTestsOUO repository...")
      EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" clone
        "${NIGHTLY_TESTS_OUO_REPO}" "${NIGHTLY_TESTS_OUO_SOURCE_DIRECTORY}"
        RESULT_VARIABLE res_nightly_tests_ouo_co)
      MESSAGE("Result of nightlyTestsOUO checkout: ${res_nightly_tests_ouo_co}")
    endif()
  endif()
  if(NOT EXISTS "${TRILINOS_SOURCE_DIRECTORY}")
    MESSAGE("Getting initial Trilinos repository...")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" clone "${TRILINOS_REPO}"
      "${TRILINOS_SOURCE_DIRECTORY}" RESULT_VARIABLE res_trico)
    MESSAGE("Result of trilinos clone: ${res_trico}")
  endif()

  # Checkout the Trilinos branch, if requested
  if (NOT "${TRIBRANCH}" STREQUAL "")
    MESSAGE("Checking out Trilinos branch...")
    EXECUTE_PROCESS(COMMAND "${CTEST_GIT_COMMAND}" checkout "${TRIBRANCH}"
      WORKING_DIRECTORY "${TRILINOS_SOURCE_DIRECTORY}"
      RESULT_VARIABLE res_trico)
    MESSAGE("Result of trilinos checkout: ${res_trico}")
  endif()
endif()

# This needs to be set because of CTestConfig.cmake in Trilinos
SET(CMAKE_MODULE_PATH
  "${TRILINOS_SOURCE_DIRECTORY}/cmake/tribits/core/utils")

SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")
SET(CTEST_PROJECT_NAME "Charon")

SET(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
SET(CTEST_CMAKE_COMMAND
  "${TCAD_CHARON_SOURCE_DIRECTORY}/scripts/build/all/build_charon.py")

# Begin setting the options controlling the build of Charon
SET(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} ${CMAKEHOSTDEFS}")

# If the invocation requests an install prefix
if(NOT "${INSTALLDIR}" STREQUAL "")
  SET(CTEST_CONFIGURE_COMMAND "${CTEST_CMAKE_COMMAND} --reset-option=CMAKE_INSTALL_PREFIX:${INSTALLDIR}")
endif()

if(NOT "${BSCRIPTARGS}" STREQUAL "")
  SET(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${BSCRIPTARGS}")
endif()

SET(PYDBGLEVEL "--debug=${DEBUGLEVEL}")

SET(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} ${PYDBGLEVEL}")
SET(CTEST_CONFIGURE_COMMAND
  "${CTEST_CONFIGURE_COMMAND} ${BLDARG} ${TCAD_CHARON_SOURCE_DIRECTORY}")

# Set defaults if not specified
if("${MODEL}" STREQUAL "")
  SET(MODEL "Nightly")
endif()

# Ignore all warnings from compilation of Trilinos packages
SET(CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  "Trilinos/packages"
  )

# Some warnings to ignore when coupling to Xyce. This looks for the
# string "Xyce" in the description to indicate a Xyce coupled build
string(FIND "${SITEDESC}" "Xyce" XYCEBUILD)

if(DEBUGLEVEL GREATER 2)
  MESSAGE("DEBUG: xyce var: ${XYCEBUILD}")
endif()

if (${XYCEBUILD} GREATER_EQUAL 0)
   SET(CTEST_CUSTOM_WARNING_EXCEPTION
       ${CTEST_CUSTOM_WARNING_EXCEPTION}
       "N_.*\\.h:.*warning:"
       "N_.*\\.h:.*note:"
       "boost/mpi/config.hpp.*warning:.*OMPI_BUILD_CXX_BINDINGS.*redefined"
       "mpi.h.*note: this is the location of the previous definition"
      )
endif()

if(DEBUGLEVEL GREATER 2)
  MESSAGE("DEBUG: warning exceptions: ${CTEST_CUSTOM_WARNING_EXCEPTION}")
  MESSAGE("DEBUG: ctest model: ${MODEL}")
  MESSAGE("DEBUG: ctest TRACK: ${TESTTRACK}")
  MESSAGE("DEBUG: ctest username: ${username}")
  MESSAGE("DEBUG: DONOTUPDATE: ${DONOTUPDATE}")
endif()

# Start the testing process
if("${TESTTRACK}" STREQUAL "")
  if(DEBUGLEVEL GREATER 2)
    MESSAGE("DEBUG: starting default test track...")
  endif()
  ctest_start(${MODEL})
else()
  if(DEBUGLEVEL GREATER 2)
    MESSAGE("DEBUG: starting test track ${TESTTRACK}...")
  endif()
  ctest_start(${MODEL} TRACK ${TESTTRACK})
endif()

# Don't do any updates of repositories if this is jenkins. Jenkins
# itself does those
if (NOT DONOTUPDATE)
  ctest_update(SOURCE "${TCAD_CHARON_SOURCE_DIRECTORY}" RETURN_VALUE
    res_tcad_charon_up)
  ctest_update(SOURCE "${SRC_SOURCE_DIRECTORY}" RETURN_VALUE res_src_up)
  ctest_update(SOURCE "${NIGHTLY_TESTS_SOURCE_DIRECTORY}" RETURN_VALUE
    res_nightly_tests_up)
  if (OUOTESTS)
    ctest_update(SOURCE "${NIGHTLY_TESTS_OUO_SOURCE_DIRECTORY}" RETURN_VALUE
      res_nightly_tests_ouo_up)
  endif()
  ctest_update(SOURCE "${TRILINOS_SOURCE_DIRECTORY}" RETURN_VALUE
    res_trilinos_up)
endif()

if(DEBUGLEVEL GREATER 2)
  MESSAGE("DEBUG: ctest_custom_warning_exception value is:")
  MESSAGE("DEBUG: ${CTEST_CUSTOM_WARNING_EXCEPTION}")
endif()

ctest_configure()
ctest_build()
if(NOT BUILDONLY)
  ctest_test(PARALLEL_LEVEL ${PROCESSORCOUNT} RETURN_VALUE res_test)
  if (DEBUGLEVEL GREATER 2)
    MESSAGE("DEBUG: ctest returned exit code ${res_test}")
  endif()

  if(COVERAGE)

    # Don't want Trilinos code in coverage, only Charon code.
    SET(CTEST_CUSTOM_COVERAGE_EXCLUDE
      ${CTEST_CUSTOM_COVERAGE_EXCLUDE}
      "Trilinos/.*"
      )

    if(DEBUGLEVEL GREATER 2)
      MESSAGE("DEBUG: CTEST_COVERAGE_COMMAND: ${CTEST_COVERAGE_COMMAND}")
      MESSAGE("DEBUG: CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")
    endif()

    ctest_coverage(CAPTURE_CMAKE_ERROR res_var)
  endif()

  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 30)
  if (NOT res_test EQUAL 0)
    MESSAGE(FATAL_ERROR "Testing failure(s).")
  endif()
endif()
