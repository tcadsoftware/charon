#!/bin/bash -xe
#######################################################################
# This script should be pasted into the "Execute shell" portion of the
# "Build" page on a Jenkins charon test item. Right now it supports
# the standard RHEL 6 and 7, DBG or OPT, builds when appropriate
# variables below are modified accordingly.

# NOTE: You need to modify this script if you modify any of the
# jenkins items that it is associated with. At this time that
# includes:
#	charon_rhel7_nightly_dbg
#	charon_rhel7_nightly_opt
#######################################################################

# Tests to run
#
# Valid values:
#    n - nightly
#    h - heavy
TYPE="n"

# The build flag
#
# Valid values:
#	DBG
#	OPT
export OPTARG="<FILL THIS IN>"

# The opts file arguments for the python build script
#
# Valid values:
#	rhel7.opts
OPTSFILEARG="-f <FILL THIS IN>"

# The OS string
#
# Valid values:
#       RHEL7
#       RHEL7 Heavy
SITEDESCARG="JBF <FILL THIS IN>"

export TRIBITS_BASE_DIR=${WORKSPACE}/TriBITS

source /etc/profile.d/modules.sh
source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-cmake

module load sems-gcc/6.1.0
module load sems-openmpi/1.10.1
module load sems-hdf5/1.8.12/base
module load sems-netcdf/4.4.1/exo
module load sems-boost/1.59.0/base
module unload sems-python/2.7.9
module load sems-python

case ${TYPE} in

  "n")
      LSTRING="nightly"
      TESTTRACK="Nightly"
      ;;

  "h")
      LSTRING="^heavy$"
      TESTTRACK="Heavy"
      ;;

  *)
      echo "Unknown test type!"
      exit 1
      ;;

esac
    
DBGEXCLUDE=""
if [ "${OPTARG}" = "DBG" ]
then
  DBGEXCLUDE="-LE debugexclude"
fi

ctest \
  -L ${LSTRING} \
  --extra-verbose \
  -j 16 \
  ${DBGEXCLUDE} \
  -S "${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/ctest_regression.cmake" \
  -DDEBUGLEVEL:INT=8 \
  -DTYPE:STRING=${OPTARG} \
  -DDISTRIB:STRING="Jenkins" \
  -DCOMPILER:STRING="MPI_GCC_6.1.x" \
  -DPROCESSORCOUNT:INT=16 \
  -DBASETESTDIR:STRING="${WORKSPACE}" \
  -DBSCRIPTARGS:STRING="${OPTSFILEARG}" \
  -DSITEDESC:STRING="${SITEDESCARG}" \
  -DTESTTRACK=${TESTTRACK}
