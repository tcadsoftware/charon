#!/bin/bash -xe
#######################################################################
# This script should be pasted into the "Execute shell" portion of the
# "Build" page on a Jenkins charon test item. Right now it supports
# the standard RHEL 6 and 7, DBG or OPT, builds when appropriate
# variables below are modified accordingly.

# NOTE: You need to modify this script if you modify any of the
# jenkins items that it is associated with. At this time that
# includes:
#	charon_rhel6_nightly_dbg
#	charon_rhel6_nightly_opt
#	charon_rhel7_nightly_dbg
#	charon_rhel7_nightly_opt
#######################################################################

# The build flag
#
# Valid values:
#	DBG
#	OPT
export OPTARG="<FILL THIS IN>"

# The opts file arguments for the python build script
#
# Valid values:
#	rhel6.opts
#	rhel7.opts
OPTSFILEARG="-f <FILL THIS IN>"

# The OS string
#
# Valid values:
#	RHEL6
#	RHEL7
SITEDESCARG="JBF <FILL THIS IN>"

export TRIBITS_BASE_DIR=${WORKSPACE}/TriBITS

source /etc/profile.d/modules.sh
source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-cmake

module load sems-gcc/6.1.0
module load sems-openmpi/1.10.1
module load sems-python
module load sems-hdf5/1.8.12/base
module load sems-netcdf/4.4.1/exo
module load sems-boost/1.59.0/base

DBGEXCLUDE=""
if [ "${OPTARG}" = "DBG" ]
then
  DBGEXCLUDE="-LE debugexclude"
fi

ctest --extra-verbose -j 16 -L nightly ${DBGEXCLUDE} -S "${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/ctest_regression.cmake" \
  -DDEBUGLEVEL:INT=8 \
  -DTYPE:STRING=${OPTARG} \
  -DDISTRIB:STRING="Jenkins" \
  -DCOMPILER:STRING="MPI_GCC_6.1.x" \
  -DPROCESSORCOUNT:INT=16 \
  -DBASETESTDIR:STRING="${WORKSPACE}" \
  -DBSCRIPTARGS:STRING="${OPTSFILEARG}" \
  -DSITEDESC:STRING="${SITEDESCARG}"
