#!/bin/bash -xe
#######################################################################
# This script should be pasted into the "Execute shell" portion of the
# "Build" page on a Jenkins charon test item. Right now this script
# supports a charon-xyce coupled build on rhel7 and skybridge.
#
# Note the first two variables need to be set appropriately.
#######################################################################

# Change this to rhel6.opts or rhel7.opts as appropriate
PYSCRIPTARG="-f rhel7.opts -f rhel7-with-xyce-tri-libs.opts"

# Change this to RHEL6 or RHEL7 as appropriate
SITEDESCARG="JBF RHEL7 Xyce Coupled"

export TRIBITS_BASE_DIR=${WORKSPACE}/TriBITS

source /etc/profile.d/modules.sh
source /projects/sems/modulefiles/utils/sems-modules-init.sh

module load sems-cmake
module load sems-gcc/6.1.0
module load sems-openmpi/1.10.1
module load sems-python
module load sems-hdf5/1.8.12/base
module load sems-netcdf/4.4.1/exo

# Build boost
cd ${WORKSPACE}/charon-boost
./buildBoost.sh
cd ${WORKSPACE}

# Add library path for boost that was installed above
[[ ":${LD_LIBRARY_PATH}:" != *":${WORKSPACE}/install/boost-1_68_0/lib:"* ]] && LD_LIBRARY_PATH="${WORKSPACE}/install/boost-1_68_0/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH

# Save the library path in order to reset it later
OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

# Do the build, including Trilinos libraries that will be needed in
# the next phase for the Xyce build
ctest -j 16 -L nightly ${DBGEXCLUDE} -S "${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/ctest_regression.cmake" \
  -DDEBUGLEVEL:INT=8 \
  -DTYPE:STRING="OPT" \
  -DDISTRIB:STRING="Jenkins" \
  -DCOMPILER:STRING="MPI_GCC_6.1.x" \
  -DPROCESSORCOUNT:INT=16 \
  -DBASETESTDIR:STRING="${WORKSPACE}" \
  -DBSCRIPTARGS:STRING="${PYSCRIPTARG}" \
  -DSITEDESC:STRING="${SITEDESCARG}" \
  -DBUILDONLY:BOOL=TRUE

cd ${WORKSPACE}/TEST.OPT
make -j16 install
cd ${WORKSPACE}

# Add the path to the charon/trilinos libraries built above to the
# path
[[ ":${LD_LIBRARY_PATH}:" != *":${WORKSPACE}/install/lib:"* ]] && LD_LIBRARY_PATH="${WORKSPACE}/install/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH

# Now build Xyce, linking to the libraries created in the previous build
cd Xyce
${WORKSPACE}/tcad-charon/scripts/build/jenkins/xyce-conf-c2.sh
cd ${WORKSPACE}/Xyce/build
make -j16 install prefix=${WORKSPACE}/install

# xyce doesn't install all the necessary headers so just copy every header
# from xyce into the install directory
cd ${WORKSPACE}/Xyce/src
find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;
cd ${WORKSPACE}/Xyce/build/src
find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;
cd ${WORKSPACE}

# Clean out the charon build directory
rm -rf ${WORKSPACE}/TEST.OPT

# Reset the library path to avoid using non-cluster libraries for the
# following cluster build
export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}

# Finally, build charon again, linking to Xyce
PYSCRIPTARG="-f rhel7.opts -f rhel7-with-xyce-tri-libs.opts -f xyce-cluster.opts"

ctest -j 16 -L nightly ${DBGEXCLUDE} -S "${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/ctest_regression.cmake" \
  -DDEBUGLEVEL:INT=8 \
  -DTYPE:STRING="OPT" \
  -DDISTRIB:STRING="Jenkins" \
  -DCOMPILER:STRING="MPI_GCC_6.1.x" \
  -DPROCESSORCOUNT:INT=16 \
  -DBASETESTDIR:STRING="${WORKSPACE}" \
  -DBSCRIPTARGS:STRING="${PYSCRIPTARG}" \
  -DSITEDESC:STRING="${SITEDESCARG}" \
  -DTESTTRACK:STRING="XyceCoupled" \
  -DBUILDONLY:BOOL=FALSE

exit 0
