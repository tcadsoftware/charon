#!/bin/bash -xe

####################################################
#### For Debugging set and uncomment the following
#**WORKSPACE=/gpfs1/glhenni/workspace
####################################################

# Set this appropriately if the stupid license servers
# change
export LM_LICENSE_FILE="28518@cee-infra009.sandia.gov"

# This is the directory used by the nightly jenkins job that does the
# build which will be installed by this script. This is the
# "WORKSPACE" of that job.
NWSDIR=/gpfs1/jenkins/skybridge-slave/workspace/charon_skybridge_xyceCoupled_heavy

# The install location
PROJINSTDIR=/projects/charon/install/skybridge.jenkins

if [ ! -d ${PROJINSTDIR} ]
then
  echo "ERROR: Project installation directory does not exist"
  exit 1
fi

if [ ! -d ${NWSDIR} ]
then
  echo "ERROR: Workspace directory does not exist"
  exit 1
fi

# Clean out the install directories
INSTDIRS="bin gtest include lib share lib64 charonInterpreter"
cd ${PROJINSTDIR}
for idir in ${INSTDIRS}
do
  if [ -d ${idir} ]
  then
    rm -rf ${idir}
  fi
done

# Load the necessary environment
source /usr/share/lmod/lmod/init/bash
module purge
source /projects/sems/modulefiles/utils/sems-modules-init.sh

MODULE_LIST="intel \
  cmake \
  mkl \
  openmpi-intel/3.0 \
  sems-git"

for mod in ${MODULE_LIST}
do
  module load ${mod}
done

git lfs install
export NO_PROXY=".sandia.gov,localhost"

cd ${NWSDIR}/charon-boost/boost_1_68_0
./b2 -j8 install --toolset=intel --prefix=${PROJINSTDIR}

# Install charon
cd ${NWSDIR}/TEST.OPT
make -j8 install

# Copy interpreter over
rsync --delete -a src/interpreter/charonInterpreter ${PROJINSTDIR}

# Make a link from bin to the interpreter
cd ${PROJINSTDIR}/bin
chmod ug+rwx ../charonInterpreter/charonInterpreter.py
ln -f -s ../charonInterpreter/charonInterpreter.py .

# Install xyce
cd ${NWSDIR}/Xyce/build
make -j8 install prefix=${PROJINSTDIR}

# xyce doesn't install all the necessary headers so just copy every header
# from xyce into the install directory
cd ${NWSDIR}/Xyce/src
find . -name '*.h' -exec cp {} ${PROJINSTDIR}/include \;
cd ${NWSDIR}/Xyce/build/src
find . -name '*.h' -exec cp {} ${PROJINSTDIR}/include \;

# Set permissions correctly for access by wg-charon-users group
cd ${PROJINSTDIR}

# Change group. Doing individual directories avoids the fact that this
# script can't change permissions on the root directory.
INSTDIRS="bin gtest include lib share charonInterpreter"
for idir in ${INSTDIRS}
do
  if [ -d ${idir} ]
  then
    chgrp --recursive wg-charon-users ${idir}
  fi
done

find . ! -path . -type d -print|xargs chmod 775
find . -type f -print|grep -v 'setup-run-env.'|xargs chmod g+rw
find . -type f -perm -700 -print|grep -v 'setup-run-env.'|xargs chmod 775

exit 0
