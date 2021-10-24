#!/bin/bash -xe

####################################################
#### For Debugging set and uncomment the following
#**WORKSPACE=/gpfs1/glhenni/workspace
####################################################

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

# Install boost
cd ${NWSDIR}/charon-boost/boost_1_68_0
./b2 -j8 install --toolset=intel --prefix=${PROJINSTDIR}

# Install xyce
cd ${NWSDIR}/Xyce/build
make -j8 install prefix=${PROJINSTDIR}

# xyce doesn't install all the necessary headers so just copy every header
# from xyce into the install directory
cd ${NWSDIR}/Xyce/src
find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;
cd ${NWSDIR}/Xyce/build/src
find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;

# Install charon
cd ${NWSDIR}/TEST.OPT
make -j8 install

# Set permissions correctly for access by wg-charon-users group
cd ${PROJINSTDIR}
find . ! -path . -type d -print|xargs chmod 775
find . -type f -print|grep -v 'setup-run-env.sh'|xargs chmod g+rw
find . -type f -perm -700 -print|grep -v 'setup-run-env.sh'|xargs chmod 775

exit 0
