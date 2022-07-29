#!/bin/bash -xe

####################################################
#### For Debugging set and uncomment the following
#**WORKSPACE=/gpfs1/glhenni/workspace
####################################################

# Set this appropriately if the stupid license servers
# change
export LM_LICENSE_FILE="28518@cee-infra009.sandia.gov"

# Set to "1" to remove the directories containing the test
# repositories so that they can be cloned from scratch
CLEANTESTREPO=0

# Use "n" for nightly, "h" for heavy, "x" for xyce heavy
TEST_TYPE="n"

# /tmp on skybridge seems to be filling up and that causes
# the compiler to crash. Use this instead
export TMP=${WORKSPACE}/tmp
mkdir -p ${TMP}

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

# Make sure the group, wg-charon-developers, has
# full access t everything created (-rw-rw----)
umask 007

#####
#####
# Clone the LFS repositories that can't be done above because of a
# lack of LFS support in the jenkins git. When git is upgraded on
# the HPC systems this should be taken out and the jenkins facility
# used to clone the repos.
cd ${WORKSPACE}/tcad-charon/test

# Sometimes you need to clean out the directories if a previous
# clone failed, for example.
if [ "${CLEANTESTREPO}" != "0" ]
then
  rm -rf nightlyTests
  rm -rf nightlyTestsOUO
  rm -rf heavyTests
  rm -rf heavyTestsOUO
fi

if [ -d nightlyTests ]
then
  cd nightlyTests
  git pull
  cd ${WORKSPACE}/tcad-charon/test
else
  git clone git@cee-gitlab.sandia.gov:Charon/nightlyTests.git
fi

if [ -d nightlyTestsOUO ]
then
  cd nightlyTestsOUO
  git pull
  cd ${WORKSPACE}/tcad-charon/test
else
  git clone git@cee-gitlab.sandia.gov:Charon/nightlyTestsOUO.git
fi

if [ -d heavyTests ]
then
  cd heavyTests
  git pull
  cd ${WORKSPACE}/tcad-charon/test
else
  git clone git@cee-gitlab.sandia.gov:Charon/heavyTests.git
fi

if [ -d heavyTestsOUO ]
then
  cd heavyTestsOUO
  git pull
  cd ${WORKSPACE}/tcad-charon/test
else
  git clone git@cee-gitlab.sandia.gov:Charon/heavyTestsOUO.git
fi
#####
#####

cd ${WORKSPACE}

if [ "${TEST_TYPE}" = "n" ]
then
  SBATCH_ARGS=""
  REPORT_NAME="Nightly"
  PART="short,batch"
  TIME="3:59:00"
elif [ "${TEST_TYPE}" = "h" ]
then
  SBATCH_ARGS="heavy"
  REPORT_NAME="WeeklyHeavy"
  PART="short,batch"
  TIME="3:59:00"
elif [ "${TEST_TYPE}" = "x" ]
then
  SBATCH_ARGS="xyce"
  REPORT_NAME="XyceHeavy"
  PART="short,batch"
  TIME="3:59:00"
else
    echo "ERROR: Uknown script option!"
    exit 1
fi

export todaysDate=$(date +%y-%m-%d)
export pathToReport='/projects/charon/Testing/Regression/Skybridge'

# Delete files older than 30 days
find ${pathToReport} -mtime +30 -exec rm {} \;

# Create todays testing report filename
export reportName=${pathToReport}/${REPORT_NAME}Report-${todaysDate}

# Replace the length of the job in the queue as appropriate for
# nightly or weekly heavy
cat ${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/skybridge/skybridge.sbatch.sh.TEMPL | sed -e s/TIMEVAL/${TIME}/g > ${WORKSPACE}/batch-job.sh

# Submit the job to the queue
chmod +x ${WORKSPACE}/batch-job.sh
srun --output=${reportName} \
       --error=${reportName} \
       --account=FY140041 \
       --time=${TIME} \
       --partition=${PART} \
       --job-name=charon_nt \
       --nodes=1 \
       ${WORKSPACE}/batch-job.sh ${SBATCH_ARGS}

# Change group to charon developers and make sure people in that group
# have full access
find ${pathToReport} -type f -print|xargs chgrp wg-charon-developers
find ${pathToReport} -type f -print|xargs chmod g+rw
