#!/bin/bash -x

#SBATCH --nodes=1                   # Number of nodes - all cores per node are allocated to the job
#SBATCH --time=TIMEVAL              # Wall clock time (HH:MM:SS) - once the job exceeds this time, the job will be terminated (default is 5 minutes)
#SBATCH --account=FY140041          # WC ID
#SBATCH --job-name=charon_nt        # Name of job

nodes=$SLURM_JOB_NUM_NODES           # Number of nodes - the number of nodes you have requested (for a list of SLURM environment variables see "man sbatch")
cores=16                             # Number MPI processes to run on each node (a.k.a. PPN)
                                     # TLCC2 has 16 cores per node

####################################################
########### BEGIN DEBUG STUFF ######################
# Uncomment and set appropriately for debugging
#***WORKSPACE=/gpfs1/glhenni/workspace
########### END DEBUG STUFF ########################
####################################################

CHARONROOT=${WORKSPACE}/tcad-charon

OPTARG="OPT"
PROC_COUNT=${cores}
BLDSCRIPTARG="-f chama.opts"
DONOTUPDATE=TRUE

export TRIBITS_BASE_DIR="${WORKSPACE}/TriBITS"

# /tmp on chama seems to be filling up and that causes
# the compiler to crash. Use this instead.
export TMP=${WORKSPACE}/tmp
mkdir -p ${TMP}

# Find out the type of run, eg., heavy, xyce
if [ "$#" -eq "1" ]
then
  if [ "$1" = "heavy" ]
  then
    SITEDESC="chama heavy"
    TRACK="Heavy"
    LARG="^heavy$"
  elif [ "$1" = "xyce" ]
  then
    BLDSCRIPTARG="${BLDSCRIPTARG} -f chama-with-xyce-tri-libs.opts"
    SITEDESC="chama xyce coupled heavy"
    TRACK="XyceCoupledHeavy"
    LARG="^xycecoupledheavy$"
  else
    echo "ERROR: Unknown script argument(s)"
    echo "  $*"
    exit 1
  fi
elif [ "$#" -eq "0" ]
then
  SITEDESC="chama"
  LARG="nightly"
  TRACK=""
else
  echo "ERROR: Unknown number of arguments to script!"
  exit 1
fi

# Set up the build environment
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


# Clean out the test directory. The ctest_empty_binary_directory()
# cmake function is too cautious to do this and fails often.
if [ -d ${WORKSPACE}/TEST.OPT ]
then

  # Try to delete it "safely"
  if [ "x${WORKSPACE}" != "x" ]
  then
    cd ${BASETESTDIR}
    if [ $? -eq 0 ]
    then
      rm -rf TEST.OPT
    fi
  else
      echo "ERROR:"
      echo "  WORKSPACE=${WORKSPACE}"
      exit 1
  fi

  mkdir -p TEST.OPT
  cd ${WORKSPACE}
fi

# For coupled xyce-charon builds
if [ "${TRACK}" = "XyceCoupledHeavy" ]
then

  DONOTUPDATE=TRUE

  # Clean out the install directory. If you don't do this and the xyce
  # build fails, for example, the test may incorrectly complete using
  # the previous install.
  rm -rf ${WORKSPACE}/install
  mkdir ${WORKSPACE}/install

  # Build and install boost in the local work space
  cd ${WORKSPACE}/charon-boost
  ./buildBoost.sh INTEL
  cd ${WORKSPACE}

  # Add library path for boost that was installed above
  [[ ":${LD_LIBRARY_PATH}:" != *":${WORKSPACE}/install/lib:"* ]] && LD_LIBRARY_PATH="${WORKSPACE}/install/lib:${LD_LIBRARY_PATH}"
  export LD_LIBRARY_PATH

  # Build charon with the libraries also required by xyce
  ctest -j${PROC_COUNT} \
        -S ${CHARONROOT}/src/cmake/ctest/machines/ctest_regression.cmake \
        -DDEBUGLEVEL:INT=8 \
        -DTYPE:STRING=${OPTARG} \
        -DTESTTRACK:STRING="${TRACK}" \
        -DDISTRIB:STRING="chama" \
        -DCOMPILER:STRING="MPI_INTEL_19.x" \
        -DPROCESSORCOUNT:INT=${PROC_COUNT} \
        -DBASETESTDIR:STRING="${WORKSPACE}" \
        -DBSCRIPTARGS:STRING="${BLDSCRIPTARG}" \
        -DSITEDESC:STRING="${SITEDESC}" \
        -DBUILDONLY:BOOL=TRUE \
        -DDONOTUPDATE:BOOL=${DONOTUPDATE}

  # Install the trilinos libraries from the previous build
  cd ${WORKSPACE}/TEST.OPT
  make -j${PROC_COUNT} install

  # Build and install xyce using the trilinos libraries from the previous build
  cd ${WORKSPACE}/Xyce
  ${WORKSPACE}/tcad-charon/scripts/build/jenkins/xyce-conf-c2.sh
  cd build
  make -j16 install prefix=${WORKSPACE}/install

  # xyce doesn't install all the necessary headers so just copy every header 
  # from xyce into the install directory
  cd ${WORKSPACE}/Xyce/src
  find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;
  cd ${WORKSPACE}/Xyce/build/src
  find . -name '*.h' -exec cp {} ${WORKSPACE}/install/include \;
  cd ${WORKSPACE}

  # set up for a coupled build
  BLDSCRIPTARG="-f chama.opts -f chama-with-xyce-tri-libs.opts -f chama-xyce-cluster.opts"

  # Clean out the charon build directory for a fresh start during the
  # coupled build
  cd ${WORKSPACE}
  rm -rf TEST.OPT
  mkdir -p TEST.OPT

fi

EXTRAOPTS=""
if [ "x${OPTARG}" = "xDBG" ]
then
  EXTRAOPTS="-LE debugexclude"
fi

# Make sure the group, wg-charon-developers, has
# full access t everything created.
umask 007

# For debugging
###EXTRAOPTS="${EXTRAOPTS} -VV"

ctest -L ${LARG} ${EXTRAOPTS} -j${PROC_COUNT} \
      -S ${CHARONROOT}/src/cmake/ctest/machines/ctest_regression.cmake \
           -DDEBUGLEVEL:INT=1 \
           -DTYPE:STRING=${OPTARG} \
           -DTESTTRACK:STRING="${TRACK}" \
           -DDISTRIB:STRING="chama" \
           -DCOMPILER:STRING="MPI_INTEL_19.x" \
           -DPROCESSORCOUNT:INT=${PROC_COUNT} \
           -DBASETESTDIR:STRING="${WORKSPACE}" \
           -DBSCRIPTARGS:STRING="${BLDSCRIPTARG}" \
           -DSITEDESC:STRING="${SITEDESC}" \
           -DDONOTUPDATE:BOOL=${DONOTUPDATE}

exit 0
