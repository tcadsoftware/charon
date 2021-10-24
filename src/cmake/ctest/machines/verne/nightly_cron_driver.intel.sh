#!/bin/bash
#
# Typical cron entries:
#
#
# 00	20	*	*	*	nightly_cron_driver.sh DBG
# 00	21	*	*	*	nightly_cron_driver.sh OPT
#
# NOTE: Overlapping, or even simultaneous, OPT and DBG runs should
# work fine

echo " "
echo "Starting nightly Charon testing on `hostname`: `date`"
echo " "

export TRIBITS_BASE_DIR=${HOME}/Projects/Charon2/TriBITS

BASETESTDIR=/home/glhenni/Nightly-Testing/Charon2
TESTDIR=TEST.INTEL.OPT

PATH=${PATH}:${HOME}/Software/bin

BTYPE="$*"

if [ "x${BTYPE}" != "xOPT" -a "x${BTYPE}" != "xDBG" ]
then
  echo "ERROR: Must specify type of nightly build, DBG or OPT"
  exit 1
fi

## These can't be set because if they are the web tools built in to
## ctest will try to go through the proxy to get to verne, which
## won't work
unset http_proxy
unset https_proxy

# Intel environment
BASEDIR=/opt/intel/2019
ARCH=intel64
unset LD_LIBRARY_PATH
# Serial compilers
source ${BASEDIR}/bin/compilervars.sh -arch ${ARCH} -platform linux

# Compilers for MPI
I_MPI_CC=icc
I_MPI_CXX=icpc
I_MPI_F77=ifort
I_MPI_F90=ifort
export I_MPI_CC I_MPI_CXX I_MPI_F77 I_MPI_F90

# Math libraries
source ${BASEDIR}/mkl/bin/mklvars.sh ${ARCH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/Software/lib

# Clean out the test directory. The ctest_empty_binary_directory()
# cmake function is too cautious to do this and fails often.
if [ -d ${BASETESTDIR}/${TESTDIR} ]
then

  # Try to delete it "safely"
  if [ "x${BASETESTDIR}" != "x" -a "x${TESTDIR}" != "x" ]
  then
    cd ${BASETESTDIR}
    if [ $? -eq 0 ]
    then
      rm -rf ${TESTDIR}
    fi
  else
      echo "ERROR:"
      echo "  BASETESTDIR=${BASESTESTDIR}"
      echo "  TESTDIR=${TESTDIR}"
      exit 1
  fi

  mkdir -p ${BASETESTDIR}/${TESTDIR}
  cd ${BASETESTDIR}
fi

# This has to be set for some Charon tests to pass
export I_MPI_SHM_LMT=shm

# Disable leak detection on nightlies on the debug build. It isn't
# done on the opt build anyway.
EXTRAOPTS=""

if [ "x${BTYPE}" = "xDBG" ]
then
  EXTRAOPTS="-LE debugexclude"
fi

PROC_COUNT="20"
export MAKEFLAGS="-j${PROC_COUNT}"

EXTRAOPTS="-L nightly ${EXTRAOPTS}"

ctest ${EXTRAOPTS} -j${PROC_COUNT} -S ${BASETESTDIR}/scripts/ctest_regression.cmake \
      -DTYPE:STRING=${BTYPE} \
      -DPROCESSORCOUNT:INT=${PROC_COUNT} \
      -DDISTRIB:STRING="Ubuntu_18.04" \
      -DCOMPILER:STRING="IntelMPI_19.x" \
      -DBASETESTDIR:STRING="${BASETESTDIR}" \
      -DTRIBRANCH:STRING="develop" \
      -DTBLDDIR:STRING="TEST.INTEL.OPT" \
      -DBSCRIPTARGS="-f linux-intel-mpi.opts"

CTEST_STAT=$?

echo " "
echo "Finished nightly Charon testing on `hostname`: `date`"
echo " "

exit $CTEST_STAT
