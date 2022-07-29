#!/bin/bash -xe
#
# Typical cron entries:
#
#
# 00	20	*	*	*	nightly_cron_driver.sh DBG
# 00	21	*	*	*	nightly_cron_driver.sh OPT
#
# NOTE: Overlapping, or even simultaneous, OPT and DBG runs should
# work fine
GNU_VERS=11

echo " "
echo "Starting nightly Charon testing on `hostname`: `date`"
echo " "

# Do a build with specific version of GNU tools
export OMPI_CXX=g++-${GNU_VERS}
export OMPI_CC=gcc-${GNU_VERS}
export OMPI_F77=gfortran-${GNU_VERS}
export OMPI_FC=gfortran-${GNU_VERS}

export TRIBITS_BASE_DIR=${HOME}/Projects/Charon2/TriBITS

BASETESTDIR=/home/glhenni/Nightly-Testing/Charon2

PATH=${PATH}:/snap/bin:${HOME}/Software/bin

BTYPE="$*"

if [ "x${BTYPE}" != "xOPT" -a "x${BTYPE}" != "xDBG" -a "x${BTYPE}" != "xCOV" ]
then
  echo "ERROR: Must specify type of nightly build, DBG, OPT or COV"
  exit 1
fi

## These can't be set because if they are the web tools built in to
## ctest will try to go through the proxy to get to verne, which
## won't work
unset http_proxy
unset https_proxy

COVFLAG="FALSE"
BSCRIPTARGS=""
if [ "x${BTYPE}" = "xDBG" ]
then
  TESTDIR="TEST.GNU${GNU_VERS}.DBG"
elif [ "x${BTYPE}" = "xCOV" ]
then
  TESTDIR="TEST.GNU${GNU_VERS}.COV"
  COVFLAG="TRUE"
  BSCRIPTARGS="-f linux-gcov.opts"
else
  TESTDIR="TEST.GNU${GNU_VERS}.OPT"
fi

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
else
  mkdir -p ${BASETESTDIR}/${TESTDIR}
fi
cd ${BASETESTDIR}

# Disable leak detection on nightlies on the debug build. It isn't
# done on the opt build anyway.
EXTRAOPTS=""

if [ "x${BTYPE}" = "xDBG" -o "x${BTYPE}" = "xCOV" ]
then
  export ASAN_OPTIONS="detect_leaks=false"
  EXTRAOPTS="-LE debugexclude"
fi

PROC_COUNT="20"
export MAKEFLAGS="-j${PROC_COUNT}"

EXTRAOPTS="-L nightly ${EXTRAOPTS}"

ctest -V ${EXTRAOPTS} -j${PROC_COUNT} -S ${BASETESTDIR}/scripts/ctest_regression.cmake \
      -DTYPE:STRING=${BTYPE} \
      -DPROCESSORCOUNT:INT=${PROC_COUNT} \
      -DDISTRIB:STRING="Ubuntu_20.04" \
      -DCOMPILER:STRING="OpenMPI_4.0.3_GCC_${GNU_VERS}.x" \
      -DBASETESTDIR:STRING="${BASETESTDIR}" \
      -DTBLDDIR:STRING="${TESTDIR}" \
      -DCOVERAGE:BOOL="${COVFLAG}" \
      -DBSCRIPTARGS:STRING="${BSCRIPTARGS}" \
      -DDEBUGLEVEL:INT="8" \
      -DTESTTRACK:STRING="Experimental" \
      -DTRIBRANCH:STRING="develop"

CTEST_STAT=$?

echo " "
echo "Finished nightly Charon testing on `hostname`: `date`"
echo " "

exit $CTEST_STAT
