#!/bin/bash -x
#
# Typical cron entries:
#
#
# 00	20	*	*	*	nightly_cron_driver.sh DBG
# 00	21	*	*	*	nightly_cron_driver.sh OPT
#
# NOTE: Overlapping, or even simultaneous, OPT and DBG runs should
# work fine

source /etc/profile

CHARONGRP="wg-charon-users"

export TRIBITS_BASE_DIR=/projects/charon/install/TriBITS

INSTDIR=/projects/charon/install/BOD
BASETESTDIR=/scratch/glhenni/Nightly

# Set permissions on this directory so debuggers and such can read the
# source code.
chgrp -R ${CHARONGRP} ${BASETESTDIR}
find ${BASETESTDIR} -type d -print0|xargs -0 chmod 750
find ${BASETESTDIR} -type f -print0|xargs -0 g+r

export TRIBITS_BASE_DIR=/projects/charon/install/TriBITS

echo " "
echo "Starting nightly Charon testing on `hostname`: `date`"
echo " "

if [ $# -eq 0 ]
then
  echo "ERROR: Usage is"
  echo "  $0 [OPT|DBG] INST"
fi

BTYPE="$1"

if [ "x${BTYPE}" != "xOPT" -a "x${BTYPE}" != "xDBG" ]
then
  echo "ERROR: Must specify type of nightly build, DBG or OPT"
  exit 1
fi

DOINST=0
if [ $# -eq 2 ]
then
  DOINST=1
fi

module purge
module add sems-env

MODLIST="sems-devpack-intel/17.0.1 \
sems-cmake \
sems-git/2.10.1"
for mod in ${MODLIST}
do
  module load ${mod}
done
module unload sems-python/2.7.9 
module load sems-python


PROC_COUNT="60"
export MAKEFLAGS="-j${PROC_COUNT}"

if [ "x${BTYPE}" = "xDBG" ]
then
  TESTDIR="TEST.DBG"
else
  TESTDIR="TEST.OPT"
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
  cd ${BASETESTDIR}
fi

EXTRAOPTS=""
if [ "x${BTYPE}" = "xDBG" ]
then
  EXTRAOPTS="-LE debugexclude"
fi

EXTRAOPTS="-L nightly ${EXTRAOPTS}"

if [ $DOINST -eq 1 -a $? -eq 0 ]
then
  if [ "${BTYPE}" = "OPT" ]
  then
    diradd="sems-intel.opt"
  else
    diradd="sems-intel.dbg"
  fi
fi

ctest ${EXTRAOPTS} -j${PROC_COUNT} -S ${BASETESTDIR}/scripts/ctest_regression.cmake \
      -DTYPE:STRING=${BTYPE} \
      -DPROCESSORCOUNT:INT=${PROC_COUNT} \
      -DDISTRIB:STRING="RHEL_7.x" \
      -DCOMPILER:STRING="SEMS_OpenMPI_1.10_Intel_17.x" \
      -DBASETESTDIR:STRING="${BASETESTDIR}" \
      -DTRIBRANCH:STRING="develop" \
      -DBSCRIPTARGS:STRING="-f cee-intel.opts" \
      -DSITEDESC:STRING="CEE Build Farm" \
      -DINSTALLDIR:PATH="/projects/charon/install/BOD/${diradd}"

CTEST_STAT=$?

# If requested, and ctest passed, install the build into a standard
# location
if [ $DOINST -eq 1 -a $CTEST_STAT -eq 0 ]
then
  cd ${BASETESTDIR}/TEST.${BTYPE}

  make -j${PROC_COUNT} install

  # Copy/install the interpreter
  rsync --delete -a ${BASETESTDIR}/${TESTDIR}/src/interpreter/charonInterpreter ${INSTDIR}/${diradd}

  # set permissions on the installation
  find ${INSTDIR}/${diradd} -type d -print|xargs chmod 775
  find ${INSTDIR}/${diradd} -type d -print|xargs chmod g+s
  find ${INSTDIR}/${diradd} -type f -print|xargs chmod ug+rw

  # Copy the source code
  rsync --delete --exclude='.gitignore' --exclude='.git/' -a ${BASETESTDIR}/tcad-charon ${INSTDIR}/src/

  # Logical link from bin to interpreter
  cd ${INSTDIR}/${diradd}/bin
  ln -f -s ../charonInterpreter/charonInterpreter.py .
  cd ${BASETESTDIR}/TEST.${BTYPE}

  # Make sure the wg-charon-users group members have access
  chgrp -R wg-charon-users ${INSTDIR}

  # set permissions on the source
  find ${INSTDIR}/src -type d -print|xargs chmod 775
  find ${INSTDIR}/src -type d -print|xargs chmod g+s
  find ${INSTDIR}/src -type f -print0|xargs -0 chmod ug+rw

fi

echo " "
echo "Finished nightly Charon testing on `hostname`: `date`"
echo " "

exit $CTEST_STAT
