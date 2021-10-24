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


BTYPE="$*"

if [ "x${BTYPE}" != "xOPT" -a "x${BTYPE}" != "xDBG" ]
then
  echo "ERROR: Must specify type of nightly build, DBG or OPT"
  exit 1
fi

. /usr/share/Modules/init/bash
module purge
MODLIST="intel/14.0 python/2.7 openmpi-intel gnu/4.8.1 mkl/14.0"
for mod in ${MODLIST}
do
  module load ${mod}
done

ctest -V -S ${WORKSPACE}/tcad-charon/src/cmake/ctest/machines/jenkins/ctest_nightly_jenkins.cmake -DTYPE:STRING=${BTYPE} -DPROCESSORCOUNT:INT=12 -DROOTDIR=${WORKSPACE} -DBOOSTINCDIR=${SEMS_BOOST_INCLUDE_PATH} -DBOOSTLIBDIR=${SEMS_BOOST_LIBRARY_PATH} -DNETCDFINCDIR=${SEMS_NETCDF_INCLUDE_PATH} -DNETCDFLIBDIR=${SEMS_NETCDF_LIBRARY_PATH} -DHDF5INCDIR=${SEMS_HDF5_INCLUDE_PATH} -DHDF5LIBDIR=${SEMS_HDF5_LIBRARY_PATH} -DLAPACKLIBRARYNAMES=${LAPACK_LIBRARY_NAMES} -DLAPACKLIBRARYDIR=${LAPACK_LIBRARY_DIRS} -DLAPACKLIBRARYNAMES=${BLAS_LIBRARY_NAMES} -DLAPACKLIBRARYDIR=${BLAS_LIBRARY_DIRS}


echo " "
echo "Finished nightly Charon testing on `hostname`: `date`"
echo " "
