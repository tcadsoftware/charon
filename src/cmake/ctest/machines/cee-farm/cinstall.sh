#!/bin/bash -xe
# Installation script for CEE jenkins builds

# Root directory of installation
PROJINSTROOT=/projects/charon/install/BOD

# Function to find out what the currently user-active install
# directory is given the link name. The link points to the user-active
# directory so the install should be performed in the other
# version. For example, if the user-active link is pointing to
# cee.cde.gnu.opt.2 then this function will set the global INSTSUFFIX
# to 1 and vice versa.
function find_install_suffix () {
    local linktarg="`readlink ${PROJINSTROOT}/$1`"
    local extension="${linktarg##*.}"
    case "$extension" in
        1)
            INSTSUFFIX=2
            ;;
        2)
            INSTSUFFIX=1
            ;;
        *)
            echo "ERROR: Unknown active target"
            exit 2
            ;;
    esac
}

# Set this appropriately if the stupid license servers
# change
export LM_LICENSE_FILE="28518@cee-infra009.sandia.gov"

if (( $# == 0 ))
then
  # Set to GNU build if no arguments given
  set -- "-g"
fi

if [ ! -d ${PROJINSTROOT} ]
then
  echo "ERROR: Installation root directory \"${PROJINSTROOT}\" doesn't exist"
  exit 1
fi

while getopts "gidk" opt; do
    case "$opt" in
        g)
            # GNU opt install
            ACTLINKNAME="cee.cde.gnu.opt.active"
            find_install_suffix $ACTLINKNAME
            INSTALLDIR="cee.cde.gnu.opt.${INSTSUFFIX}"
            MODULES="cde/v1/cmake/3.17.1 \
              cde/v1/anaconda3/4.6.14 \
              cde/v1/compiler/gcc/7.2.0 \
              cde/v1/gcc/7.2.0/openmpi/3.1.6 \
              cde/v1/gcc/7.2.0/hdf5/1.10.6 \
              cde/v1/gcc/7.2.0/netcdf-c/4.7.3 \
              cde/v1/gcc/7.2.0/boost/1.73.0"
            BLDSUBDIR="TEST.OPT"
            ;;
        d)
            # GNU debug install
            ACTLINKNAME="cee.cde.gnu.dbg.active"
            find_install_suffix $ACTLINKNAME
            INSTALLDIR="cee.cde.gnu.dbg.${INSTSUFFIX}"
            MODULES="cde/v1/cmake/3.17.1 \
              cde/v1/anaconda3/4.6.14 \
              cde/v1/compiler/gcc/7.2.0 \
              cde/v1/gcc/7.2.0/openmpi/3.1.6 \
              cde/v1/gcc/7.2.0/hdf5/1.10.6 \
              cde/v1/gcc/7.2.0/netcdf-c/4.7.3 \
              cde/v1/gcc/7.2.0/boost/1.73.0"
            BLDSUBDIR="TEST.DBG"
            ;;
        i)
            # Intel opt install
            ACTLINKNAME="cee.cde.intel.opt.active"
            find_install_suffix $ACTLINKNAME
            INSTALLDIR="cee.cde.intel.opt.${INSTSUFFIX}"
            MODULES="cde/v1/cmake/3.17.1 \
              cde/v1/anaconda3/4.6.14 \
              cde/v1/compiler/intel/19.0.5.281 \
              cde/v1/intel/19.0.5.281/openmpi/3.1.6 \
              cde/v1/intel/19.0.5.281/hdf5/1.10.6 \
              cde/v1/intel/19.0.5.281/netcdf-c/4.7.3 \
              cde/v1/intel/19.0.5.281/boost/1.73.0"
            BLDSUBDIR="TEST.OPT"
            ;;
        k)
            # Intel dbg install
            ACTLINKNAME="cee.cde.intel.dbg.active"
            find_install_suffix $ACTLINKNAME
            INSTALLDIR="cee.cde.intel.dbg.${INSTSUFFIX}"
            MODULES="cde/v1/cmake/3.17.1 \
              cde/v1/anaconda3/4.6.14 \
              cde/v1/compiler/intel/19.0.5.281 \
              cde/v1/intel/19.0.5.281/openmpi/3.1.6 \
              cde/v1/intel/19.0.5.281/hdf5/1.10.6 \
              cde/v1/intel/19.0.5.281/netcdf-c/4.7.3 \
              cde/v1/intel/19.0.5.281/boost/1.73.0"
            BLDSUBDIR="TEST.DBG"
            ;;
        \?)
            echo "Valid options are -g, -i, -d or -k"
            exit 1
            ;;
    esac
done

# Set up module environment
source /etc/profile
module purge
for mod in ${MODULES}
do
    module load ${mod}
done

# Where to install
PROJINSTDIR=${PROJINSTROOT}/${INSTALLDIR}
if [ ! -d ${PROJINSTDIR} ]
then
  echo "ERROR: Project installation directory, ${PROJINSTDIR}, doesn't exist"
  exit 2
fi

BLDDIR=${WORKSPACE}/${BLDSUBDIR}

# Install from
if [ ! -d ${BLDDIR} ]
then
  echo "ERROR: Build directory, ${BLDDIR}, doesn't exist"
  exit 2
fi

# Clean out directories in the target installation directory,
# individually. This will leave any files in there, like the one to
# setup the environment, intact.
cd ${PROJINSTDIR}
files="`ls`"
for file in ${files}
do
    if [ -d ${file} ]
    then
      rm -rf ${file}
    fi
done

# Go to the build directory for installation
cd ${BLDDIR}
make -j8 install

# Set permissions on installation directory
cd ${PROJINSTDIR}
chgrp --recursive wg-charon-users .
find . -type d -print|xargs chmod 755
find . -type f -print|xargs chmod g+r,g-w,o+r,o-w
find . -type f -perm -700 -print|xargs chmod 755

exit 0
