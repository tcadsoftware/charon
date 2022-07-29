#!/bin/bash -xe

## The following are only valid when script is submitted via
## sbatch. For salloc these must be specified on the command line.
#SBATCH --time=TIMEVAL              # Wall clock time (HH:MM:SS)
#SBATCH --account=FY200136          # WC ID
#SBATCH --job-name=charon_nt        # Name of job

nodes=$SLURM_JOB_NUM_NODES          # Number of nodes

# Number MPI processes to run on each node (a.k.a. PPN)
#
# manzaon: 2 slots per node x 18 cores per slot = 36 cores per node
cores=$(( ${nodes}*36 ))

####################################################
########### BEGIN DEBUG STUFF ######################
# Uncomment and set appropriately for debugging when 
# not running via Jenkins 
##export WORKSPACE=/nscratch/glhenni/workspace
#
# Always run heavy tests to avoid forgetting it
# on the command line
##set -- "-c" "-x"
########### END DEBUG STUFF ########################
####################################################

# Installation directory used during the build process
BLDINSTDIR=${WORKSPACE}/installs

# Directory for subsequent jenkins build job to install. This isn't
# used by anything here, it simply embeds the location into the
# makefiles in the final phase of the build. There is no "make
# install" done after that phase in this script though. That is done
# by a separate Jenkins job that runs after the one that submits this
# job. (attaway and manzano are functionally equivalent so you will
# see "attaway" in a few places within the script).
PROJINSTROOT=/projects/charon/install/attaway.jenkins

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

# Turn on dakota and pymesh for manzano builds
BLDSCRIPTARG="${BLDSCRIPTARG} --reset-option=tcad-charon_ENABLE_DAKOTA_DRIVERS:ON"
BLDSCRIPTARG="${BLDSCRIPTARG} --reset-option=tcad-charon_ENABLE_PYMESH:ON"

# To find cython
export PYTHONPATH="/projects/charon/pythonTools:${PYTHONPATH}"

CHARONROOT=${WORKSPACE}/tcad-charon

BLDTYPE="OPT"

export TRIBITS_BASE_DIR="${WORKSPACE}/TriBITS"

# /tmp on manzano seems to be filling up and that causes
# the compiler to crash. Use this instead.
export TMP=${WORKSPACE}/tmp
mkdir -p ${TMP}

# Initialize
trilinosUpdate=false
developToMaster=false
heavyTest=false
xyceTest=false
nightlyTest=false
cdeGnuBuild=false
cdeIntelBuild=false

# Default if no arguments
if [ $# = 0 ]
then
  nightlyTest=true
  cdeGnuBuild=true
fi

while getopts "tdxhnig" opt; do
    case "$opt" in
        t)
            trilinosUpdate=true
            developToMaster=false
            ;;
        d)
            trilinosUpdate=false
            developToMaster=true
            ;;
        x)
            xyceTest=true
            heavyTest=false
            ;;
        h)
            xyceTest=false
            heavyTest=true
            ;;
        n)
            nightlyTest=true
            ;;
        i)
            cdeIntelBuild=true
            cdeGnuBuild=false
            ;;
        g)
            cdeIntelBuild=false
            cdeGnuBuild=true
            ;;
        \?)
            exit 1
            ;;
    esac
done

shift $(( $OPTIND - 1 ))

if ${heavyTest}
then
  SITEDESC="manzano heavy"
  TRACK="Heavy"
  LARG="^heavy$"
elif ${xyceTest}
then
  BLDSCRIPTARG="${BLDSCRIPTARG} -f with-xyce-tri-libs.opts"
  SITEDESC="manzano Xyce Coupled heavy"
  TRACK="XyceCoupled"
  LARG="^xycecoupledheavy$"
else
  SITEDESC="manzano"
  LARG="nightly"
  TRACK=""
fi

if ${trilinosUpdate}
then
  SITEDESC="${SITEDESC} TrilinosUpdate"
  TRACK="TrilinosUpdate"
fi

if ${developToMaster}
then
  SITEDESC="${SITEDESC} DevelopToMaster"
  TRACK="DevelopToMaster"
fi

# Set up the build environment
. /etc/profile
module purge
module load cde/v2/cmake/3.19.2

if ${cdeIntelBuild}
then
  module load cde/v2/cmake/3.19.2
  module load cde/v2/compiler/intel/19.1.2
  module load cde/v2/compiler/gcc/7.2.0
  module load cde/v2/intel/19.1.2/openmpi/4.0.5
  module load cde/v2/intel/19.1.2/intel-mkl/2020.3.279
  module load cde/v2/intel/19.1.2/boost/1.73.0
  module load cde/v2/intel/19.1.2/hdf5/1.10.6
  module load cde/v2/intel/19.1.2/netcdf-c/4.7.3
  module load cde/v2/intel/19.1.2/parallel-netcdf/1.12.1

  COMPILER="CDE_V2_OpenMPI_4.0.5_Intel_19.1.2"

  # CDE options
  BLDSCRIPTARG="${BLDSCRIPTARG} -f attaway-cde-intel.opts"
  
  ACTLINKNAME="cde.intel.active"
  find_install_suffix $ACTLINKNAME
  PROJINSTDIR="${PROJINSTROOT}/cde.intel.${INSTSUFFIX}"
else
  module load cde/v2/cmake/3.19.2
  module load cde/v2/compiler/gcc/7.2.0
  module load cde/v2/gcc/7.2.0/openmpi/4.0.5
  module load cde/v2/intel/19.1.2/intel-mkl/2020.3.279
  module load cde/v2/gcc/7.2.0/boost/1.73.0
  module load cde/v2/gcc/7.2.0/hdf5/1.10.6
  module load cde/v2/gcc/7.2.0/netcdf-c/4.7.3
  module load cde/v2/gcc/7.2.0/parallel-netcdf/1.12.1
    
  COMPILER="CDE_V2_OpenMPI_4.0.5_GNU_7.2.0"

  # CDE options
  BLDSCRIPTARG="${BLDSCRIPTARG} -f attaway-cde-gnu.opts"

  ACTLINKNAME="cde.gnu.active"
  find_install_suffix $ACTLINKNAME
  PROJINSTDIR="${PROJINSTROOT}/cde.gnu.${INSTSUFFIX}"
fi

# Clean out the test and install directories. The
# ctest_empty_binary_directory() cmake function is too cautious to do
# this and fails often.
if [ -d ${WORKSPACE}/TEST.OPT ]
then
  rm -r ${WORKSPACE}/TEST.OPT
  mkdir ${WORKSPACE}/TEST.OPT
fi
if [ -d ${BLDINSTDIR} ]
then
  rm -r ${BLDINSTDIR}
fi
cd ${WORKSPACE}

# For coupled xyce-charon builds
if [ "${TRACK}" = "XyceCoupled" ]
then

  # Build charon with the libraries also required by xyce
  ctest \
    -j${cores} \
    -S ${CHARONROOT}/src/cmake/ctest/machines/ctest_regression.cmake \
    -DDEBUGLEVEL:INT=8 \
    -DTYPE:STRING=${BLDTYPE} \
    -DTESTTRACK:STRING="${TRACK}" \
    -DDISTRIB:STRING="manzano" \
    -DCOMPILER:STRING="${COMPILER}" \
    -DPROCESSORCOUNT:INT=${cores} \
    -DBASETESTDIR:STRING="${WORKSPACE}" \
    -DBSCRIPTARGS:STRING="${BLDSCRIPTARG}" \
    -DSITEDESC:STRING="${SITEDESC}" \
    -DBUILDONLY:BOOL=TRUE \
    -DINSTALLDIR:STRING="${BLDINSTDIR}/charon-for-xyce" \
    -DDONOTUPDATE:BOOL=TRUE \
    -DDONOTCLONE:BOOL=TRUE

  # Install the trilinos libraries from the previous build
  cd ${WORKSPACE}/TEST.OPT
  make -j${cores} install

  # Build and install xyce using the trilinos libraries from the previous build
  cd ${WORKSPACE}/Xyce
  rm -rf ${WORKSPACE}/Xyce/build
  mkdir ${WORKSPACE}/Xyce/build
  cd ${WORKSPACE}/Xyce/build
  ${WORKSPACE}/tcad-charon/scripts/build/charonops/xyce-for-charon-cmake.sh
  make -j${cores} install

  # xyce doesn't install all the necessary headers so just copy every header
  # from xyce into the install directory
  cd ${WORKSPACE}/Xyce/src
  find . -follow -name '*.h' -exec cp {} ${BLDINSTDIR}/xyce-for-charon/include \;
  cd ${WORKSPACE}/Xyce/build/src
  find . -follow -name '*.h' -exec cp {} ${BLDINSTDIR}/xyce-for-charon/include \;
  cd ${WORKSPACE}

  # set up for a coupled build
  BLDSCRIPTARG="${BLDSCRIPTARG} -f xyce-cluster.cts1.opts"

  # Reset the install directory to a globally accessible location for
  # general use to be used by a separate jenkins job later. This
  # essentially just embeds the location for that install into the
  # makefiles. It's not actually used here since there is no "make
  # install" after the final phase of the build.
  BLDINSTDIR=${PROJINSTDIR}

  # Clean up build directory for next phase of testing
  cd ${WORKSPACE}
  rm -rf TEST.OPT
  mkdir -p TEST.OPT

fi

EXTRAOPTS=""
if [ "x${BLDTYPE}" = "xDBG" ]
then
  EXTRAOPTS="-LE debugexclude"
fi

# Avoid proxy problems with cdash
export no_proxy=".sandia.gov"

# Make sure the group, wg-charon-developers, has
# full access t everything created.
umask 007

# For debugging
###EXTRAOPTS="${EXTRAOPTS} -VV"

ctest \
  -L ${LARG} ${EXTRAOPTS} -j${cores} \
  -S ${CHARONROOT}/src/cmake/ctest/machines/ctest_regression.cmake \
  -DDEBUGLEVEL:INT=8 \
  -DTYPE:STRING=${BLDTYPE} \
  -DTESTTRACK:STRING="${TRACK}" \
  -DDISTRIB:STRING="manzano" \
  -DCOMPILER:STRING="${COMPILER}" \
  -DPROCESSORCOUNT:INT=${cores} \
  -DBASETESTDIR:STRING="${WORKSPACE}" \
  -DBSCRIPTARGS:STRING="${BLDSCRIPTARG}" \
  -DSITEDESC:STRING="${SITEDESC}" \
  -DINSTALLDIR:STRING="${PROJINSTDIR}" \
  -DDONOTUPDATE:BOOL=TRUE \
  -DDONOTCLONE:BOOL=TRUE

