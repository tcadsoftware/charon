#! /usr/bin/env python3

import sys
import os
import shutil
import stat
import threading
import time

######################################################
# Set build and install paths
######################################################
buildDir = sys.argv[1]
installDir = sys.argv[2]

chirpDir = installDir+"/charonInterpreter"
binDir = installDir+"/bin"
buildSrc = buildDir+"/src/interpreter/charonInterpreter"

######################################################
# Set  scripts
######################################################
chirp = binDir+"/chirp"
charoninterpreterpy = binDir+"/charonInterpreter.py"
charoninterpreter = binDir+"/charonInterpreter"
chirpsource = chirpDir+"/charonInterpreter.py"

######################################################
# clean up install dir
######################################################
print ("CLEANING INTERPRETER INSTALL IN  ",installDir)

if os.path.isfile(chirp):
    os.remove(chirp)

if os.path.isfile(charoninterpreterpy):
    os.remove(charoninterpreterpy)

if os.path.isfile(charoninterpreter):
    os.remove(charoninterpreter)

if os.path.isdir(chirpDir):
    shutil.rmtree(chirpDir)

if not os.path.isdir(installDir):
    os.mkdir(installDir)

if not os.path.isdir(binDir):
    os.mkdir(binDir)

######################################################
# copy the build
######################################################
print ("INSTALLING INTERPRETER IN  ",installDir)
shutil.copytree(buildSrc,chirpDir)

######################################################
# create symbolic links
######################################################

if os.path.isfile(chirpsource):
    os.symlink(chirpsource,chirp)
    os.symlink(chirpsource,charoninterpreter)
    os.symlink(chirpsource,charoninterpreterpy)

