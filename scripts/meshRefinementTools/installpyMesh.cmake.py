#! /usr/bin/env python3

import sys
import os
import shutil

######################################################
# Set build and install paths
######################################################
buildDir = sys.argv[1]
installDir = sys.argv[2]

meshDir = installDir+"/meshRefinementTools"
binDir = installDir+"/bin"
buildSrc = buildDir+"/src/meshRefinementTools"

######################################################
# Set  scripts
######################################################
pymesh = binDir+"/pyMesh"
pymeshpy = binDir+"/pyMesh.py"
pymeshsource = meshDir+"/pyMesh.py"

######################################################
# clean up install dir
######################################################
print ("CLEANING PYMESH INSTALL IN  ",installDir)

if os.path.isfile(pymesh):
    os.remove(pymesh)

if os.path.isfile(pymeshpy):
    os.remove(pymeshpy)

if os.path.isdir(meshDir):
    shutil.rmtree(meshDir)

if not os.path.isdir(installDir):
    os.mkdir(installDir)

if not os.path.isdir(binDir):
    os.mkdir(binDir)

######################################################
# copy the build
######################################################
print ("INSTALLING PYMESH IN  ",installDir)
shutil.copytree(buildSrc,meshDir)

######################################################
# create symbolic links
######################################################

os.symlink(pymeshsource,pymesh)
os.symlink(pymeshsource,pymeshpy)

