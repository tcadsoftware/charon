#! /usr/bin/env python3

import sys
import os
import shutil
import stat

######################################################
# Set build and install paths
######################################################
buildDir = sys.argv[1]
installDir = sys.argv[2]

dakotaDir = installDir+"/Dakota"
binDir = installDir+"/bin"
buildSrc = buildDir+"/src/Dakota"

######################################################
# Set  scripts
######################################################
dakotadriverpy = binDir+"/charonDriver.py"
dakotasource = dakotaDir+"/driver/driver.py"

vtdriverpy = binDir+"/thresholdVoltage"
vtsource = dakotaDir+"/driver/thresholdVoltage.py"

cavdriverpy = binDir+"/currentAtVoltage"
cavsource = dakotaDir+"/driver/currentAtVoltage.py"

######################################################
# clean up install dir
######################################################
print ("CLEANING DAKOTA DRIVER INSTALL IN  ",installDir)

if os.path.isfile(dakotadriverpy):
    os.remove(dakotadriverpy)

if os.path.isfile(vtdriverpy):
    os.remove(vtdriverpy)

if os.path.isfile(cavdriverpy):
    os.remove(cavdriverpy)

if os.path.isdir(dakotaDir):
    shutil.rmtree(dakotaDir)

if not os.path.isdir(installDir):
    os.mkdir(installDir)

if not os.path.isdir(binDir):
    os.mkdir(binDir)

######################################################
# copy the build
######################################################
print ("INSTALLING DAKOTA DRIVER IN  ",installDir)
shutil.copytree(buildSrc,dakotaDir)
changeMode = "chmod -R g+rX "+dakotaDir
os.system(changeMode)

######################################################
# create symbolic links
######################################################

os.symlink(dakotasource,dakotadriverpy)
os.symlink(vtsource,vtdriverpy)
os.symlink(cavsource,cavdriverpy)
