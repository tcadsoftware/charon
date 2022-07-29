#! /usr/bin/env python3

import sys
import os
import shutil
import glob

print ("arg1 = ",sys.argv[1])
print ("arg2 = ",sys.argv[2])

######################################################
# Set source  and build paths
######################################################
srcDir =  sys.argv[1]
buildDir = sys.argv[2]

meshDir = srcDir+"/scripts//meshRefinementTools"
buildSrc = buildDir+"/src/meshRefinementTools"
refineModDir = buildSrc+"/refineModule"

######################################################
# Set  scripts
######################################################
pymesh = buildDir+"/pyMesh"
pymeshpy = buildDir+"/pyMesh.py"
pymeshsource = buildSrc+"/pyMesh.py"

######################################################
# clean up links
######################################################

if os.path.isfile(pymesh):
    os.remove(pymesh)

if os.path.isfile(pymeshpy):
    os.remove(pymeshpy)


######################################################
# build refine module
######################################################

os.chdir(refineModDir)

os.system("python3 setupMeshRefineP3.py build_ext -i --inplace")

file = glob.glob('meshRefine*.so')
shutil.move(file[0],'../meshRefine.so')
os.chdir(buildDir)

######################################################
# create symbolic links
######################################################

os.symlink(pymeshsource,pymesh)
os.symlink(pymeshsource,pymeshpy)

