#!/usr/bin/env python3 
####!/projects/sems/install/rhel6-x86_64/sems/compiler/python/3.4.2/bin/python3
from __future__ import print_function


#
# Make a class of all the junk we need to refine the mesh
#

class refinementParameters:
    "object containing refinement parameters"
    activeBlocks = []
    Dimension = 0
    meshType = "noMeshType"
    numProcs = 1
    maximumRecursion = 0
    maxSurfaceRecursions = 0
    guaranteedSurfaceRecursions = 0
    minCellSize = 0.0
    minCellNSize = 0.0
    minCellPSize = 0.0
    minCellILSize = 0.0
    refinementDistance = 0.0
    refinementDistanceN = 0.0
    refinementDistanceP = 0.0
    refinementDistanceIL = 0.0
    junctionSurfacesFile = "noFileName"
    refinementMethod = "noMethodSpecified"
    xPlane = 0.0
    sizeMeasure = "average"
    xNodes = []
    yNodes = []
    zNodes = []
    BBxNodes = []
    BByNodes = []
    BBzNodes = []
    refineToLine = False
    refineToSurface = False
    ILThickness = 0
    analyticalSurfaces = 0
    journalSurfaces = "false"
    dopingFunctions = []
    dopingBounds = []
    dopingMin = []
    dopingMax = []
    dopingFunctionWidth = []
    dopingDirection = []
    dopingLocation = []
    listDopingFunctions = False
    dopingXBounds = []
    dopingYBounds = []
    dopingZBounds = []
    writeJunctions = False
    refinementFactor = 3
    autoRefine = False




#
# Informational crap
#
#-------------------------------------------------------------------
#
def info(title):
    print (title)
    print ('module name:', __name__)
    if hasattr(os, 'getppid'):  # only available on Unix
        print ('parent process:', os.getppid())
    print ('process id:', os.getpid())
#
#-------------------------------------------------------------------
#
# Hardcore refine tet decision making
#

def getRefinedTets(tets):

    localRefineTets = []
    doIrefine = [False] * (len(tets)+1)

    zcoord = []

    for tetNum in tets:
        #doIrefine = False

        if localParams.meshType == "tet":
            tetNodes = cubitMesh.get_connectivity("tet",tetNum)  #3D
        elif localParams.meshType == "tri":
            tetNodes = cubitMesh.get_connectivity("tri",tetNum)  #2D

        refine.setTetNum(tetNum)

        zcoord[:] = []
        for nodeNum in tetNodes:
            nodeCoords = cubitMesh.get_nodal_coordinates(nodeNum)
            refine.fillCoordinates(nodeCoords)
            zcoord.append(nodeCoords[2])

        if localParams.refinementMethod.lower() == "centroidsided":
            doIrefine[tetNum] = refine.doIrefineCentroidSided()
        elif localParams.refinementMethod.lower() == "2d":
            doIrefine[0] = refine.doIrefine2D()  ##deprecated
        elif localParams.refinementMethod.lower() == "xplane":
            doIrefine[tetNum] = refine.doIrefineXplane(localParams.xPlane)
        elif localParams.refinementMethod.lower() == "3d":
            doIrefine[0] = False
            doIrefine[0] = refine.doIrefine3D()
        else:
            print("No valid refinement method named",localParams.refinementMethod.lower())

        if doIrefine[0]:
            localRefineTets.append(tetNum)
            #print (localRefineTets)
        refine.resetNodeCounter()

    #print ("local refined tets on process: ", os.getppid())
    #print (localRefineTets)
    #print ()

    return localRefineTets

#
#-------------------------------------------------------------------
#
#
# Test the finish metric
#

def getRefinedTetsMetric(tets):

    localMetric = 0.0
    localMetricCounter=0

    for tetNum in tets:
        tetNodes = cubitMesh.get_connectivity("tet",tetNum)
        for nodeNum in tetNodes:
            nodeCoords = cubitMesh.get_nodal_coordinates(nodeNum)
            refine.fillCoordinates(nodeCoords)

        tetMetric = refine.doIrefineCentroidSidedFinishMetric()
        if localMetric > 0:
            localMetricCounter = localMetricCounter + 1
        localMetric = localMetric + tetMetric
        refine.resetNodeCounter()

#    localMetric = localMetric/len(tets)
    localMetric = localMetric/localMetricCounter

    return localMetric

#
#-------------------------------------------------------------------
#
# parse the refinement directive
#

def parseRefinementDirective(line):

    dimension = 0
    try:
        if line[1].lower() == "dimension":
            if line[2].lower() == "3d":
                localParams.Dimension = 3
                refine.setDimension(3)
            elif line[2].lower() == "analytical3d":
                localParams.Dimension = 3
                refine.setDimension(3)
                localParams.analyticalSurfaces = 1;
            elif line[2].lower() == "2d":
                localParams.Dimension = 2
                refine.setDimension(2)
            elif line[2].lower() == "analytical2d":
                localParams.Dimension = 2
                refine.setDimension(2)
                localParams.analyticalSurfaces = 1;
            else:
                raise Exception("Bad dimension in refinement directive. "+line[2]+" is invalid.")
        elif line[1].lower() == "setblocks":
            localParams.activeBlocks = [int(i) for i in line[2:]]
        elif line[1].lower() == "refinementfactor":
            localParams.refinementFactor = float(line[2])
        elif line[1].lower() == "setmesh":
            localParams.meshType = line[2]
        elif line[1].lower() == "setnumprocs":
            localParams.numProcs = int(line[2])
        elif line[1].lower() == "writejunctions":
            localParams.writeJunctions = True
        elif line[1].lower() == "setmaximumrecursion":
            localParams.maximumRecursion = int(line[2])
        elif line[1].lower() == "setmaximumsurfacerecursion":
            localParams.maxSurfaceRecursion = int(line[2])
        elif line[1].lower() == "setguaranteedsurfacerecursion":
            localParams.guaranteedSurfaceRecursion = int(line[2])
        elif line[1].lower() == "setmincellsize":
            localParams.minCellSize = float(line[2])
        elif line[1].lower() == "setmincellnsize":
            localParams.minCellNSize = float(line[2])
        elif line[1].lower() == "setmincellpsize":
            localParams.minCellPSize = float(line[2])
        elif line[1].lower() == "setmincellilsize":
            localParams.minCellILSize = float(line[2])
        elif line[1].lower() == "setrefinementdistance":
            if line[2].lower() == "autorefine":
                localParams.autoRefine = True
            else:
                localParams.refinementDistance = float(line[2])
        elif line[1].lower() == "setrefinementdistancen":
            localParams.refinementDistanceN = float(line[2])
        elif line[1].lower() == "setrefinementdistancep":
            localParams.refinementDistanceP = float(line[2])
        elif line[1].lower() == "setrefinementdistanceil":
            localParams.refinementDistanceIL = float(line[2])
        elif line[1].lower() == "journalsurfaces":
            localParams.journalSurfaces = "true"
        elif line[1].lower() == "refinetolinesegment":
            localParams.refineToLine = True
            localParams.xNodes.append(float(line[2]))
            localParams.yNodes.append(float(line[3]))
            localParams.xNodes.append(float(line[4]))
            localParams.yNodes.append(float(line[5]))
            localParams.ILThickness = float(line[6])
        elif line[1].lower() == "refinetosurfaceelement":
            localParams.refineToSurface = True
            localParams.xNodes.append(float(line[2]))
            localParams.yNodes.append(float(line[3]))
            localParams.zNodes.append(float(line[4]))
            localParams.xNodes.append(float(line[5]))
            localParams.yNodes.append(float(line[6]))
            localParams.zNodes.append(float(line[7]))
            localParams.xNodes.append(float(line[8]))
            localParams.yNodes.append(float(line[9]))
            localParams.zNodes.append(float(line[10]))
            localParams.xNodes.append(float(line[11]))
            localParams.yNodes.append(float(line[12]))
            localParams.zNodes.append(float(line[13]))
            localParams.ILThickness = float(line[14])
        elif line[1].lower() == "refinementboundingbox":
            if len(line) == 6:
                localParams.BBxNodes.append(float(line[2]))
                localParams.BByNodes.append(float(line[3]))
                localParams.BBxNodes.append(float(line[4]))
                localParams.BByNodes.append(float(line[5]))
                localParams.BBzNodes.append(0)
                localParams.BBzNodes.append(0)
            elif len(line) == 8:
                localParams.BBxNodes.append(float(line[2]))
                localParams.BByNodes.append(float(line[3]))
                localParams.BBzNodes.append(float(line[4]))
                localParams.BBxNodes.append(float(line[5]))
                localParams.BByNodes.append(float(line[6]))
                localParams.BBzNodes.append(float(line[7]))
        elif line[1].lower() == "sizemeasure":
            localParams.sizeMeasure = line[2]
        elif line[1].lower() == "setjunctionfilename":
            localParams.junctionSurfacesFile = line[2]
        elif line[1].lower() == "setrefinementmethod":
            localParams.refinementMethod = line[2]
            if len(line) > 3 and line[2].lower() == "xplane":
                localParams.xPlane = float(line[3])
        elif line[1].lower() == "setdopingfunction":  # set doping function
            localParams.dopingFunctions.append([line[2].lower(),line[3].lower(),line[4].lower()])
        elif line[1].lower() == "setdopingxbounds":  # set coordinate bounds
            if line[3].lower() == "x":
                localParams.dopingXBounds.append([line[2].lower(),line[3].lower(),float(line[4]),float(line[5])])
            elif line[3].lower() == "y":
                localParams.dopingYBounds.append([line[2].lower(),line[3].lower(),float(line[4]),float(line[5])])
            elif line[3].lower() == "z":
                localParams.dopingZBounds.append([line[2].lower(),line[3].lower(),float(line[4]),float(line[5])])
            else:
                raise Exception("Illegal coordinate direction in ",line)
        elif line[1].lower() == "setbendlocation":  # set coordinate bounds
            if localParams.Dimension == 2:
                localParams.dopingLocation.append([line[2].lower(),float(line[3]),float(line[4])])
            elif localParams.Dimension == 3:
                localParams.dopingLocation.append([line[2].lower(),float(line[3]),float(line[4]),float(line[5])])
            else:
                raise Exception("Unable to set function locations from line, ",line)
        elif line[1].lower() == "setdopingwidth":  # set coordinate bounds
            if localParams.Dimension == 2:
                localParams.dopingFunctionWidth.append([line[2].lower(),float(line[3]),float(line[4])])
            elif localParams.Dimension == 3:
                localParams.dopingFunctionWidth.append([line[2].lower(),float(line[3]),float(line[4]),float(line[5])])
            else:
                raise Exception("Unable to set function widths from line, ",line)
        elif line[1].lower() == "setdopingbounds":  # set doping function bounds
            localParams.dopingBounds.append([line[2].lower(),float(line[3]),float(line[4])])
        elif line[1].lower() == "setfunctiondirection":  # set doping function bounds
            localParams.dopingDirection.append([line[2].lower(),line[3]])
        elif line[1].lower() == "listdopingfunctions":  # set doping function bounds
            localParams.listDopingFunctions = True
        elif line[1].lower() == "executerefinement":
            executeRefinement()
        else:
            raise Exception("Unknown directive following refinementDirective keyword.")
    except Exception as e:
        print ("refinement directive failed:")
        print(e)
        print(traceback.format_exc())
        #print(sys.exc_info()[0])...        
        print("Terminatng!")
        print ()
        sys.exit(1)

#
#-------------------------------------------------------------------
#
# Execute refinement
#

def executeRefinement():

    print ("Initializing refine ")

    refine.init(localParams.Dimension)
    if localParams.writeJunctions == True:
        refine.setWriteJunctions(True)



    try:
        if localParams.Dimension == 2 and localParams.meshType == "tri":
            refine.sizeCellNodes(3) #2D
        elif localParams.Dimension == 3 and localParams.meshType == "tet":
            refine.sizeCellNodes(4)  #3D
        else:
            raise Exception(" Illegal combination of mesh dimension and mesh type is specified.")
    except Exception as e:
        print(e)
        sys.exit(1)

    # Turn on auto refinement around junctions
    if localParams.autoRefine == True:
        refine.setAutoRefine()

# Create refinement bounding box
    refine.setRefinementBoundingBox(localParams.BBxNodes[0],localParams.BByNodes[0],localParams.BBzNodes[0],localParams.BBxNodes[1],localParams.BByNodes[1],localParams.BBzNodes[1])
    if localParams.maxSurfaceRecursion != 0:
        refine.setMaxSurfaceRecursions(localParams.maxSurfaceRecursion)
    if localParams.guaranteedSurfaceRecursion != 0:
        refine.setGuaranteedSurfaceRecursions(localParams.guaranteedSurfaceRecursion)
# Create analytical Doping functions
    for line in localParams.dopingFunctions:
        refine.setDopingFunction(line[0].encode(),line[1].encode(),line[2].encode())
# set the bounds of the doping
    for line in localParams.dopingBounds:
        refine.setDopingBounds(line[0].encode(),float(line[1]),float(line[2]))
# set the X bounds of the doping
    for line in localParams.dopingXBounds:
        refine.setXBounds(line[0].encode(),line[1].encode(),line[2],line[3])
# set the Y bounds of the doping
    for line in localParams.dopingYBounds:
        refine.setXBounds(line[0].encode(),line[1].encode(),line[2],line[3])
# set the Z bounds of the doping
    for line in localParams.dopingZBounds:
        refine.setXBounds(line[0].encode(),line[1].encode(),line[2],line[3])
# set the bounds of the doping
    for line in localParams.dopingFunctionWidth:
        if localParams.Dimension == 2:
            refine.setDopingWidth(line[0].encode(),line[1],line[2])
        else:
            refine.setDopingWidth(line[0].encode(),line[1],line[2],line[3])
# set the locations of the doping
    for line in localParams.dopingFunctionWidth:
        if localParams.Dimension == 2:
            refine.setDopingLocation(line[0].encode(),line[1],line[2])
        else:
            refine.setDopingLocation(line[0].encode(),line[1],line[2],line[3])
# refine to a specified line segment
    if localParams.refineToLine == True:
        for i in range(int((len(localParams.xNodes)/2))):
            refine.addRefineToLine(localParams.xNodes[i*2],localParams.yNodes[i*2],localParams.xNodes[i*2+1],localParams.yNodes[i*2+1],localParams.ILThickness)

# refine to a specified surface element
    if localParams.refineToSurface == True:
        numRefineSurfs = 0
        refine.addRefineToSurface(localParams.xNodes[numRefineSurfs*4],localParams.yNodes[numRefineSurfs*4],localParams.zNodes[numRefineSurfs*4],
                                  localParams.xNodes[numRefineSurfs*4+1],localParams.yNodes[numRefineSurfs*4+1],localParams.zNodes[numRefineSurfs*4+1],
                                  localParams.xNodes[numRefineSurfs*4+2],localParams.yNodes[numRefineSurfs*4+2],localParams.zNodes[numRefineSurfs*4+2],
                                  localParams.xNodes[numRefineSurfs*4+3],localParams.yNodes[numRefineSurfs*4+3],localParams.zNodes[numRefineSurfs*4+3],
                                  localParams.ILThickness)
        

# List the functions if requested
    if localParams.listDopingFunctions == True:
        refine.listFunctions()

# Read in or create surfaces for refinement
    try:
        if localParams.refinementMethod.lower() != "xplane" and localParams.junctionSurfacesFile == "noFileName" and localParams.analyticalSurfaces != 1:
            raise Exception ("You have specified a refinement method which requires a data file of junction surfaces.  Cannot continue.")
        if localParams.junctionSurfacesFile != "noFileName":
            refine.readSurfaces(localParams.junctionSurfacesFile)
#        if localParams.Dimension == 3 and localParams.analyticalSurfaces == 1:
        if localParams.analyticalSurfaces == 1:
            refine.createSurfaces()
    except Exception as e:
        print(e)
        sys.exit(1)

# create a journal file of the surfaces if requested

    if localParams.journalSurfaces == "true":
        refine.journalSurfaces()


# set the metric type for a size measure
    refine.setSizeMeasure(localParams.sizeMeasure.lower().encode())


#set cell minimums
    refine.setCellMinimum(localParams.minCellSize)
    refine.setCellMinimumN(localParams.minCellNSize)  #3D
    refine.setCellMinimumP(localParams.minCellPSize)    #3D
    refine.setCellMinimumIL(localParams.minCellILSize)  #3D
#set refinement distances
    refine.setRefinementDistance(localParams.refinementDistance)
    refine.setRefinementDistanceN(localParams.refinementDistanceN)
    refine.setRefinementDistanceP(localParams.refinementDistanceP)
    refine.setRefinementDistanceIL(localParams.refinementDistanceIL)
#Set refinement factor
    refine.setRefinementFactor(localParams.refinementFactor)

#write the xml equivalent doping parameters
    refine.printXMLDopingFunctions()

# Initialize some timing variables and a recursion counter
    times = []
    times.append(time())
    recursionCounter = 0
    numProcs = localParams.numProcs

    while recursionCounter < localParams.maximumRecursion:
        print ("executing recursion layer",recursionCounter)
        refineTets = []

        for block in localParams.activeBlocks:

            if localParams.Dimension == 2:
                vols = cubitMesh.get_block_surfaces(block)  #2D
            elif localParams.Dimension == 3:
                vols = cubitMesh.get_block_volumes(block)  #3D

            for volNum in vols:
                if localParams.meshType == "tet":
                    tets = cubitMesh.get_volume_tets(volNum)  #3Dtets
                elif localParams.meshType == "tri":
                    tets = cubitMesh.get_surface_tris(volNum)  #2D

                tetLength = len(tets)
                tetList = []
                offset = int(tetLength/numProcs)
                tetList.append(list(tets[:offset]))
                for nP in range(1,numProcs-1):
                    tetList.append(list(tets[offset*nP:offset*(nP+1)]))
                if(numProcs > 1):
                    tetList.append(list(tets[offset*(numProcs-1):]))

                if __name__ == '__main__':
                    p = Pool(numProcs)
                    tmpRefinedTets = []
                    tmpRefinedTets = p.map(getRefinedTets,tetList)
                    p.terminate()
                    refineTets.extend(tmpRefinedTets)

        tetCount = 0
        for tC in refineTets:
            #print ("Counting tets ",tC)
            tetCount = tetCount + len(tC)

        print (tetCount," ",len(refineTets)," in list ")

##Create the refine command and fill it
        times.append(time())
        if tetCount > 0:
            print ("refining ",tetCount," tetrahedra ")
            if localParams.meshType == "tet":
                refineCommand = "refine tet "  #3D tet
            elif localParams.meshType == "tri":
                refineCommand = "refine tri "
            for refTet in refineTets:
                for refTetInner in refTet:
                    refineCommand = refineCommand + str(refTetInner) + " "
            refineCommand = refineCommand + " numsplit 1 Depth 1"
            print (" refine command \n", refineCommand)

            cubitMesh.cmd(refineCommand)
        times.append(time())


        recursionCounter = recursionCounter + 1



#
#-------------------------------------------------------------------
#
#
#
##Main part of script
#
#


#
#Import all modules
#
import sys
import os
from os import path
import sysconfig
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import meshRefine
from time import time, sleep
from pyMeshPrepro import *

import multiprocessing
from multiprocessing import Pool

from processInputFile import *

#Get the cubit modules
#    sys.path.append('/projects/cubit/claro.Lin64.15.5/bin/')

if 'CUBIT_PATH' in os.environ:
    cubitPath = os.environ['CUBIT_PATH']
    sys.path.append(cubitPath)
else:
    print ("Error.  No path to Cubit specied.  Must set the CUBIT_PATH environment variable to the location of Cubit.")
    sys.exit(1)

import cubit as cubitMesh

#for fun in dir(cubitMesh):
 #   print(fun)


#
# initialize the refinement module and input file processor
#

refine = meshRefine.PyMeshRefine()

pI = processInputFile()

#
#Initialize cubit
#

cubitMesh.init(['help'])
cubitMesh.cmd('record "pyMesh.jou"')

pMPrepro = pyMeshPrepro(cubitMesh)

#
# Initialize the Dimension Variable
#

#Dimension = 0

#
# parse the journal file
#

print ("Reading from input file ",sys.argv[1])

meshJournalLines = list(open(sys.argv[1]))  # list of strings

meshJournalLines = pI.processIncludes(meshJournalLines)
meshJournalLines = pI.apreproInputFile(meshJournalLines)
logFile = sys.argv[1]+".log"
pI.writeConsolidatedInputFile(meshJournalLines,logFile)

localParams = refinementParameters

for line in meshJournalLines:
    lineTokens = line.split()
    apreproList = list(cubitMesh.get_aprepro_vars())
    if len(lineTokens) != 0 and lineTokens[0].lower() == "refinementdirective":
        lineProcessed = pMPrepro.processLine(line)
        lineTokensProcessed = lineProcessed.split()
        print (lineProcessed)
        parseRefinementDirective(lineTokensProcessed)
    else:
        cubitMesh.cmd(line)

# As of Cubit v15.4, destroy() causes a fault and dumps core.
# It's not clear that this must absolutely be done as it comes at the end of the job.
# Follow up with Cubit team or remove altogether. --LCM
cubitMesh.cmd('record stop')
cubitMesh.cmd('exit')


if localParams.writeJunctions == True:
    cubitMesh.cmd('reset')
    cubitMesh.cmd('record "pySat.jou"')
    junctionJournalLines = list(open("junctionSurfaces.jou"))  # list of strings
    for line in junctionJournalLines:
        lineTokens = line.split()
        cubitMesh.cmd(line)
    cubitMesh.cmd('record stop')
    cubitMesh.cmd('exit')

#print localParams

print ("Dimension = ", refine.getDimension())
print ()
print ("blocks  = ",localParams.activeBlocks)
print ()
print ("Mesh Type =  ",localParams.meshType)
