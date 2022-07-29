#! /usr/bin/env python3

import os
from os import path
from os import listdir
from os.path import exists, join
import subprocess
import argparse

from interpreterParser import *
from domainDecomposition import *
from specialInformationHandler import *
from processInputFile import *

#######################################################################################################
##  create the interpreter and argument parser objects
#######################################################################################################

iP = interpreterParser()
parser = argparse.ArgumentParser()

#######################################################################################################
##  create the domain decomposition objects
#######################################################################################################

domainDecomp = domainDecomposition()

#######################################################################################################
##  create the special information handler object
#######################################################################################################

sIHandler = specialInformationHandler()

#######################################################################################################
##  Get the path to the charon executable and the executable name if they exist in the environment
#######################################################################################################

def checkForCharon():
    searchPath = os.environ['PATH'].split(":")
    for sP in searchPath:
        if exists(join(sP,charon_executable)):
            return sP

    return ""

#######################################################################################################
##  Initialize some variables
#######################################################################################################


generateHelp = "None"
readFileName = False

numProcs = 1

readFilename = False
filename = ""
plFilename = "NoName"
logFilename = ""
logFilenamePlus = ""
silentOutput = False
logRun = True


verboseOutput = False
verbosity = 0

executeRun = False
path_to_charon = ""
writeCurrents = " --current"

decomposeDomain = False
resizeDomain = False
oldNumProcs = 0

cleanTextData = False

specialInformation = []

#######################################################################################################
##  add command line options 
#######################################################################################################

#parser.add_argument("filename",help="Specify the intpreter input file name.")
parser.add_argument("-i","--input",help="Specify the interpreter input file")
parser.add_argument("-r","--run",help="Execute Charon",action="store_true")
parser.add_argument("--silent",help="Supresses screen output.  Screen output will be stored in log.",action="store_true")
parser.add_argument("--no_current",help="Disable text output of current",action="store_true")
parser.add_argument("--np",help="Specify the number of processors to use in the Charon run.",type=int)
parser.add_argument("-s","--syntax",help="Get help for syntax of interpreter. Append \"| less -i\" to scrolling and search.",action="store_true")
parser.add_argument("-S","--longsyntax",help="Get expanded help for syntax of interpreter. Append \"| less -i\" to scrolling and search.",action="store_true")
parser.add_argument("-d","--decomp",help="Decompose the state file for parallel execution.  Will abort if decomposition of same size already exists.",action="store_true")
parser.add_argument("-R","--resize_from_nprocs",help="Resizes the decomposition from nprocs domains to np domains.  Will abort if decomposition of same size is detected.",type=int)

parser.add_argument("-v","--verbosity",help="Specify verbose output of the interpreter ranging from 0 to 25",type=int)
parser.add_argument("-p","--path_to_charon",help="Specify path to Charon exectutable")
parser.add_argument("--cleanTextData",help="This option will remove sweep and transient text data files prior to execution.",action="store_true")
parser.add_argument("--mpitile",help="uses Dakota's mpitile on large computing systems with high concurrency.",action="store_true")
parser.add_argument("-l","--label",help="Specify an optional label for the run.")
parser.add_argument("--nolog",help="Opts out of logging the run.",action="store_true")
parser.add_argument("--mpipath",help="Specify path to an mpi run command.  This path may also be set through the MPIRUN_COMMAND_PATH environment variable.")
parser.add_argument("--mpiexec",help="Specify an mpi executable name.  This name may also be set through the MPIRUN_COMMAND environment variable.")

args = parser.parse_args()

#######################################################################################################
##  process command line options 
#######################################################################################################

if args.run:
    executeRun = True
else:
    executeRun = False

if args.silent:
    silentOutput = True
else:
    silentOutput = False

mpitile = False
if args.mpitile:
    mpitile = True
else:
    mpitile = False

if args.decomp:
    decomposeDomain = True
else:
    decomposeDomain = False

if args.resize_from_nprocs != None:
    resizeDomain = True
    oldNumProcs = args.resize_from_nprocs
else:
    resizeDomain = False

if args.no_current:
    writeCurrents = ""

if args.nolog:
    runLog = False
else:
    runLog = True

if args.input != None:
    readFileName = True
    filename = args.input
    logFilename = filename+".log"
    if runLog:
        if silentOutput == True and runLog == True:
            logFilenamePlus = " |& tee -a "+logFilename+" >& /dev/null"
        else:
            logFilenamePlus = " |& tee -a "+logFilename


if args.np != None:
    numProcs = args.np
else:
    numProcs = 1

if args.syntax == True:
    generateHelp = "Short"

if args.longsyntax == True:
    generateHelp = "Long"

if args.verbosity != None:
    verboseOutput = True
    verbosity = args.verbosity

if args.path_to_charon != None:
    path_to_charon = args.path_to_charon+"/"

if args.mpipath != None:
    path_to_mpi = args.mpipath+"/"
else:
    path_to_mpi = ""

if args.mpiexec != None:
    mpi_exec = args.mpiexec
else:
    mpi_exec = "mpirun"



#otion to remove old loca and transient data files
if args.cleanTextData == True:
    cleanTextData = True

#Optional label for the run
runLabel = ""
runLabelSpecified = False
if args.label != None:
    runLabelSpecified = True
    runLabel = "  "+args.label

#IF a request for help has been made, do nothing else.
if generateHelp != "None":
    executeRun = False

#######################################################################################################
##  Generate syntax help
#######################################################################################################

if (generateHelp == "Short") or (generateHelp == "Long"):
    iP.generateHelp(generateHelp)

#######################################################################################################
##  open the input and pull into a list
#######################################################################################################

if readFileName == True:
    inputLines = list(open(filename))
    pI = processInputFile()
    inputLines = pI.processIncludes(inputLines)
    inputLines = pI.apreproInputFile(inputLines)
 
    plFilename = filename+".xml"
    specialInformation = iP.parseInterpreterFile(inputLines, plFilename, verbosity)
    # Now that the consolidated input file is in hand and logging is on, write it to the log file
    if os.path.isfile(logFilename):
        os.remove(logFilename)
    if runLog == True:
        logFile = open(logFilename,'w+')
        for line in inputLines:
            logFile.write(line)
        logFile.close()

#######################################################################################################
##  special instruction handling
#######################################################################################################

if cleanTextData:
    sIHandler.cleanTextData(specialInformation)

#######################################################################################################
##  Get the path to the charon executable and the executable name if they exist in the environment
#######################################################################################################

if 'CHARON_EXECUTABLE' in os.environ:
    charon_executable = os.environ['CHARON_EXECUTABLE']
else:
    charon_executable = "charon_mp.exe"

if path_to_charon == "":
    if 'CHARON_EXECUTABLE_PATH' in os.environ:
        path_to_charon = os.environ['CHARON_EXECUTABLE_PATH']+"/"
        if not exists(join(path_to_charon,charon_executable)) and executeRun == True:
            print("Interpreter cannot find ",join(path_to_charon,charon_executable),". in your CHARON_EXECUTABLE_PATH. Check the path again.")
    else:
        charon_search = checkForCharon()
        if charon_search == "" and executeRun == True:
            print ("Interpreter cannot find ",charon_executable,". set the CHARON_EXECUTABLE_PATH environment variable.")
            path_to_charon = ""
        else:
            path_to_charon = charon_search+"/"

#######################################################################################################
##  Get the path to the mpirun executable and the executable name if they exist in the environment
#######################################################################################################

if 'MPIRUN_COMMAND' in os.environ:
    if args.mpiexec == None:  ##Always defer to the command line option
        mpi_exec = os.environ['MPIRUN_COMMAND']

if path_to_mpi == "":
    if 'MPIRUN_COMMAND_PATH' in os.environ:
        if args.mpipath == None:
            path_to_mpi = os.environ['MPIRUN_COMMAND_PATH']+"/" ##Always defer to the command line option

#######################################################################################################
##  Decompose the domain
#######################################################################################################

doIContinue = True

if decomposeDomain:
    #Check if decomp was found
    if domainDecomp.decompExists():
        print("Decomposing the domain into ",str(numProcs))
        doIContinue = domainDecomp.decompose(sIHandler.getValue(specialInformation,"importStateFileName")[0],numProcs)


#######################################################################################################
##  Resize the domain
#######################################################################################################

if resizeDomain:
    #Check if decomp was found
    if domainDecomp.decompExists() and domainDecomp.epuExists():
        print("Resizing the domain from ",str(oldNumProcs)," to ",str(numProcs)," processors")
        decompFile = sIHandler.getValue(specialInformation,"importStateFileName")[0]
        doIContinue = domainDecomp.resize(decompFile,numProcs,oldNumProcs)


#######################################################################################################
##  Exit if resize or decompose fails
#######################################################################################################

if not doIContinue:
    sys.exit(1)


#######################################################################################################
##  Execute Charon
#######################################################################################################


if executeRun == True:
    try:
        if mpitile == False:
            runCommand = path_to_mpi+mpi_exec+" "
        if mpitile == True:
             runCommand = "mpitile   --bind-to none "
        arguments = " -np "+str(numProcs)+" "+path_to_charon+charon_executable+" --i="+plFilename+writeCurrents
        totalCommand = runCommand+arguments+logFilenamePlus
        notification = "executing:  "
        if runLabelSpecified == True:
            notification = "executing "+runLabel+":  "
        print (notification,totalCommand)
        returnValue = subprocess.check_call(totalCommand,shell=True)

    except subprocess.CalledProcessError as e:
        print ("Interpreter aborting.",e)

