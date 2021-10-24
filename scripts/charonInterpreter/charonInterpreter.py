#! /usr/bin/env python3
from __future__ import print_function


import sys
import os
from os import path
from os import listdir
from os.path import exists, join
import subprocess
import argparse

from interpreterParser import *

#######################################################################################################
##  create the interpreter and argument parser objects
#######################################################################################################

iP = interpreterParser()
parser = argparse.ArgumentParser()

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

verboseOutput = False
verbosity = 0

executeRun = False
path_to_charon = ""
writeCurrents = " --current"

#######################################################################################################
##  add command line options 
#######################################################################################################

#parser.add_argument("filename",help="Specify the intpreter input file name.")
parser.add_argument("-i","--input",help="Specify the interpreter input file")
parser.add_argument("-r","--run",help="Execute Charon",action="store_true")
parser.add_argument("--no_current",help="Disable text output of current",action="store_true")
parser.add_argument("--np",help="Specify the number of processors to use in the Charon run.",type=int)
parser.add_argument("-s","--syntax",help="Get help for syntax of interpreter. Append \"| less -i\" to scrolling and search.",action="store_true")
parser.add_argument("-S","--longsyntax",help="Get expanded help for syntax of interpreter. Append \"| less -i\" to scrolling and search.",action="store_true")
parser.add_argument("-v","--verbosity",help="Specify verbose output of the interpreter ranging from 0 to 25",type=int)
parser.add_argument("-p","--path_to_charon",help="Specify path to Charon exectutable")
args = parser.parse_args()

#######################################################################################################
##  process command line options 
#######################################################################################################


if args.run:
    executeRun = True
else:
    executeRun = False

if args.no_current:
    writeCurrents = ""

if args.input != None:
    readFileName = True
    filename = args.input

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
    plFilename = filename+".xml"
    iP.parseInterpreterFile(inputLines, plFilename, verbosity)

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


if executeRun == True:
    try:
        runCommand = "mpirun "
        arguments = " -np "+str(numProcs)+" "+path_to_charon+charon_executable+" --i="+plFilename+writeCurrents
        totalCommand = runCommand+arguments
        print ("executing:  ",runCommand+arguments)
        returnValue = subprocess.check_call(totalCommand,shell=True)

    except subprocess.CalledProcessError as e:
        print ("Interpreter aborting.",e)

