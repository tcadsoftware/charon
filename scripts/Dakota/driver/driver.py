#! /usr/bin/env python3


import sys
import os
import argparse
import time
import shutil
import math

from configureDriver import *
from preProcessCommand import *
from executeCommand import *
from extractInputParameter import *
from driverHelper import *

############################################################################### 
## Import modules to generate responses
###############################################################################

from modules.thresholdVoltage.thresholdVoltage import *
from modules.currentAtVoltage.currentAtVoltage import *
from modules.betaGain.betaGain import *
from modules.IVCurve.IVCurve import *
from modules.compositeFunction.compositeFunction import *

############################################################################### 
## Take a pause for the cause
###############################################################################

time.sleep(1)



############################################################################### 
## Initialize execute methods list.
###############################################################################

executeMethods = []

############################################################################### 
## Generate a list of all possible responses
###############################################################################

selectedResponses = []
downSelectedResponses = []
availableResponses = []
availableResponses.append(thresholdVoltage())
availableResponses.append(currentAtVoltage())
availableResponses.append(betaGain())
availableResponses.append(IVCurve())
availableResponses.append(compositeFunction())


############################################################################### 
## Create the helper object
###############################################################################

helpObj = driverHelper(availableResponses)


############################################################################### 
## If the argument is listResponses, print the possibles
###############################################################################

if sys.argv[1].lower() == "--listresponses":
    for index,response in enumerate(availableResponses):
        print ("Response ",str(index)," is ",response.myResponse())
    sys.exit(0)


############################################################################### 
## If the argument is listResponses, print the possibles
###############################################################################

if sys.argv[1].lower() == "--help":
    if len(sys.argv) == 2:
        print ("charonDriver help:")
        print ("for a list of available responses (QOIs), charonDriver.py --listresponses")
        print ("For a description of how to use a particular response: charonDriver.py --help <response>")
    else:
        helpObj.getHelp(sys.argv[2])
    sys.exit(0)


############################################################################### 
# Arguments from Dakota are special variables that contain the 
# Dakota parameters and results files, respectively.
############################################################################### 

inputParameters = sys.argv[1]
outputFile = sys.argv[2]

############################################################################### 
#Collect responses responses
############################################################################### 

############################################################################### 
## Save the input parameters in case they're wanted later
###############################################################################

shutil.copy(inputParameters,"pythonParameters.inp")
eIP = extractInputParameter()
evalID = eIP.extract("pythonParameters.inp","eval_id",1,0)


############################################################################### 
# Read the driver.config file to configure driver responses
############################################################################### 

executProcessors = 1
inputTemplate = "NoName"
derivedInput = "NoName"
responseList = []
respArgList = []
calibrationConfig = []
parametersDict = {}

cd = configureDriver()
(parametersDict,executeMethods,responseList,respArgList,calibrationConfig) = cd.configure()
selectedResponses = cd.setResponses(availableResponses)

parametersDict['label'] = evalID


for index,resp in enumerate(selectedResponses):
    print ("EvalID ",evalID," Selected Response ",index,":  ",resp.myResponse())

############################################################################### 
# set the paramter dict for each selected repsonse
############################################################################### 
for index,sr in enumerate(selectedResponses):
    sr.createDict(respArgList[index])

############################################################################### 
# read the dakota input if it's present and make sure responses are consistent
# in both type and order
############################################################################### 

#selectedResponses = cd.readDakotaInputAndOrder(selectedResponses)
downSelectedResponses = cd.validateAndOrderResponses(selectedResponses,inputParameters)

for index,dSR in enumerate(downSelectedResponses):
    print ("EvalID ",evalID,"down selected response ",index,dSR.myResponse())


############################################################################### 
##
## Pre-processing Phase -- Generate/configure an input file for your simulation 
##  by substituting in parameter values from the Dakota parameters file.
##
###############################################################################

print ("PreProcessing")

ppC = preProcessCommand()

for eM in executeMethods:
    preProcessingCommand = ppC.preProcess(inputParameters,eM)
    for prep in preProcessingCommand:
        print ("pre ",prep)
        os.system(prep)


print ("Done processing")

############################################################################### 
##
## Execution Phase -- Run your simulation
##
###############################################################################


print ("Executing")

eC = executeCommand()

for eM in executeMethods:
    executeCommand = eC.execute(parametersDict,eM)
    print ("exec ",executeCommand)
    os.system(executeCommand)

print ("Finished execution")

############################################################################### 
##
## Post-processing Phase -- Extract (or calculate) quantities of interest
##  from your simulation's output and write them to a properly-formatted
##  Dakota results file.
##
###############################################################################
output = open(outputFile,'w')

writeString = ""
responseVector = [0.0]*len(downSelectedResponses)
compSucceeded = True

for index,resp in enumerate(downSelectedResponses):
    print ("Response at evalID ",evalID,":   ",resp.myResponse(),resp.getResponse())
    if resp.getMyType() == "compositeFunction":
        compResp = resp.getResponse()
        compRespVal = 0.0
        for cR in compResp:
            for sR in selectedResponses:
                if cR == sR.myResponse():
                    response = sR.getResponse()
                    if response == "FAIL":
                        compSucceeded = False
                        break
                    compRespVal += float(response)*float(response)
        responseVector[index] = str(math.sqrt(compRespVal))
    else:
        responseVector[index] = str(resp.getResponse())
        if responseVector[index] == "FAIL":
            compSucceeded = False

if compSucceeded == True:
    for rV in responseVector:
        writeString += str(rV)+" "
else:
    writeString = "FAIL"

output.write(writeString)

output.close()

###############################################################################
# Copy results file for convenient access later
###############################################################################
shutil.copy(outputFile,"myParameters.dat")
