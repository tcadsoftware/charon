#! /usr/bin/env python3

import sys
import os
from os import path
from os import listdir
from ..commonFunctions.fileManager import *
from ..commonFunctions.findCurrentAtParameter import *
from ..commonFunctions.parse.parse import *



class currentAtVoltage:
    "calculate the current at a fiven voltage"


    def __init__(self):
        self.myName = "currentAtVoltage"
        self.myType = "currentAtVoltage"
        self.readFilename = False
        self.filename = "NoName"

        self.readCurrent = False
        self.targetVoltage = 1.0  # Default quantity

        self.voltColumnFlag = True  # This will be column 0 almost every time
        self.voltColumn = 0
        self.currentColumnFlag = False
        self.currentColumn = -1
        self.responseDict = {'filename':"NoName",'targetVoltage':1.0}
        self.context = "singular"
        self.targetIsSet = False
        self.targetValue = 0.0
        self.weight = 1.0
        self.weighted=False
        self.voltColumnType = "string"
        self.voltColumnString = "gate,voltage"
        self.currentColumnType = "string"
        self.currentColumnString = "source,current"

    def getMyType(self):
        return self.myType
        
    def myResponse(self):
        return self.myName
        
    def myContext(self):
        return self.context

    def renameResponse(self,newName):
        self.myName = newName
        
    def createDict(self,respArgs):
        args = " ".join(respArgs)
        args = " ".join(args.replace("="," = ").split())
        filenameD = search("filename = {filename:>S}",args)
        voltageD = search("voltage = {voltage:>g}",args)
        responseNameD = search(" responsename = {responseName:>S}",args)
        contextD = search ("context = {context:>S}",args)
        targetValueD = search("target = {target:>g}",args)
        weightD = search (" weighted = {weighted:>S}",args)
        vColumnD = search(" voltColumn = {vColumn:>S}", args)
        cColumnD = search(" currentColumn = {cColumn:>S}", args)


        if vColumnD != None:
            self.voltColumnType = "index"
            self.voltColumnString = vColumnD['vColumn']
        else:
            self.voltColumnType = "string"
            self.voltColumnString = "gate,voltage"

        if cColumnD != None:
            self.currentColumnType = "index"
            self.currentColumnString = cColumnD['cColumn']
        else:
            self.currentColumnType = "string"
            self.currentColumnString = "source,current"

        if filenameD == None:
            print("Error!  No filename specified in currentAtVoltage response")
            sys.exit(1)
        else:
            filename = filenameD['filename']

        if contextD != None:
            localContext = contextD['context']
            if localContext.lower() != "singular" and localContext.lower() != "compound":
                print("Error! The context in ",self.myResponse()," must be singular or compound")
                sys.exit(1)
            else:
                self.context = localContext.lower()

        if targetValueD != None:
            self.targetValue = targetValueD['target']
            self.targetIsSet = True

        if voltageD == None:
            print("Error!  No voltage specified in thresholdVoltage response")
            sys.exit(1)
        else:
            targetVoltage = voltageD['voltage']

        if responseNameD != None:
            responseName = responseNameD['responseName']
            self.renameResponse(responseName)

        if weightD != None:
            weighting = weightD['weighted']
            if weighting.lower() == "yes":
                self.weighted = True
            elif weighting.lower() == "no":
                self.weighted = False
            else:
                print ("Weighting option for threshold voltage must be set to yes or no!") 
                sys.exit(1)

        self.responseDict['filename'] = filename
        self.responseDict['targetVoltage'] = targetVoltage

    def getDict(self):
        return self.responseDict

    def getResponse(self,responseDict=None):  #passing in argument dict
        if responseDict is not None:
            self.responseDict = responseDict
        currentAtV = self.getResponseInternal()
        return currentAtV


    def getResponseInternal(self):
        self.filename = self.responseDict['filename']
        self.targetVoltage = self.responseDict['targetVoltage']
        if self.weighted == True:
            self.weight = self.targetValue
        else:
            self.weight = 1.0

        #######################################################################################################
        ##  create the file manager object
        #######################################################################################################
        
        self.fM = fileManager()

        #######################################################################################################
        ##  set the response
        #######################################################################################################

        success = self.fM.setFileName(self.filename)
        if not success:
            return "FAIL"


        #self.fM.setColumns("gate,voltage","source,current")
        self.fM.setColumns(self.voltColumnType,self.voltColumnString,self.currentColumnType,self.currentColumnString)


        (volts,current) = self.fM.readFile()

        #######################################################################################################
        ##  create the find threshold object
        #######################################################################################################
        
        self.fC = findCurrentAtParameter()

        self.fC.setData(volts,current)
        self.fC.setTargetParameter(self.targetVoltage)

        success = False
        (currentAtV,success) = self.fC.calculateCurrent()

        if not success:
            return "FAIL"

        returnValue = currentAtV
        if self.targetIsSet:
            returnValue = (returnValue - self.targetValue)/self.weight

        if self.targetIsSet == False:
            targetOutput = "N/A"
            residualOutput = "N/A"
            rawResidualOutput = "N/A"
        else:
            targetOutput = str(self.targetValue)

        print ()
        print ("-------------------------------------------------------------------")
        print ("Response:   ",self.myResponse())
        print ("-------------------------------------------------------------------")
        print ("value:                     ",str(currentAtV))
        print ("target:                    ",targetOutput)
        print ("-------------------------------------------------------------------")
        print()

        return returnValue




        #######################################################################################################
        ##  create help
        #######################################################################################################

    def myHelp(self):

        args = []
        details = []

        print ("reponse currentAtVoltage argument1=value1 argument2=value2...")
        print ()
        print ("This response returns the current at a specified voltage from a voltage sweep simulation.")

        args.append("filename")
        details.append("Name of the file where the currents and voltages are expect to be found.")

        args.append("voltage")
        details.append("Set a numerical value of the voltage at which the current is calculated.")

        args.append(" target")
        details.append("The target is used for calibration.  It is the expected value of the current at the specified voltage.  When the target is specified, a residual of the expected and calculated values is returned.")

        args.append("weighted")
        details.append("When the residual is calculated with a target value, setting this to yes or no will opt into weighting the residual by the expected value.  The default is NOT to weight the residual.")

        args.append("responsename")
        details.append("The default name of the response is currentAtVoltage.  However, if different currents are to be calculated such as weighted and unweighted or at multiple voltages values, a unique name must be given to each reponse.")

        args.append("voltColumn")
        details.append("Specify a custom column index for the volts in the I-V data file.  Note that the index starts from 0.  I SAY AGAIN!!  THE FIRST COLUMN OF THE DATA IS INDEX 0. THE SECOND COLUMN IS 1.  And so on.")

        args.append("currentColumn")
        details.append("Specify a custom column index for the current in the I-V data file.  Note that the index starts from 0.  I SAY AGAIN!!  THE FIRST COLUMN OF THE DATA IS INDEX 0. THE SECOND COLUMN IS 1.  And so on.")

        
        for index,arg in enumerate(args):
            print()
            print (arg,"--",details[index])
            print ()
