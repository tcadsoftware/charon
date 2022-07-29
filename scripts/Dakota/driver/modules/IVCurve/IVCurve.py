#! /usr/bin/env python3

import sys
import os
import math
from os import path
from os import listdir
from ..commonFunctions.fileManager import *
from ..commonFunctions.findCurrentAtParameter import *
from ..commonFunctions.parse.parse import *




class IVCurve:
    "calculate the current at a fiven voltage"


    def __init__(self):
        self.myName = "IVCurve"
        self.myType = "IVCurve"
        self.readFilename = False
        self.filename = "NoName"
        self.coordinates = "NoName"

        self.readCurrent = False
        self.targetVoltage = 1.0  # Default quantity

        self.voltColumnFlag = True  # This will be column 0 almost every time
        self.voltColumn = 0
        self.currentColumnFlag = False
        self.currentColumn = -1
        self.responseDict = {'filename':"NoName",'targetVoltage':1.0}
        self.targetIsSet = False
        self.targetFile = ""
        self.weight = 1.0
        self.weighted=False
        self.voltColumnType = "string"
        self.voltColumnString = "gate,voltage"
        self.currentColumnType = "string"
        self.currentColumnString = "source,current"

    def myResponse(self):
        return self.myName
        
    def getMyType(self):
        return self.myType
        
    def myContext(self):
        return self.context

    def renameResponse(self,newName):
        self.myName = newName
        
    def createDict(self,respArgs):
        args = " ".join(respArgs)
        args = " ".join(args.replace("="," = ").split())
        filenameD = search("filename = {filename:>S}",args)
        responseNameD = search(" responsename = {responseName:>S}",args)
        coordinatesD = search(" coordinates = {coordinates:>S}",args)
        targetFileD = search(" targetFile = {targetFile:>S}",args)
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

        if responseNameD != None:
            responseName = responseNameD['responseName']
            self.renameResponse(responseName)

        if coordinatesD == None:
            print("Error!  No coordinates file specified in IVCurve response")
            sys.exit(1)
        else:
            coordsFile = coordinatesD['coordinates']

 
        if targetFileD != None:
            self.targetFile = targetFileD['targetFile']
            self.targetIsSet = True

        if weightD != None:
            weighting = weightD['weighted']
            if weighting.lower() == "yes":
                self.weighted = True
            elif weighting.lower() == "no":
                self.weighted = False
            else:
                print ("Weighting option for the I-V curve must be set to yes or no!") 
                sys.exit(1)

        self.responseDict['filename'] = filename
        self.responseDict['coordinates'] = coordsFile


    def getDict(self):
        return self.responseDict

    def getResponse(self,responseDict=None):  #passing in argument dict
        if responseDict is not None:
            self.responseDict = responseDict
        currentAtV = self.getResponseInternal()
        return currentAtV


    def getResponseInternal(self):
        self.filename = self.responseDict['filename']

        self.coordinates = self.responseDict['coordinates']

        #######################################################################################################
        ##  Read in the coordinates 
        #######################################################################################################
        
        coordinates = list(open(self.coordinates))


        #######################################################################################################
        ##  Read in the target values if specified
        #######################################################################################################
        
        if self.targetIsSet == True:
            targets = list(open(self.targetFile))


        self.weight = [1.0]*len(coordinates)

        if self.weighted == True and self.targetIsSet == True:
            self.weight = targets

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

        currents = []
        residuals = []
        for index,volts in enumerate(coordinates):
            self.fC.setTargetParameter(float(volts))
            success = False
            (currentAtV,success) = self.fC.calculateCurrent()
            if not success:
                return "FAIL"
            else:
                if self.targetIsSet == True:
                    residual = (float(currentAtV) - float(targets[index]))/float(self.weight[index])
                    residuals.append(residual)
                else:
                    currents.append(str(currentAtV))

        if self.targetIsSet == True:
            for res in residuals:
                currents.append(str(res))

        returnValue = ""
        for curr in currents:
            returnValue += str(curr)+"\n"

        return returnValue




        #######################################################################################################
        ##  create help
        #######################################################################################################

        def myHelp(self):

            args = []
            details = []

            print ("reponse IVCurve argument1=value1 argument2=value2...")
        print ()
        print ("This response returns the currents at a series of specified voltage from a voltage sweep simulation.  See Dakota documentation for details.")

        args.append("filename")
        details.append("Name of the file where the currents and voltages are expect to be found.")

        args.append("voltage")
        details.append("Set a numerical value of the voltage at which the current is calculated.")

        args.append(" targetFile")
        details.append("The targetFile is used for calibration.  It is the name of a file that has expected values of the current at specified voltages.  When the target file is specified, a set of residuals of the expected and calculated values is returned.")

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
