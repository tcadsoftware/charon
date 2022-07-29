#! /usr/bin/env python3

import sys
import os
from os import path
from os import listdir
from ..commonFunctions.fileManager import *
from ..commonFunctions.findCurrentAtParameter import *
from ..commonFunctions.parse.parse import *



class betaGain:
    "calculate the gain either at a specific voltage or from a standalone log"


    def __init__(self):
        self.myName = "betaGain"
        self.myType = "betaGain"
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
        self.voltColumnString = "emitter,voltage"
        self.baseCurrentColumnType = "string"
        self.baseCurrentColumnString = "base,current"
        self.collectorCurrentColumnType = "string"
        self.collectorCurrentColumnString = "collector,current"

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
        dataTypeD = search("data = {dataType:>S}",args)
        filenameD = search("filename = {filename:>S}",args)
        voltageD = search("voltage = {voltage:>g}",args)
        responseNameD = search(" responsename = {responseName:>S}",args)
        contextD = search ("context = {context:>S}",args)
        targetValueD = search("target = {target:>g}",args)
        weightD = search (" weighted = {weighted:>S}",args)
        vColumnD = search(" voltContact = {vColumn:>S}", args)
        bCColumnD = search(" baseCurrentContact = {bCColumn:>S}", args)
        cCColumnD = search(" collectorCurrentContact = {cCColumn:>S}", args)

        if vColumnD != None:
            self.voltColumnType = "string"
            self.voltColumnString = vColumnD['vColumn']
        else:
            self.voltColumnType = "string"
            self.voltColumnString = "emitter,voltage"

        if bCColumnD != None:
            self.baseCurrentColumnType = "string"
            self.baseCurrentColumnString = bCColumnD['bCColumn']+",current"
        else:
            self.baseCurrentColumnType = "string"
            self.baseCurrentColumnString = "base,current"

        if cCColumnD != None:
            self.collectorCurrentColumnType = "string"
            self.collectorCurrentColumnString = cCColumnD['cCColumn']+",current"
        else:
            self.collectorCurrentColumnType = "string"
            self.collectorCurrentColumnString = "collector,current"


        if dataTypeD == None:
            print("Error!  No data type specified in gain response; single or continuous?")
            sys.exit(1)
        else:
            dataType = dataTypeD['dataType']

        if filenameD == None:
            print("Error!  No filename specified in gain response")
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

        if voltageD != None:
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
                print ("Weighting option for betagain must be set to yes or no!") 
                sys.exit(1)

        self.responseDict['dataType'] = dataType
        self.responseDict['filename'] = filename
        if dataType.lower() == "continuous":
            if voltageD == None:
                print("Error! If gain data is continuous, you must specify a voltage at which gain is calculated")
                sys.exit(1)
            self.responseDict['targetVoltage'] = targetVoltage

    def getDict(self):
        return self.responseDict

    def getResponse(self,responseDict=None):  #passing in argument dict
        if responseDict is not None:
            self.responseDict = responseDict
        gain = self.getResponseInternal()
        return gain


    def getResponseInternal(self):
        self.dataType = self.responseDict['dataType']
        self.filename = self.responseDict['filename']
        if self.weighted == True and self.targetIsSet == True:
            self.weight = self.targetValue
        else:
            self.weight = 1.0
        if self.dataType.lower() == "continuous":
            gain = self.getResponseContinuous()
        else:
            gain = self.getResponseSinglePoint()

        return gain/self.weight



    def getResponseSinglePoint(self):
        logFile = list(open(self.filename))
        collectorCurrent = 0.0
        baseCurrent = 0.0
        baseCurrentName = self.baseCurrentColumnString.split(",")
        collectorCurrentName = self.collectorCurrentColumnString.split(",")
        for line in logFile:
            lineTokens = line.split()
            if len(lineTokens) != 0:
                if "emitter" in lineTokens[0].lower() and "current" in lineTokens[0].lower():
                    emitterCurrent = lineTokens[2]
                if baseCurrentName[0] in lineTokens[0].lower() and baseCurrentName[1] in lineTokens[0].lower():
                    baseCurrent = lineTokens[2]
                if collectorCurrentName[0] in lineTokens[0].lower() and collectorCurrentName[1] in lineTokens[0].lower():
                    collectorCurrent = lineTokens[2]

        gain = float(collectorCurrent)/float(baseCurrent)
        returnValue = gain
        if self.targetIsSet:
            returnValue -= self.targetValue

        return returnValue

        

    def getResponseContinuous(self):
        self.targetVoltage = self.responseDict['targetVoltage']
        
        #######################################################################################################
        ##  create the file manager object
        #######################################################################################################
        
        self.fM = fileManager()

        #######################################################################################################
        ##  get the currents
        #######################################################################################################

        success = self.fM.setFileName(self.filename)
        if not success:
            return "FAIL"
        self.fM.setColumns(self.voltColumnType,self.voltColumnString,self.baseCurrentColumnType,self.baseCurrentColumnString)
        (volts,baseCurrent) = self.fM.readFile()

        self.fM.setColumns("string","base","string","collector")
        self.fM.setColumns(self.voltColumnType,self.voltColumnString,self.collectorCurrentColumnType,self.collectorCurrentColumnString)
        (volts,collectorCurrent) = self.fM.readFile()


        #######################################################################################################
        ##  create the find current object
        #######################################################################################################
        
        self.fC = findCurrentAtParameter()
        self.fC.setTargetParameter(self.targetVoltage)

        self.fC.setData(volts,baseCurrent)

        success = False
        (baseCurrent,success) = self.fC.calculateCurrent()

        if not success:
            return "FAIL"

        self.fC.setData(volts,collectorCurrent)
        (collectorCurrent,success) = self.fC.calculateCurrent()

        if not success:
            return "FAIL"

        gain = collectorCurrent/baseCurrent

        returnValue = gain
        if self.targetIsSet:
            returnValue -= self.targetValue

        return returnValue




        #######################################################################################################
        ##  create help
        #######################################################################################################

    def myHelp(self):

        args = []
        details = []

        print ("reponse betaGain argument1=value1 argument2=value2...")
        print ()
        print ("This response returns the beta gain of a transistor.")

        args.append("filename")
        details.append("Name of the file where the currents and voltages are expect to be found.")

        args.append("datatype")
        details.append("<singular>,<continuous> -- if a single steady state solution is calculated, specify singular and the name of the run log for filename and the log will be searched for the currents to calculate the gain.  If a sweep was done, specify continuos and the name of the sweep data file and a voltage at whcih the gain is to be calculated.")

        args.append("voltage")
        details.append("Set a numerical value of the voltage at which the gain is calculated if the data type is continuous.")

        args.append(" target")
        details.append("The target is used for calibration.  It is the expected value of the gain at the specified voltage.  When the target is specified, a residual of the expected and calculated values is returned.")

        args.append("weighted")
        details.append("When the residual is calculated with a target value, setting this to yes or no will opt into weighting the residual by the expected value.  The default is to weight the residual.")

        args.append("responsename")
        details.append("The default name of the response is betaGain.  However, if different gains are to be calculated such as weighted and unweighted or at multiple voltage values, a unique name must be given to each response.")

        args.append("voltColumn")
        details.append("Specify a custom column index for the volts in the I-V data file.  Note that the index starts from 0.  I SAY AGAIN!!  THE FIRST COLUMN OF THE DATA IS INDEX 0. THE SECOND COLUMN IS 1.  And so on.")


        args.append("baseCurrentContact")
        details.append("Specify a custom column header name for the base contact in the I-V data file.  This is the contact name only as specified in Charon input file.")

        args.append("collectorCurrentContact")
        details.append("Specify a custom column header name for the collector contact in the I-V data file.  This is the contact name only as specified in Charon input file.")



        for index,arg in enumerate(args):
            print()
            print (arg,"--",details[index])
            #            print (details[index])
            print ()


        bCColumnD = search(" baseCurrentContact = {bCColumn:>S}", args)
        cCColumnD = search(" collectorCurrentContact = {cCColumn:>S}", args)
