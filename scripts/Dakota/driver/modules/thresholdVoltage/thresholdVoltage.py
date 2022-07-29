#! /usr/bin/env python3

import sys
import os
from os import path
from os import listdir
from ..commonFunctions.fileManager import *
from ..commonFunctions.findParameterAtLogCurrent import *
from ..commonFunctions.parse.parse import *



class thresholdVoltage:
    "calculate threshold voltage based on tabulated current"


    def __init__(self):
        self.myName = "thresholdVoltage"
        self.myType = "thresholdVoltage"
        self.readFilename = False
        self.filename = "NoName"

        self.readCurrent = False
        self.targetCurrent = 0.01  # Default quantity

        self.voltColumnFlag = True  # This will be column 0 almost every time
        self.voltColumn = 0
        self.currentColumnFlag = False
        self.currentColumn = -1
        self.responseDict = {'filename':"NoName",'targetCurrent':0.001,'standalone':False}
        self.UseSlope = False
        self.context = "singular"
        self.targetIsSet = False
        self.targetValue = 0.0
        self.icap = 1e10
        self.weight = 1.0
        self.weighted=False
        self.voltColumnType = "string"
        self.voltColumnString = "gate,voltage"
        self.currentColumnType = "string"
        self.currentColumnString = "drain,current"

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
        currentD = search("current = {current:>g}",args)
        methodD = search("method = {method:>S}",args)
        responseNameD = search("responsename = {responseName:>S}",args)
        contextD = search (" context = {context:>S}",args)
        targetValueD = search(" target = {target:>g}",args)
        standaloneD = search (" standalone = {standalone:>S}",args)
        weightD = search (" weighted = {weighted:>S}",args)
        icapD = search(" icap = {icap:>g}",args)
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
            self.currentColumnString = "drain,current"


        if filenameD == None:
            print("Error!  No filename specified in thresholdVoltage response")
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


        if currentD == None and methodD == None:
            print("Error!  Must specify either current or method  in thresholdVoltage response")
            sys.exit(1)
        elif currentD != None and methodD != None:
            print("Error!  Must specify either current or method BUT NOT BOTH  in thresholdVoltage response")
            sys.exit(1)
        else:
            if currentD != None:
                targetCurrent = currentD['current']
            if methodD != None:
                if methodD['method'].lower() == "slope":
                    self.UseSlope = True
                else:
                    print (methodD['method']," is not a valid method in thresholdVoltage response")
                    sys.exit(1)
 
        if responseNameD != None:
            responseName = responseNameD['responseName']
            self.renameResponse(responseName)

        if standaloneD != None:
            if standaloneD['standalone'] == "True":
                self.responseDict['standalone'] = True

        if weightD != None:
            weighting = weightD['weighted']
            if weighting.lower() == "yes":
                self.weighted = True
            elif weighting.lower() == "no":
                self.weighted = False
            else:
                print ("Weighting option for threshold voltage must be set to yes or no!") 
                sys.exit(1)

        if icapD != None:
            self.icap = float(icapD['icap'])

        self.responseDict['filename'] = filename

        if not self.UseSlope:
            self.responseDict['targetCurrent'] = targetCurrent

    def getDict(self):
        return self.responseDict

    def getResponse(self,responseDict=None):  #passing in argument dict
        if responseDict is not None:
            self.responseDict = responseDict
        vT = self.getResponseInternal()
        return vT


    def getResponseInternal(self):
        self.filename = self.responseDict['filename']
        standalone = self.responseDict['standalone']
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

        #self.fM.setColumns("gate,voltage","drain,current")
        self.fM.setColumns(self.voltColumnType,self.voltColumnString,self.currentColumnType,self.currentColumnString)

        (volts,current) = self.fM.readFile()

        if current[-1] < 0:
            for index,curr in enumerate(current):
                current[index] = -curr
            
        #######################################################################################################
        ##  create the find voltage at current object
        #######################################################################################################
        if not self.UseSlope:
            self.targetCurrent = self.responseDict['targetCurrent']
            self.fT = findParameterAtLogCurrent()

            self.fT.setData(volts,current)
            self.fT.setTargetCurrent(self.targetCurrent)

            success = False
            (vT,success) = self.fT.calculateParameter()

            if success == False:
                return "FAIL"

            returnValueF = vT
            if self.targetIsSet:
                returnValueF = (returnValueF - self.targetValue)/self.weight
 
            returnValue = str(returnValueF) +"\n"
            if standalone == True:
                returnValue += str(self.targetCurrent)+"\n"


        if self.UseSlope:

            slopes = []
            if len(volts) != len(current):
                print ("Something has gone haywire")

            if self.icap < current[0]:
                print ("WARNING: icap too small.  Taking first sequence. Vt may not be accurate.")
                self.icap = current[1]
                       

            for index,vl in enumerate(volts[:-1]):
                if current[index] > self.icap:
                    continue
                mySlope = (current[index+1] - current[index]) / (volts[index+1] - volts[index])
                slopes.append(mySlope)

            maxSlope = max(slopes)
            maxSlopeIdx = slopes.index(maxSlope)

            currentIntercept = current[maxSlopeIdx] - maxSlope*volts[maxSlopeIdx]

            vT = -currentIntercept/maxSlope

            returnValueF = vT
            if self.targetIsSet:
                returnValueF = (returnValueF - self.targetValue)/self.weight

            returnValue = str(returnValueF)+"\n"
            if standalone == True:
                returnValue += str(current[maxSlopeIdx])+"\n"
 
        if self.targetIsSet == False:
            targetOutput = "N/A"
            residualOutput = "N/A"
            rawResidualOutput  = "N/A"
        else:
            targetOutput = str(self.targetValue)
 
        print ()
        print ("-------------------------------------------------------------------")
        print ("Response:   ",self.myResponse())
        print ("-------------------------------------------------------------------")
        print ("value:                     ",str(vT))
        print ("target:                    ",targetOutput)
        print ("-------------------------------------------------------------------")
        print ()

        return returnValue



        #######################################################################################################
        ##  create help
        #######################################################################################################

    def myHelp(self):

        args = []
        details = []

        print ()
        print ("reponse thresholdVoltage argument1=value1 argument2=value2...")
        print ()
        print ("This response returns the threshold voltage.")


        args.append("method")
        details.append("<slope>,<current>   If the method is set to slope, the function will calculate the slope of the IV curve and use it's maximum to interpolate the threshold voltage.  If method is not specified or is specified as current, a current must be specified.")

        args.append("current")
        details.append("Set a numerical value of the contact current at which the threshold voltage is calculated.  This is an alternative to the slope method.")

        args.append("filename")
        details.append("Name of the file where the currents and voltages are expect to be found.")

#        args.append("context")
 #       details.append("")

        args.append(" target")
        details.append("The target is used for calibration.  It is the expected value of the threshold voltage.  When the target is specified, a residual of the expected and calculated values is returned.")

        args.append(" icap")
        details.append("Numerical cap on the current at which the slope method can be used to calculate the threshold voltage")

        args.append("weighted")
        details.append("When the residual is calculated with a target value, setting this to yes or no will opt into weighting the residual by the expected value.  The default is to weight the residual.")

        args.append("responsename")
        details.append("The default name of the response is thresholdVoltage.  However, if differnt threshold voltages are to be calculated such as weighted and unweighted or at multiple current values, a unique name must be given to each reponse.")

        
        args.append("voltColumn")
        details.append("Specify a custom column index for the volts in the I-V data file.  Note that the index starts from 0.  I SAY AGAIN!!  THE FIRST COLUMN OF THE DATA IS INDEX 0. THE SECOND COLUMN IS 1.  And so on.")

        args.append("currentColumn")
        details.append("Specify a custom column index for the current in the I-V data file.  Note that the index starts from 0.  I SAY AGAIN!!  THE FIRST COLUMN OF THE DATA IS INDEX 0. THE SECOND COLUMN IS 1.  And so on.")

        
  
        for index,arg in enumerate(args):
            print()
            print (arg,"--",details[index])
#            print (details[index])
            print ()
