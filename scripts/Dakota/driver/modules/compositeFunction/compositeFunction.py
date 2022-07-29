#! /usr/bin/env python3

import sys
import os
from os import path
from os import listdir
from ..commonFunctions.fileManager import *
from ..commonFunctions.findCurrentAtParameter import *
from ..commonFunctions.parse.parse import *



class compositeFunction:
    "calculate the current at a fiven voltage"


    def __init__(self):
        self.myName = "compositeFunction"
        self.myType = "compositeFunction"
        self.myResponses = []

        self.responseDict = {'method':"NoMethod",'responses':"NoResponses"}

    def getMyType(self):
        return self.myType
        
    def myResponse(self):
        return self.myName
        
    def renameResponse(self,newName):
        self.myName = newName
        
    def createDict(self,respArgs):

        args = " ".join(respArgs)
        args = " ".join(args.replace("="," = ").split())

        methodD = search("method = {method:>S}",args)
        responsesD = search("responses = {responses:>S}",args)
        responseNameD = search(" responsename = {responseName:>S}",args)
 
        if methodD == None:
            print("Error!  No method specified in compositeFunction response")
            sys.exit(1)
        else:
            method = methodD['method']

        if responsesD == None:
            print("Error!  No responses specified in compositeFunction response")
            sys.exit(1)
        else:
            responses = responsesD['responses']
            self.myResponses = responses.split(",")

        if responseNameD != None:
            responseName = responseNameD['responseName']
            self.renameResponse(responseName)


        self.responseDict['method'] = method

    def getDict(self):
        return self.responseDict

    def getResponse(self,responseDict=None):  #passing in argument dict
        if responseDict is not None:
            self.responseDict = responseDict
        compositeValue = self.getResponseInternal()


        return compositeValue


    def getResponseInternal(self):


        returnValue = 0.0

        #For a composite, the return value is a list of the desired component responses
        returnValue =  self.myResponses


        print ()
        print ("-------------------------------------------------------------------")
        print ("Response:   ",self.myResponse())
        print ("-------------------------------------------------------------------")
        print ("value:                     ",str(returnValue))
        print ("-------------------------------------------------------------------")
        print()

        return returnValue




        #######################################################################################################
        ##  create help
        #######################################################################################################

    def myHelp(self):

        args = []
        details = []

        print ("reponse compositeFunction argument1=value1 argument2=value2...")
        print ()
        print ("This response returns the current at a specified voltage from a voltage sweep simulation.")

        args.append("responses")
        details.append("Name of the response functions that comprise this composite.")

        args.append("method")
        details.append("Set the calculation method for the composite: Least Squares. ")

        args.append("responsename")
        details.append("The default name of the response is compositeFunction.  However, if different composites are to be calculated, a unique name must be given to each reponse.")

        
        for index,arg in enumerate(args):
            print()
            print (arg,"--",details[index])
            #            print (details[index])
            print ()
