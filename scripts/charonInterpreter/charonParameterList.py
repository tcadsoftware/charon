#! /usr/bin/env python3
from __future__ import print_function

from charonParameter import *
import sys

class charonParameterList:

    def __init__(self, name='NoName'):
        self.lists = []
        self.parameters = []
        self.nestLevel = 0
        self.listName = name
        self.Initiator = name
        self.Terminator = ""
        self.spacing = "  "
        self.satisfied = False
        self.unsatisfiedThings = []
        self.success = False

    def addBlankList(self, nestedPath, newListName,newListInitiator):
        for l in self.lists:
            self.foundList = False
            if l.listName == newListName:
                self.foundList = True

        if self.foundList != True:
            self.localParameterList = charonParameterList(newListName)
            self.localParameterList.Initiator = newListInitiator
            self.lists.append(self.localParameterList)


    def addList(self, nestedPath, newList):
        self.listName = newList.listName
        self.success = False
        # First check if this is where the parameter goes
        # loop over main sub lists
        for l in self.lists:
            if l.listName == nestedPath[0]:
                for ll in newList.lists:
                    l.lists.append(ll)
                for p in newList.parameters:
                    l.parameters.append(p)
                return True

        return(self.success)

    def addSubList(self, newListName, index):
        self.newList = charonParameterList(newListName)
        self.newList.nestLevel = index
        self.lists.append(self.newList)
        return self.lists[len(self.lists)-1]

    def printList(self,indention):
        self.printLine = indention+"<ParameterList name=\""+self.Initiator+"\" >"
        if self.Initiator.find("ANONYMOUS") >= 0:
            self.printLine = indention+"<ParameterList>"
        print(self.printLine,end='\n')
        for p in self.parameters:
            p.printParameter(indention + self.spacing)
        for pL in self.lists:
            pL.printList(indention + self.spacing)
        self.printLine = indention+"</ParameterList>"
        print(self.printLine,end='\n')


    def writeXMLList(self,indention,filehandle):
        self.printLine = indention+"<ParameterList name=\""+self.Initiator+"\" >"
        if self.Initiator.find("ANONYMOUS") >= 0:
            self.printLine = indention+"<ParameterList>"
        filehandle.write(self.printLine+"\n")
        for p in self.parameters:
            self.printLine = p.getParameterString(indention + self.spacing)
            filehandle.write(self.printLine+"\n")
        for pL in self.lists:
            pL.writeXMLList(indention + self.spacing,filehandle)
        self.printLine = indention+"</ParameterList>"
        filehandle.write(self.printLine+"\n")

#
# Insert a parameter into the list
#

    def insertParameter(self,nestedPath,parameterName,parameterType,parameterValue,index):
        (nestedPath,parameterName,parameterType,parameterValue) = self.cleanParameterList(nestedPath,parameterName,parameterType,parameterValue)
        self.index = index
        success = False
        self.parameterName = parameterName
        # First check if this is where the parameter goes
        if nestedPath[len(nestedPath)-1] == self.listName:
            #First check to see if the parameter has already been defined in this list
            replaced = False
            for p in self.parameters:
                if p.parameterName == parameterName:
                    p.replaceParameterValue(parameterType,parameterValue)
                    replaced = True
                    success  = True
            if replaced == False:
                self.localParameter = charonParameter()
                self.localParameter.addMe(parameterName,parameterType,parameterValue)
                self.parameters.append(self.localParameter)
                success = True

        #if not found locally, loop over sublists
        if success != True:
            for pL in self.lists:
                nestedSuccess = pL.insertParameter(nestedPath,parameterName,parameterType,parameterValue,self.index+1)
                if nestedSuccess == True:
                    return(True)

        return(success)


#
# Insert a parameter into the list and create a list if necessary
#

    def insertParameterWithListCreation(self,nestedPath,parameterName,parameterType,parameterValue,index):
        (nestedPath,parameterName,parameterType,parameterValue) = self.cleanParameterList(nestedPath,parameterName,parameterType,parameterValue)
        self.index = index
        self.parameterName = parameterName
        # First check if this is where the parameter goes
        self.parameterInserted = False
        if nestedPath[len(nestedPath)-1] == self.listName and self.nestLevel == len(nestedPath)-1:
#        if nestedPath[len(nestedPath)-1] == self.listName:
            #First check to see if the parameter has already been defined in this list
            replaced = False
            for p in self.parameters:
                if p.parameterName == parameterName:
                    p.replaceParameterValue(parameterType,parameterValue)
                    replaced = True
                    self.parameterInserted  = True
                    #return (self.success)
            if replaced == False:
                self.localParameter = charonParameter()
                self.localParameter.addMe(parameterName,parameterType,parameterValue)
                self.parameters.append(self.localParameter)
                self.parameterInserted = True
                #return (self.success)
        if self.parameterInserted == True:
            return True

        #if not found locally, loop over sublists
        self.nestedSuccess = False
        #if self.success != True:
        if self.parameterInserted != True:
            for pL in self.lists:
                if pL.listName == nestedPath[self.index+1]:
                    self.nestedSuccess = pL.insertParameterWithListCreation(nestedPath,parameterName,parameterType,parameterValue,index+1)
                    if self.nestedSuccess == True:
                        return self.nestedSuccess
            # If the list wasn't found, create it
            # Only get to this line if it wasn't found--no need for if block
            if self.nestedSuccess == False:
                self.newSubList = self.addSubList(nestedPath[index+1],index+1)
                self.nestedSuccess = self.newSubList.insertParameterWithListCreation(nestedPath,parameterName,parameterType,parameterValue,index+1)
                return self.nestedSuccess

        return(self.success)


#
#  Perform an error check of the parameter list
#

    def errorCheckLists(self):
        self.allIsWell = True
        if self.satisfied == True:
            self.allIsWell = True
        for pL in self.lists:
            (self.localCheck,self.returnName) = self.allIsWell = pL.errorCheckLists()
            print("Cjec=== ",self.localCheck,end='\n')
            if self.localCheck == False:
                self.unsatisfiedStuff.append(self.returnName)
        for p in self.parameters:
            (self.localCheck,self.returnName) = p.errorCheckParameter()
            if self.localCheck == False:
                self.unsatisfiedStuff.append(self.returnName)
        print("Returning ",self.satisfied,end='\n')
        return (self.allIsWell,self.listName)


#
# clean the parameter list of leading and trailing whitespace--it screws up Charon
#

    def cleanParameterList(self,nestedPath,parameterName,parameterType,parameterValue):
        for index, value in enumerate(nestedPath):
            nestedPath[index] = value.strip()
        parameterName = parameterName.strip()
        parameterType = parameterType.strip()
        parameterValue = parameterValue.strip()

        return (nestedPath,parameterName,parameterType,parameterValue)
