#! /usr/bin/env python3
from __future__ import print_function

class charonDAG:
    "Poor man's charon DAG"


    def __init__(self,parent):
        self.childNodes = []
        self.myValue = ""
        self.myType = ""
        self.myNodeName = ""
        self.myIndent = "  "
        self.myParent = parent
        self.delimiter = ""

    #############################################################################
    #  Create an equal operator
    #############################################################################
    def __eq__(self,other):
        return isinstance(other,self.__class__) and self.getName() == other.getName()


    #############################################################################
    #  Create a not equal operator
    #############################################################################
    def __ne__(self, other):
        return not self.__eq__(other)

    #############################################################################
    #  Create a subtraction operator
    #############################################################################
    def __sub__(self,other):
        return self.__class__(*[item for item in self if item not in other]).getName()

    #############################################################################
    #  Create a hash
    #############################################################################
    def __hash__(self):
        return hash((self.getName(), self.getName()))

    #############################################################################
    #  Overload the string for printing
    #############################################################################
    def __str__(self):
        if self.myParent != None:
            print (self.myParent,end="")
        print(self.delimiter+self.getName())
        return (self.delimiter+self.getName())

    #############################################################################
    #  Convert the DAG into a string for sorted output
    #############################################################################
    def convertToString(self,myString):
        myGlobalString = ""
        if self.myParent != None:
            myGlobalString = self.myParent.convertToString(self.delimiter+self.getName()+myString)
        else:
            myGlobalString = self.getName()+myString

        return myGlobalString

      #############################################################################
    #  get this DAG node name
    #############################################################################
    def getName(self):
        return self.myNodeName

    #############################################################################
    #  get this DAG node name
    #############################################################################
    def setName(self,name):
        self.myNodeName = name

    #############################################################################
    #  Insert a parameter list parameter
    #############################################################################
    def insertPL(self,pLList,pLDelimiters):
        self.insertPLWithCreation(pLList,pLDelimiters,1)

    #############################################################################
    #  Insert a parameter list parameter
    #############################################################################
    def insertPLWithCreation(self,pLList,pLDelimiters,index):
        foundNode = False
        for child in self.childNodes:
            if pLList[index] == child.getName():
                ##check the next level
                foundNode = True
                child.insertPLWithCreation(pLList,pLDelimiters,index+1)
                continue

        if index == len(pLList):
            return

        if foundNode == False:
            self.childNodes.append(charonDAG(self))
            self.childNodes[-1].setName(pLList[index])
            self.childNodes[-1].delimiter = pLDelimiters[index]
            self.childNodes[-1].insertPLWithCreation(pLList,pLDelimiters,index+1)

    #############################################################################
    #  print the DAG
    #############################################################################
    def printDAG(self,indent):
        print (indent+self.myNodeName)
        for cN in self.childNodes:
            cN.printDAG(indent+self.myIndent)

    #############################################################################
    #  Check if a particular node exists at this level
    #############################################################################
    def DoesThisNodeExist(self,nodeName):
        for localNodes in self.childNodes:
            if localNodes == nodeName:
                return True

        return False


    #############################################################################
    #  return the child list of a specific index
    #############################################################################
    def getChildList(self):
        return self.childNodes

    #############################################################################
    #  return the child object with a requested name
    #############################################################################
    def getChildObject(self,name):
        for child in self.childNodes:
            if child.getName() == name:
                return child

        return None

