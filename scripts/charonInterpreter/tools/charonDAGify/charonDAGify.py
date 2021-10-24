#! /usr/bin/env python3
from __future__ import print_function
from .charonDAG import *

class charonDAGify:
    "Create a poor-man's DAG of the parameter list"

    def __init__(self,parameterList,includeParameterValues=False):
        self.nodes = []
        self.myNodeName = ""
        self.retainParameterValues = False
        self.myDAG = charonDAG(None)  # This is the root and has no parent
        self.parameterList = parameterList
        self.includeParameterValues = includeParameterValues


#############################################################################
#  create and return the DAG
#############################################################################

    def createDAG(self):
        rootSet = False
        for pL in self.parameterList:
            pLList,pLDelimiters = self.processForInsertion(pL)
            if rootSet == False:
                self.myDAG.setName(pLList[0])
                self.myDAG.delimiter = ""
                rootSet = True
            self.myDAG.insertPL(pLList,pLDelimiters)

        return self.myDAG


    #############################################################################
    #  process PL into a list for insertion
    #############################################################################
    def processForInsertion(self,pL):
        entryList = pL.split("->")
        DAGDelimiterList = []
        parameterSpecifics = entryList[-1].split(",")
        DAGEntryList = entryList[:-1]
        for dg in DAGEntryList:
            DAGDelimiterList.append("->")

        if self.includeParameterValues == True:
            DAGEntryList.append(parameterSpecifics[0]) #Terminal node before parameter parameter specifics
            DAGDelimiterList.append("->")
            DAGEntryList.append(parameterSpecifics[1]) #parameter name
            DAGDelimiterList.append(",")
            DAGEntryList.append(parameterSpecifics[2]) #parameter type
            DAGDelimiterList.append(",")
            DAGEntryList.append(parameterSpecifics[3]) #parameter value
            DAGDelimiterList.append(",")
        else:
            DAGEntryList.append(parameterSpecifics[0]) #Terminal node before parameter parameter specifics
            DAGDelimiterList.append(",")
            DAGEntryList.append(parameterSpecifics[1]) #parameter name
            DAGDelimiterList.append(",")

        return DAGEntryList,DAGDelimiterList

