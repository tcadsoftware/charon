#! /usr/bin/env python3
from __future__ import print_function

import sys
import _ctypes
from charonDAGify.charonDAGify import *
from charonDAGify.charonDAG import *
from xmlToLCMConverter.xmlLCMConverter import *
from charonDAGify.DAGCompare import *

####################################################
## Set up dereferencing for the lists
####################################################

def di(obj_id):
    """ Inverse of id() function. """
    return _ctypes.PyObj_FromPtr(obj_id)

####################################################
## Set up dereferencing for the lists
####################################################

#if __name__ == '__main__':

if len(sys.argv) != 3:
    print("Usage: compareParameterLists parameterList1 parameterList2")
    sys.exit()

compareParameterValues = True

pLFile1 = sys.argv[1]
pLFile2 = sys.argv[2]

#First, convert the two xml parameter list to LCMese
x2lcmPL1 = xmlLCMConverter(pLFile1)
x2lcmPL1.convertFile()
parameterList1 = x2lcmPL1.getLCMParameters()

x2lcmPL2 = xmlLCMConverter(pLFile2)
x2lcmPL2.convertFile()
parameterList2 = x2lcmPL2.getLCMParameters()

#Create the dagifier object
pLDagifyPL1 = charonDAGify(parameterList1,compareParameterValues)
plDAGPL1 = pLDagifyPL1.createDAG()
#plDAGPL1.printDAG("")

pLDagifyPL2 = charonDAGify(parameterList2,compareParameterValues)
plDAGPL2 = pLDagifyPL2.createDAG()
#plDAGPL2.printDAG("")

#Perform the Forward comparison
dCompForward = DAGCompare(plDAGPL1,plDAGPL2,compareParameterValues)
forwardDiffList = []
forwardDiffList = dCompForward.compareDAGs(forwardDiffList)
if len(forwardDiffList) > 0:
    print()
    print ("The following is a list of items conained in "+pLFile1+", but not in "+pLFile2)
    print()
    forwardDiffList = sorted(forwardDiffList)
    for diff in forwardDiffList:
        print (diff)
    print()

dCompBackward = DAGCompare(plDAGPL2,plDAGPL1,compareParameterValues)
backwardDiffList = []
dCompBackward.compareDAGs(backwardDiffList)

if len(backwardDiffList) > 0:
    print()
    print()
    print ("The following is a list of items conained in "+pLFile2+", but not in "+pLFile1)
    print()
    backwardDiffList = sorted(backwardDiffList)
    for diff in backwardDiffList:
        print (diff)
    print ()

if  len(forwardDiffList) == 0 and len(backwardDiffList) == 0:
    print("Files are the same.")
else:
    print("Files are different.")
