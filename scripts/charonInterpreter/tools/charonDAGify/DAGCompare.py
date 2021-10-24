#! /usr/bin/env python3
from __future__ import print_function

import sys
from .charonDAG import *

class DAGCompare:
    "Compare two DAGS to see if they are equivalent"

    def __init__(self,dag1,dag2,compareParameterValues):
        self.compareParameterValues = compareParameterValues
        self.dag1 = dag1
        self.dag2 = dag2

    ####################################################
    ## Set up dereferencing for the lists
    ####################################################

    def di(self,obj_id):
        """ Inverse of id() function. """
        return _ctypes.PyObj_FromPtr(obj_id)



    def compareDAGs(self,diffList):
        #Compare dags

        if self.dag1.getChildList() != self.dag2.getChildList():
            for child in list(set(self.dag1.getChildList()) - set(self.dag2.getChildList())):
                diffList.append(child.convertToString(""))

        for child in list(set(self.dag1.getChildList()) - (set(self.dag1.getChildList()) - set(self.dag2.getChildList()))):
            newCompare = DAGCompare(child, self.dag2.getChildObject(child.getName()),self.compareParameterValues )
            diffList = newCompare.compareDAGs(diffList)

        # send the diff back
        return diffList
