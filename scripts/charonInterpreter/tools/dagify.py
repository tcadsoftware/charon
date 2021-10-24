#! /usr/bin/env python3
from __future__ import print_function

from .charonDAGify.charonDAGify import *
from .xmlToLCMConverter.xmlLCMConverter import *


#First, convert the xml parameter list to LCMese
x2lcm = xmlLCMConverter("test.xml")
x2lcm.convertFile()
parameterList = x2lcm.getLCMParameters()

#Create the dagifier object
pLDagify = charonDAGify(parameterList,False)
plDAG = pLDagify.createDAG()

