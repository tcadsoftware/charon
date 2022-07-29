#! /usr/bin/env python3
#from __future__ import print_function

import sys
import os
from os import path
from os import listdir

from createLineParser import *
from createBlockParser import *
from createParserLibrary import *
from createModifiers import *

verbosity = 0
if len(sys.argv) > 1:
    for opt in range(len(sys.argv)-1):
        if sys.argv[opt+1][:11] == "--verbosity":
            optionTokens = sys.argv[opt+1].split("=")
            verbosity = optionTokens[1]

inputFiles = []

sourcePath = "parseInputs"
targetPath = "../parsers"
modifierTargetPath = "../modifiers"
libName = "parsers"

localBlockParserList = []
localBlockParsingKeyList = []
localParserList = []
localParsingKeyList = []

#Create the parsers directory
if path.isdir(targetPath) == False:
    os.mkdir(targetPath)
#Create the modifiers directory
if path.isdir(modifierTargetPath) == False:
    os.mkdir(modifierTargetPath)

for filename in listdir(sourcePath):
    inputFiles.append(filename)

for inputs in inputFiles:
    if inputs[-4:].lower() == ".inp":
        cP = createLineParser(inputs,sourcePath,targetPath,verbosity)
        cP.parseInputs()
        localParserList.append(cP.getParserName())
        localParsingKeyList.append(cP.getParsingKey())
        cP.createLineParserSourceFile()
        del cP
    elif inputs[-9:].lower() == ".blockinp":
        cBP = createBlockParser(inputs,sourcePath,targetPath,verbosity)
        cBP.parseInputs()
        localBlockParserList.append(cBP.getBlockParserName())
        localBlockParsingKeyList.append(cBP.getBlockParsingKey())
        cBP.createBlockParserSourceFile()
        del cBP

############################################################################
# Check for duplicate parser names and parsing keywords in this scope
############################################################################

for index,lpl in enumerate(localParserList):
    for lpl2 in localParserList[index+1:]:
        if lpl == lpl2:
            print ("ERROR!! There are duplicate parser names in the <<",sourcePath,">> scope")
            print ("Duplicated name is ",lpl,"\n")

############################################################################
# Check for duplicate block parser names and block parsing keywords in this scope
############################################################################

for index,lbpl in enumerate(localBlockParserList):
    for lbpl2 in localBlockParserList[index+1:]:
        if lbpl == lbpl2:
            print ("ERROR!! There are duplicate BLOCK parser names in the <<",sourcePath,">> scope")
            print ("Duplicated name is ",lbpl,"\n")

for index,lbpkl in enumerate(localBlockParsingKeyList):
    for lbpkl2 in localBlockParsingKeyList[index+1:]:
        if lbpkl == lbpkl2:
            print ("ERROR!! There are duplicate BLOCK parsing keys in the <<",sourcePath,">> scope")
            print ("Duplicated parsing key is ",lbpkl,"\n")


#Create the parser libraries in the primary parsers space
cPL = createParserLibrary()
cPL.writeParserLibraryClass(targetPath,libName,localParserList,localBlockParserList)

#Create the modifier library in the modifier space
cM = createModifiers(sourcePath,modifierTargetPath,verbosity)
cM.generateModifiers()
cM.createModifierLibrary()

