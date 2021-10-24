#! /usr/bin/env python3
from __future__ import print_function

import sys
import os
from os import path
from os import listdir

from createLineParser import *
from createParserLibrary import *

class createBlockParser:
    "create block parser from interpreter input and xml"



#######################################################################################################
##  Create block parser constructor
#######################################################################################################

    def __init__(self,filename,sourcePathName,targetPathName,verbosity):
        self.blockFileName = filename
        self.parserBlockName = "noName"
        self.targetPath = targetPathName
        self.sourcePath = sourcePathName
        self.toOpen = self.sourcePath+"/"+self.blockFileName
        self.lines = list(open(self.toOpen))
        self.interpreterBlockHelpLine = "No help section exists"
        self.parsingBlockKey = "Undefined"

        self.verbosity = int(verbosity)
        ###############################################################################################
        # create lists to manage sub block parsers and xml default lines
        ###############################################################################################
        self.xmlDefaultLines = []
        self.xmlExtras = []
        self.xmlExtraValues = []

        if self.verbosity > 0:
            class_name = self.__class__.__name__
            print (class_name, self.blockFileName," created")

#######################################################################################################
##  Destruct block parser destructor
#######################################################################################################

    def __del__(self):
        if self.verbosity > 3:
            class_name = self.__class__.__name__
            print (class_name, self.blockFileName," destroyed")

#######################################################################################################
##  Return the block parser name
#######################################################################################################

    def getBlockParserName(self):
        return self.parserBlockName;

#######################################################################################################
##  Return the block parsing key
#######################################################################################################

    def getBlockParsingKey(self):
        return self.parsingBlockKey;

#######################################################################################################
##  parseInputs:  parses the input file that defines the new block parser
#######################################################################################################

    def parseInputs(self):

        #Read in the interpreter lines
        for line in self.lines:
            lineTokens = line.split()
            foundLineToProcess = False
            if len(lineTokens) != 0 and lineTokens[0][:1] != "#":
                if lineTokens[0].lower() == "interpreterblock":
                    if lineTokens[1].lower() == "name":
                        self.parserBlockName = lineTokens[2]
                        foundLineToProcess = True
                    else:
                        self.parsingBlockKey = line.partition('(')[-1].rpartition(')')[0].lower()
                        if self.parsingBlockKey == "":
                             print("Error!  Cannot extract a block parsing key for \n",line," in ",self.blockFileName,". Check the input file.")
                        foundLineToProcess = True

                        #extract the arguments
                        (self.blockArgument,self.blockArgumentIndexes) = self.extractArguments(line)
                        #Get help line
                        self.interpreterBlockHelpLine = ""
                        for lT in lineTokens[1:]:
                            self.interpreterBlockHelpLine += lT+" "
                        # Remove the trailing line break
                        self.interpreterBlockHelpLine = self.interpreterBlockHelpLine.rstrip('\n')
                        # Strip out the ()s
                        self.interpreterBlockHelpLine = self.interpreterBlockHelpLine.replace(")","").replace("(","")

            #Read in the xml lines
            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "xmlrequired":
                    self.xmlRequiredLines.append((line.replace(lineTokens[0]+" ","")).rstrip())
                    foundLineToProcess = True
            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "xmldefault":
                    self.xmlDefaultLines.append((line.replace(lineTokens[0]+" ","")).rstrip())
                    foundLineToProcess = True

            if len(lineTokens) != 0:
                if lineTokens[0][0] != "#":
                    if foundLineToProcess == False:
                        print ("\nError:  I found nothing to process on this line in ",self.blockFileName,":\n",line," generateInterpreter is exiting.\n")
                        exit(1)


        ###############################################################################################
        #  Create a parser subdirectory where line parsers associated with this block parser live
        ###############################################################################################
        self.blockParserSubDir = self.targetPath + "/" + self.parserBlockName
        if path.isdir(self.blockParserSubDir) == False:
            os.mkdir(self.blockParserSubDir)
 
        ###############################################################################################
        #  Create line parsers that belong in this block subdirectory
        ###############################################################################################
        inputFiles = []
        subLineSource = self.sourcePath+"/"+self.parserBlockName
        subLineTarget = self.targetPath+"/"+self.parserBlockName
        for filename in listdir(subLineSource):
            inputFiles.append(filename)

        self.localParserList = []
        self.localBlockParserList = []
        self.localBlockParsingKeyList = []
        self.localParsingKeyList = []
        for inputs in inputFiles:
            if inputs[len(inputs)-4:].lower() == ".inp":
                cP = createLineParser(inputs,subLineSource,subLineTarget,self.verbosity)
                cP.parseInputs()
                self.localParserList.append(cP.getParserName())
                self.localParsingKeyList.append(cP.getParsingKey())
                cP.createLineParserSourceFile()
                del cP
            elif inputs[len(inputs)-9:].lower() == ".blockinp":
                cBP = createBlockParser(inputs,subLineSource,subLineTarget,self.verbosity)
                cBP.parseInputs()
                self.localBlockParserList.append(cBP.getBlockParserName())
                self.localBlockParsingKeyList.append(cBP.getBlockParsingKey())
                cBP.createBlockParserSourceFile()
                del cBP

        ############################################################################
        # Check for duplicate parser names and parsing keywords in this scope
        ############################################################################

        for index,lpl in enumerate(self.localParserList):
            for lpl2 in self.localParserList[index+1:]:
                if lpl == lpl2:
                    print ("ERROR!! There are duplicate parser names in the <<",subLineSource,">> scope")
                    print ("Duplicated name is ",lpl,"\n")

        for index,lpkl in enumerate(self.localParsingKeyList):
            for lpkl2 in self.localParsingKeyList[index+1:]:
                if lpkl == lpkl2:
                    print ("ERROR!! There are duplicate parsing keys in the <<",subLineSource,">> scope")
                    print ("Duplicated parsing key is ",lpkl,"\n")

        ############################################################################
        # Check for duplicate block parser names and block parsing keywords in this scope
        ############################################################################

        for index,lbpl in enumerate(self.localBlockParserList):
            for lbpl2 in self.localBlockParserList[index+1:]:
                if lbpl == lbpl2:
                    print ("ERROR!! There are duplicate BLOCK parser names in the <<",subLineSource,">> scope")
                    print ("Duplicated name is ",lbpl,"\n")

        for index,lbpkl in enumerate(self.localBlockParsingKeyList):
            for lbpkl2 in self.localBlockParsingKeyList[index+1:]:
                if lbpkl == lbpkl2:
                    print ("ERROR!! There are duplicate BLOCK parsing keys in the <<",subLineSource,">> scope")
                    print ("Duplicated parsing key is ",lbpkl,"\n")

        #Create the parser libraries in the parimary parsers space
        self.cPL = createParserLibrary()
        self.cPL.writeParserLibraryClass(self.targetPath+"/"+self.parserBlockName,self.parserBlockName,self.localParserList,self.localBlockParserList)


#######################################################################################################
##  extractArguments:  this extracts the required and optional arguments from the interpreter lines 
##                     that will later map to xml parameter lists 
#######################################################################################################

    def extractArguments(self,line):
        #tokenize the required line
        lineTokens = line.split()
        xmlRequiredArgument = []
        xmlRequiredArgumentIndex = []
        tokenCounter = 0
        for token in lineTokens:
            if token.find("{") > -1:
                xmlRequiredArgument.append(token)
                xmlRequiredArgumentIndex.append(tokenCounter-1)
            tokenCounter += 1
        return (xmlRequiredArgument,xmlRequiredArgumentIndex)

#######################################################################################################
##  createBlockParserSourceFile:  this generates the new parser python source code and writes it to a 
##                                new charon block parser.
#######################################################################################################

    def createBlockParserSourceFile(self):
        # First, create the parser class file
        nextLine = "\n"
        self.indent = "    "
        self.indent2 = self.indent+self.indent
        self.indent3 = self.indent2+self.indent
        self.indent4 = self.indent3+self.indent
        self.indent5 = self.indent4+self.indent
        self.indent6 = self.indent5+self.indent
        self.indent7 = self.indent6+self.indent
        self.indent8 = self.indent7+self.indent
        self.indent9 = self.indent8+self.indent
        self.indent10 = self.indent9+self.indent
        if self.parserBlockName == "noName":
            print("Error!  No parser name for "+self.blockFileName+"  hase been prescribed")
        self.filename = self.targetPath+"/charonBlockParser"+self.parserBlockName+".py"
        parserFile = open(self.filename,"w+")

        ###################################################################
        #class definition block
        ###################################################################
        #start creating the contents of the file
        fileContents = nextLine
        fileContents += "class charonBlockParser"+self.parserBlockName+":"+nextLine
        fileContents += self.indent+"\""+self.parserBlockName+" parser\""+nextLine
        fileContents += nextLine

        ###################################################################
        #constructor block
        ###################################################################
        fileContents += self.indent+"def __init__(self):"+nextLine
        fileContents += self.indent2+"# Register the parsing keys"+nextLine

         #Add the parsing keys
        fileContents += self.indent2+"self.parsingBlockKey = \""+self.parsingBlockKey+"\""+nextLine
         #Add the name
        fileContents += self.indent2+"self.parserBlockName = \""+self.parserBlockName+"\""+nextLine

        #Add block arguments
        fileContents += nextLine+self.indent2+"# Register the block arguments"+nextLine
        fileContents += self.indent2+"self.blockArgument = []"+nextLine
        for bA in self.blockArgument:
            fileContents += self.indent2+"self.blockArgument.append(\""+bA+"\")"+nextLine
        fileContents += self.indent2+"self.blockArgumentIndexes = []"+nextLine
        for bAI in self.blockArgumentIndexes:
            fileContents += self.indent2+"self.blockArgumentIndexes.append("+str(bAI)+")"+nextLine

        fileContents += self.indent2+"self.interpreterBlockHelpLine = \""+self.interpreterBlockHelpLine+"\""+nextLine

        #Add the xml default lines
        fileContents += nextLine+self.indent2+"# Register the xml default lines"+nextLine
        fileContents += self.indent2+"self.xmlDefaultLines = []"+nextLine
        for xmlDef in self.xmlDefaultLines:
            fileContents += self.indent2+"self.xmlDefaultLines.append(\""+xmlDef+"\")"+nextLine

        #Create a list for the returned xml content and priority codes
        fileContents += nextLine+self.indent2+"self.xmlReturned = []"+nextLine
        fileContents += self.indent2+"self.xmlPriorityCode = []"+nextLine

        ###################################################################
        #isThisMe block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def isThisMe(self,tokenizer,line):"+nextLine
        fileContents += self.indent2+"# Tokenize the line"+nextLine
#        fileContents += self.indent2+"lineTokens = line.split()"+nextLine
        fileContents += self.indent2+"lineTokens = tokenizer.tokenize(line)"+nextLine
        fileContents += self.indent2+"# Tokenize the parsing key"+nextLine
        fileContents += self.indent2+"parsingTokens = self.parsingBlockKey.split()"+nextLine
        fileContents += self.indent2+"returnType = True"+nextLine
        fileContents += self.indent2+"for itoken in range(len(parsingTokens)):"+nextLine
        fileContents += self.indent3+"if itoken+1 > len(lineTokens):"+nextLine
        fileContents += self.indent4+"return False"+nextLine
        fileContents += self.indent3+"if lineTokens[itoken].lower() != parsingTokens[itoken].lower():"+nextLine
        fileContents += self.indent4+"returnType = False"+nextLine
        fileContents += self.indent2+"return returnType"+nextLine

        #######################################################################################################
        ##  Return the block parser name
        #######################################################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def getName(self):"+nextLine
        fileContents += self.indent2+"# Return block parser name"+nextLine
        fileContents += self.indent2+" return self.parserBlockName"+nextLine

        ###################################################################
        #Generate block help line
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def getHelpLine(self):"+nextLine
        fileContents += self.indent2+"return self.interpreterBlockHelpLine"+nextLine

        ###################################################################
        #GenerateXML block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def generateXML(self,line):"+nextLine
        fileContents += self.indent2+"for xmlLine in self.xmlDefaultLines:"+nextLine
        fileContents += self.indent3+"self.xmlReturned.append(xmlLine)"+nextLine
        fileContents += self.indent3+"self.xmlPriorityCode.append(1) #optional lines have priority code 1"+nextLine

        fileContents += nextLine+self.indent2+"return (self.xmlReturned,self.xmlPriorityCode)"+nextLine

 
        ###################################################################
        #Generate bulk replacements block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def generateBulkReplacements(self,tokenizer,line):"+nextLine
#        fileContents += self.indent2+"lineTokens = line.split()"+nextLine
        fileContents += self.indent2+"lineTokens = tokenizer.tokenize(line)"+nextLine
        fileContents += self.indent2+"self.ArgReturned = []"+nextLine
        fileContents += self.indent2+"self.ArgReturnedValue = []"+nextLine
        fileContents += self.indent2+"for xmlLine in range(len(self.blockArgument)):"+nextLine
        fileContents += self.indent3+"self.ArgReturned.append(self.blockArgument[xmlLine])"+nextLine
        fileContents += self.indent3+"self.ArgReturnedValue.append(lineTokens[int(self.blockArgumentIndexes[xmlLine])])"+nextLine

        fileContents += nextLine+self.indent2+"return (self.ArgReturned,self.ArgReturnedValue)"+nextLine


        parserFile.write(fileContents)
