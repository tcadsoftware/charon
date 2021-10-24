from __future__ import print_function

#import _gv
#from _gv import digraph
#from _gv import graph
from os import system

from parsers.parsersParserLib import *
from charonParameterList import *
from modifiers.charonModifierLibrary import *
from charonTokenize import *


class interpreterParser:
    "Charon interpreter parser"


    #######################################################################################################
    ##  parser constructor
    #######################################################################################################

    def __init__(self):
        self.parsersLibrary = parsersParserLib()
        self.parameterLists = []
        self.parameterListsPriority = []
        self.parserLibraryObjectPointer = []


    #######################################################################################################
    ##  Get the listed lines of the interpreter file
    #######################################################################################################
    def parseInterpreterFile(self, interpreterLines,plFilename,verbosity):
        self.verbosity = verbosity
        self.parserLibraryObjectPointer.append(self.parsersLibrary)
        ## Proceed through each line and query the parsers for a match.
        self.tempPL = []
        self.tempPr = []
        self.parameterListReplacements = []
        self.parameterListReplacementValues = []
        self.PLModifiers = []
        self.setErrorCode = False
        self.inCommentBlock = False
        tokenizer = charonTokenize()
        for line in interpreterLines:
            lineTokens = line.split() 
            if len(lineTokens) != 0:           # ignore blank lines
                if self.verbosity >= 10:
                    print("    ",line.strip())
                if lineTokens[0] == "/#":
                    if self.inCommentBlock == True:
                        print ("Misplaced open comment block /#.  Fix please!")
                    elif self.inCommentBlock == False:
                        self.inCommentBlock = True
                        continue
                if lineTokens[0] == "#/":
                    if self.inCommentBlock == False:
                        print ("Misplaced end comment block #/.  Fix please!")
                        exit(1)
                    elif self.inCommentBlock == True:
                        self.inCommentBlock = False
                        continue

                           
                if lineTokens[0][0] != "#" and self.inCommentBlock == False:    # ignore # commented lines
                    #check line keys first
                    self.isItHere = False
                    self.parserObjectPointer = None
                    (self.isItHere,self.parserObjectPointer) = self.parserLibraryObjectPointer[-1].isThisMyLine(tokenizer,line)
                    if self.isItHere == True:
                        if len(self.tempPL) != 0:
                            self.tempPL[:] = []  #I really hate this.  Python 3 would be so much better
                        if len(self.tempPr) != 0:
                            self.tempPr[:] = []  #I really hate this.  Python 3 would be so much better
                        (self.tempPL,self.tempPr) = self.parserObjectPointer.generateXML(tokenizer,line)
                        self.checkForModifiers(self.parserObjectPointer.getName())
                        if self.verbosity >= 20:
                            print (line.strip()," returns the following PLs")
                            for plLine in self.tempPL:
                                print ("     ",plLine)

                        # Do a quick sanity check on line parser return
                        if self.isItHere == True and len(self.tempPL) == 0:
                            self.setErrorCode = True
                            print ("Something has gone wrong. \n I found a parser, but got no parameters for \n",line,"It could be malformed.")

                        self.replaceReplacements()
                        self.parameterLists.extend(self.tempPL)
                        self.parameterListsPriority.extend(self.tempPr)
                        continue
 

                    # Check for block parsers
                    (self.isItHere,self.parserObjectPointer,self.tempPLPointer) = self.parserLibraryObjectPointer[-1].isThisMyBlock(tokenizer,line)
                    if self.isItHere == True:
                        if self.verbosity >= 10:
                            print ("Entering block ",self.parserObjectPointer.getName())
                        self.parserLibraryObjectPointer.append(self.tempPLPointer)
                        #print (line,"Entering block ",self.parserLibraryObjectPointer[-1].getName())

                        (self.tempPLR,self.tempPLRV) = self.parserObjectPointer.generateBulkReplacements(tokenizer,line)
                        self.searchAndReplaceReplacements()
                        (self.tempPL,self.tempPr) = self.parserObjectPointer.generateXML(line)
                        #print self.tempPL
                        self.replaceReplacements()
                        #print self.tempPL

                        self.parameterLists.extend(self.tempPL)
                        self.parameterListsPriority.extend(self.tempPr)
                        continue

                    if self.isItHere == False and line.split()[0].lower() == "end":
                        if self.verbosity >= 10:
                            print ("Exiting block ")
                        self.parserLibraryObjectPointer.pop()

                    if self.isItHere == False and line.split()[0].lower() != "end":
                        self.setErrorCode = True
                        print ("Something has gone wrong. \n I can't find a parser for \n",line)

            if self.setErrorCode == True:
                print ("An error condition has been set.  Terminate interpreter.")
                quit()

        #Run through the modifiers if there are any
        plModifier = modifierLibrary()
        plModifier.executeModifiers(self.PLModifiers,self.parameterLists)

        paramList = self.constructXMLParameterList()

        if self.verbosity >= 15:
            for pl in self.parameterLists:
                print (pl)

        #paramList.printList("")
        if plFilename != "NoName":
            filehandle = open(plFilename,"w+")
            paramList.writeXMLList("",filehandle)

        #help(_gv)
        #help(digraph)
        #help(graph)

        #self.createGraph()




    #######################################################################################################
    ##  search and replace block replacement values returned by block parser
    #######################################################################################################
    def searchAndReplaceReplacements(self):
        for tempPLI in range(len(self.tempPLR)):
            foundAndReplaced = False
            for replI in range(len(self.parameterListReplacements)):
                if self.parameterListReplacements[replI] == self.tempPLR[tempPLI]:
                    self.parameterListReplacementValues[replI] = self.tempPLRV[tempPLI]
                    foundAndReplaced = True
            if foundAndReplaced == False:  # add it to the list
                self.parameterListReplacements.append(self.tempPLR[tempPLI])
                self.parameterListReplacementValues.append(self.tempPLRV[tempPLI])


    #######################################################################################################
    ##  replace block replacement values
    #######################################################################################################
    def replaceReplacements(self):
        for replI in range(len(self.parameterListReplacements)):
            for pli in range(len(self.tempPL)):
                self.tempPL[pli] = self.tempPL[pli].replace(self.parameterListReplacements[replI],self.parameterListReplacementValues[replI])


    #######################################################################################################
    ##  Construct XML formatted parameter Lists
    #######################################################################################################
    def constructXMLParameterList(self):
        #First order of business is to find the depth of priority
        localParamList = charonParameterList("Charon")
        maxPriority = 1
        for pLP in self.parameterListsPriority:
            if int(pLP) > maxPriority:
                maxPriority = int(pLP)

        for pLP in range(maxPriority):
            for lPLP in range(len(self.parameterListsPriority)):
                if int(self.parameterListsPriority[lPLP]) == pLP+1:
                    pLTokens = self.parameterLists[lPLP].split(",")
                    #This is a little dumb, but when there are comma delimited parameter values, we need to rejoin part of the list
                    if len(pLTokens) > 4:
                        temp = pLTokens[:3]
                        pLTokens = ",".join(pLTokens[3:])
                        temp.append(pLTokens)
                        pLTokens = temp
                                 
                    try: 
                        localParamList.insertParameterWithListCreation(pLTokens[0].split("->"),pLTokens[1],pLTokens[2],pLTokens[3],0)
                    except IndexError:
                        print ("bad index in line ",self.parameterLists[lPLP])
                        sys.exit(1)
                    except Exception as e:
                        print ("Failed on line ",self.parameterLists[lPLP])
                        sys.exit(1)

        #localParamList.resolveAnonymous()

        return localParamList

    #######################################################################################################
    ##  check for modifiers
    #######################################################################################################
    def checkForModifiers(self,modifierName):
        for plNum, pl in enumerate(self.tempPL):
            plTokens = pl.split()
            if len(plTokens) >= 3:
                if plTokens[-3].lower() == "use" and plTokens[-2].lower() == "modifier":
                    # Add this modifier to the list
                    nameToAdd = "charon"+modifierName+"Modifier"+plTokens[-1]
                    self.PLModifiers.append(nameToAdd)

        #Remove all modifiers from the parameter list and parameter priority list
        for plNum, pl in reversed(list(enumerate(self.tempPL))):
            plTokens = pl.split()
            if len(plTokens) >= 3:
                if plTokens[-3].lower() == "use" and plTokens[-2].lower() == "modifier":
                    del self.tempPL[plNum]
                    del self.tempPr[plNum]


    #######################################################################################################
    ##  Genreate the short help output
    #######################################################################################################
    def generateHelp(self,generateHelpFlag):
        #self.parserLibraryObjectPointer.append(self.parsersLibrary)
        self.parsersLibrary.generateHelp(generateHelpFlag,"")


    #######################################################################################################
    ##  create graph
    #######################################################################################################
    def createGraph(self):
        graph = digraph("input")
        graph.edge('Hello', 'World')

        graph.view()
