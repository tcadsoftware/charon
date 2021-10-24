from __future__ import print_function

class createLineParser:
    "create line parser from interpreter input and xml"

#######################################################################################################
##  Create line parser constructor
#######################################################################################################

    def __init__(self,filename,sourcePathName,targetPathName,verbosity):
        self.filename = filename
        self.sourcePath = sourcePathName
        self.targetPath = targetPathName
        self.parserName = "noName"
        self.lines = list(open(self.sourcePath+"/"+filename))
        self.parsingKey = "Undefined"
        self.parsingKeyOptional = []
        self.verbosity = int(verbosity)
        self.interpreterHelpLine = "No help section exists"
        self.interpreterQuickHelp = ""
        self.interpreterLongHelp = ""
        #######################################################################################
        # Required lines and arguments are ones which will ALWAYS appear in the parameter list 
        # and must be enabled by a line in the interpreter file
        #######################################################################################
        self.xmlRequiredLines = []
        self.xmlRequiredLinePriority = []
        self.xmlRequiredArgument = []
        self.xmlRequiredArgumentIndexes = []
        #######################################################################################
        # Optional lines and arguments are ones which will SOMETIMES appear in the parameter 
        # list and will insert WITH REPLACEMENT all requried and default lines based on a 
        # line in the interpreter file
        #######################################################################################
        self.xmlOptionalLines = []
        self.xmlOptionalLinePriority = []
        self.xmlOptionalLinesTemp = []
        self.xmlOptionalLinePriorityTemp = []
        self.xmlOptionalArgument = []
        self.xmlOptionalArgumentIndexes = []
        #######################################################################################
        # Default lines are ones which will ALWAYS appear in the parameter list and will 
        # insert WITHOUT REPLACEMENT all requried and optional lines based on a line in the 
        # interpreter file
        #######################################################################################
        self.xmlDefaultLines = []
        self.xmlDefaultLinePriority = []
        if self.verbosity > 0:
            class_name = self.__class__.__name__
            print (class_name, self.filename,"created")


#######################################################################################################
##  Create line parser destructor
#######################################################################################################

    def __del__(self):
        if self.verbosity > 3:
            class_name = self.__class__.__name__
            print (class_name, self.filename,"destroyed")

#######################################################################################################
##  Return the parser name
#######################################################################################################

    def getParserName(self):
        return self.parserName;

#######################################################################################################
##  Return the parsing key
#######################################################################################################
 
    def getParsingKey(self):
        return self.parsingKey;
 
#######################################################################################################
##  parseInputs:  parses the input file that defines the new line parser
#######################################################################################################

    def parseInputs(self):
        for line in self.lines:
            lineTokens = line.split()
            foundLineToProcess = False
            #Read in the interpreter lines
            if len(lineTokens) != 0:

                if lineTokens[0].lower() == "interpreter":
                    line = line.replace(","," ")
                    lineTokens = line.split()
               #If this is an interpreter line, take it apart
                if lineTokens[0].lower() == "interpreter":
                    #If this line includes the parser name, grab it and continue on
                    if lineTokens[1].lower() == "name":
                        if len(lineTokens) > 3:
                            print ("Error: parser name cannot contain spaces. \n Bad file is ",self.filename)
                        self.parserName = lineTokens[2]
                        foundLineToProcess = True
                        continue
                    elif lineTokens[1].lower() == "shorthelp":
                        self.interpreterQuickHelp = line.partition('{')[-1].rpartition('}')[0]
                        foundLineToProcess = True
                        continue
                    elif lineTokens[1].lower() == "longhelp":
                        self.interpreterLongHelp = line.partition('{')[-1].rpartition('}')[0]
                        foundLineToProcess = True
                        continue
                    elif lineTokens[1].lower() == "inputline":
                        #Save off the line for help content
                        self.interpreterHelpLine = ""
                        for lT in lineTokens[2:]:
                            self.interpreterHelpLine += lT+" "
                        # Remove the trailing line break
                        self.interpreterHelpLine = self.interpreterHelpLine.rstrip('\n')
                        # Strip out hte ()s
                        self.interpreterHelpLine = self.interpreterHelpLine.replace(")","").replace("(","")
                        foundLineToProcess = True
                    else:
                        print ("\nError: I got the interpreter keyword, but don't understand the following keyword in,",self.filename," on line: \n",line," generateInterpreter is exiting.\n")
                        exit(1)


                    #If not a name, extract any options
                    interpreterOptionalLines = []
                    (interpreterLine,interpreterOptionalLines) = self.extractOptions(line)
                     
                    #extract the parsing keys from the main line and the options
                    (self.parsingKey,self.parsingKeyOptional) = self.extractKeys(interpreterLine,interpreterOptionalLines)

                    #extract the arguments
                    (self.xmlRequiredArgument,self.xmlRequiredArgumentIndexes,self.xmlOptionalArgument,self.xmlOptionalArgumentIndexes) = self.extractArguments(interpreterLine,interpreterOptionalLines)


            #Read in the xml lines
            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "xmlrequired":
                    if lineTokens[-2].lower() == "priority":
                        self.xmlRequiredLinePriority.append(lineTokens[-1])
                        line = line.replace(lineTokens[-1],"")
                        line = line.replace(lineTokens[-2],"")
                    else:
                        self.xmlRequiredLinePriority.append("2")
                    self.xmlRequiredLines.append((line.replace(lineTokens[0]+" ","")).rstrip())
                    foundLineToProcess = True

            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "xmldefault":
                    if lineTokens[-2].lower() == "priority":
                        self.xmlDefaultLinePriority.append(lineTokens[-1])
                        line = line.replace(lineTokens[-1],"")
                        line = line.replace(lineTokens[-2],"")
                    else:
                        self.xmlDefaultLinePriority.append("1")
                    self.xmlDefaultLines.append((line.replace(lineTokens[0]+" ","")).rstrip())
                    foundLineToProcess = True

            if len(lineTokens) != 0:
                if lineTokens[0].lower() == "xmloptional":
                    if lineTokens[-2].lower() == "priority":
                        self.xmlOptionalLinePriorityTemp.append(lineTokens[-1])
                        line = line.replace(lineTokens[-1],"")
                        line = line.replace(lineTokens[-2],"")
                    else:
                        self.xmlOptionalLinePriorityTemp.append("2")
                    self.xmlOptionalLinesTemp.append((line.replace(lineTokens[0]+" ","")).rstrip())
                    
                    foundLineToProcess = True

            if len(lineTokens) != 0:
                if lineTokens[0][0] != "#":
                    if foundLineToProcess == False:
                        print ("\nError:  I found nothing to process on this line in ",self.filename,":\n",line," generateInterpreter is exiting.\n")
                        exit(1)
        self.organizeOptionalLinesAndStripOptionalKeys()



#######################################################################################################
##  OptionalLineReorganize:  this moves a list of optional lines into list of lists of optional lines.
#######################################################################################################


    def optionalLineReorganize(self):
        optTemp = self.xmlOptionalLines
        



#######################################################################################################
##  extractOptions:  this extracts the optional parts of an interpreter line from the whole and returns
##                   the main required line and a list of the options
#######################################################################################################


    def extractOptions(self,line):
        #Extract the compound options line
        compoundOptions = line.partition('[')[-1].rpartition(']')[0]
        interpreterLine = line.replace("["+compoundOptions+"]","")
        #strip off the interpreter directive token
        interpreterLine = interpreterLine.replace(interpreterLine.split()[0]+" ","")

        optionLines = []
        while len(compoundOptions) > 0:
            remainingOptions = compoundOptions.partition('[')[-1].rpartition(']')[0]
            optionLines.append(compoundOptions.replace("["+remainingOptions+"]",""))
            compoundOptions = remainingOptions
            
        return (interpreterLine,optionLines)


#######################################################################################################
##  extractKeys:  this extracts the parsing keys from the required and optional interpreter lines
#######################################################################################################

    def extractKeys(self,interpreterLine,interpreterOptional):
        parsingKey = interpreterLine.partition('(')[-1].rpartition(')')[0].lower()
        if parsingKey == "":
            print("Error!  Cannot extract a parsing key for \n",interpreterLine," in ",self.filename,". Check the input file.")
        parsingKeyOptional = []
        for option in interpreterOptional:
            parsingKeyOptional.append(option.partition('(')[-1].rpartition(')')[0].lower())
        return (parsingKey,parsingKeyOptional)

#######################################################################################################
##  stripOptionalKeys:  this strips the parsing keys from the optional xml lines
#######################################################################################################

    def organizeOptionalLinesAndStripOptionalKeys(self):
        maxOptions = len(self.parsingKeyOptional)
        #initialize the lists
        for iopt in range(maxOptions):
            self.xmlOptionalLines.append([])
            self.xmlOptionalLinePriority.append([])
        oLineCounter=-1
        for xmlOLine in self.xmlOptionalLinesTemp:
            oLineCounter += 1
            parsingKey = "("+xmlOLine.partition('(')[-1].partition(')')[0]+")"
            for ioptKey in range(len(self.parsingKeyOptional)):
                if parsingKey.replace("(","").replace(")","").lower() == self.parsingKeyOptional[ioptKey]:
                    self.xmlOptionalLines[ioptKey].append(xmlOLine.replace(parsingKey,""))
                    self.xmlOptionalLinePriority[ioptKey].append(self.xmlOptionalLinePriorityTemp[oLineCounter])


#######################################################################################################
##  extractArguments:  this extracts the required and optional arguments from the interpreter lines 
##                     that will later map to xml parameter lists 
#######################################################################################################

    def extractArguments(self,interpreterLine,interpreterOptional):
        #tokenize the required line
        lineTokens = interpreterLine.split()
        xmlRequiredArgument = []
        xmlRequiredArgumentIndex = []
        tokenCounter = 0
        for token in lineTokens:
            if token.find("{") > -1:
                xmlRequiredArgument.append(token)
                xmlRequiredArgumentIndex.append(tokenCounter-1)
            tokenCounter += 1

        #Loop over the list of optional lines and extract the arguments
        xmlOptionalArgument = []
        xmlOptionalArgumentIndex = []
        for oline in interpreterOptional:
            #tokenize the option line
            lineTokens = oline.split()
            localOptionalArgument = []
            localOptionalArgumentIndex = []
            tokenCounter = 0
            for token in lineTokens:
                if token.find("{") > -1:
                    localOptionalArgument.append(token)
                    localOptionalArgumentIndex.append(tokenCounter)
                tokenCounter += 1
            xmlOptionalArgument.append(localOptionalArgument)
            xmlOptionalArgumentIndex.append(localOptionalArgumentIndex)

        #Return the lists
        return (xmlRequiredArgument,xmlRequiredArgumentIndex,xmlOptionalArgument,xmlOptionalArgumentIndex)


#######################################################################################################
##  createLineParserSourceFile:  this generates the new parser python source code and writes it to a 
##                               new charon line parser.
#######################################################################################################

    def createLineParserSourceFile(self):
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
        if self.parserName == "noName":
            print("Error!  No parser name for "+self.filename+" hase been prescribed")
        self.filename = self.targetPath+"/charonLineParser"+self.parserName+".py"
        parserFile = open(self.filename,"w+")

        ###################################################################
        #import the copy module for deep copies
        ###################################################################
        fileContents = nextLine
        fileContents += "from __future__ import print_function"+nextLine
        fileContents += "import copy"+nextLine
        fileContents += nextLine


        ###################################################################
        #class definition block
        ###################################################################
        #start creating the contents of the file
        fileContents += nextLine
        fileContents += "class charonLineParser"+self.parserName+":"+nextLine
        fileContents += self.indent+"\""+self.parserName+" parser\""+nextLine
        fileContents += nextLine

        ###################################################################
        #constructor block
        ###################################################################
        fileContents += self.indent+"def __init__(self):"+nextLine
        fileContents += self.indent2+"# Register the parsing keys"+nextLine

        # Add parser name
        fileContents += self.indent2+"self.parserName = \""+self.parserName+"\""+nextLine
        #Add the parsing keys
        fileContents += self.indent2+"self.parsingKey = \""+self.parsingKey+"\""+nextLine
        fileContents += self.indent2+"self.parsingKeyOptional = []"+nextLine
        for optionalKey in self.parsingKeyOptional:
            fileContents += self.indent2+"self.parsingKeyOptional.append(\""+optionalKey+"\")"+nextLine

        #Add help related strings
        fileContents += self.indent2+"self.interpreterHelpLine = "+"\""+self.interpreterHelpLine+"\""+nextLine
        fileContents += self.indent2+"self.interpreterQuickHelp = "+"\""+self.interpreterQuickHelp+"\""+nextLine
        fileContents += self.indent2+"self.interpreterLongHelp = "+"\""+self.interpreterLongHelp+"\""+nextLine


        #Add the xml required lines
        fileContents += nextLine+self.indent2+"# Register the xml required lines"+nextLine
        fileContents += self.indent2+"self.xmlRequiredLines = []"+nextLine
        fileContents += self.indent2+"self.xmlRequiredLinePriority = []"+nextLine
        for xmlReq in self.xmlRequiredLines:
            fileContents += self.indent2+"self.xmlRequiredLines.append(\""+xmlReq+"\")"+nextLine
        for xmlReqP in self.xmlRequiredLinePriority:
            fileContents += self.indent2+"self.xmlRequiredLinePriority.append("+xmlReqP+")"+nextLine

        fileContents += self.indent2+"self.xmlNewRequiredLines = []"+nextLine

        #Add the xml required arguments
        fileContents += nextLine+self.indent2+"# Register the xml required arguments and their indexes"+nextLine
        fileContents += self.indent2+"self.xmlRequiredArgument = []"+nextLine
        for xmlReqArg in self.xmlRequiredArgument:
            fileContents += self.indent2+"self.xmlRequiredArgument.append(\""+xmlReqArg+"\")"+nextLine
        fileContents += self.indent2+"self.xmlRequiredArgumentIndexes = []"+nextLine
        for xmlReqArg in self.xmlRequiredArgumentIndexes:
            fileContents += self.indent2+"self.xmlRequiredArgumentIndexes.append(\""+str(xmlReqArg)+"\")"+nextLine

        #Add the xml optional lines
        fileContents += nextLine+self.indent2+"# Register the xml optional lines"+nextLine
        fileContents += self.indent2+"self.xmlOptionalLines = [[]]"+nextLine
        fileContents += self.indent2+"self.xmlOptionalLinePriority = [[]]"+nextLine
        optCounter=0
        for xmlOpt in self.xmlOptionalLines:
            for xmlOptInner in xmlOpt:
                fileContents += self.indent2+"self.xmlOptionalLines["+str(optCounter)+"].append(\""+xmlOptInner+"\")"+nextLine
            optCounter += 1
            if optCounter < len(self.xmlOptionalLines):
                fileContents += self.indent2+"self.xmlOptionalLines.append([])"+nextLine
        #Add priorities
        optCounter=0
        for xmlOptP in self.xmlOptionalLinePriority:
            for xmlOptPInner in xmlOptP:
                fileContents += self.indent2+"self.xmlOptionalLinePriority["+str(optCounter)+"].append("+xmlOptPInner+")"+nextLine
            optCounter += 1
            if optCounter < len(self.xmlOptionalLinePriority):
                fileContents += self.indent2+"self.xmlOptionalLinePriority.append([])"+nextLine

        #Add the xml optional arguments
        fileContents += nextLine+self.indent2+"# Register the xml optional arguments and their indexes"+nextLine
        fileContents += self.indent2+"self.xmlOptionalArgument = "+str(self.xmlOptionalArgument)+nextLine
        fileContents += self.indent2+"self.xmlOptionalArgumentIndexes = "+str(self.xmlOptionalArgumentIndexes)+nextLine


        #Add the xml default lines
        fileContents += nextLine+self.indent2+"# Register the xml default lines"+nextLine
        fileContents += self.indent2+"self.xmlDefaultLines = []"+nextLine
        fileContents += self.indent2+"self.xmlDefaultLinePriority = []"+nextLine
        for xmlDef in self.xmlDefaultLines:
            fileContents += self.indent2+"self.xmlDefaultLines.append(\""+xmlDef+"\")"+nextLine
        for xmlDefP in self.xmlDefaultLinePriority:
            fileContents += self.indent2+"self.xmlDefaultLinePriority.append("+xmlDefP+")"+nextLine

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
        fileContents += self.indent2+"parsingTokens = self.parsingKey.split()"+nextLine
        fileContents += self.indent2+"returnType = True"+nextLine
        fileContents += self.indent2+"for itoken in range(len(parsingTokens)):"+nextLine
        fileContents += self.indent3+"if itoken+1 > len(lineTokens):"+nextLine
        fileContents += self.indent4+"return False"+nextLine
        fileContents += self.indent3+"if lineTokens[itoken].lower() != parsingTokens[itoken].lower():"+nextLine
        fileContents += self.indent4+"returnType = False"+nextLine
        fileContents += self.indent2+"return returnType"+nextLine


        ###################################################################
        #getName block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def getName(self):"+nextLine
        fileContents += self.indent2+"# Return parser name"+nextLine
        fileContents += self.indent2+" return self.parserName"+nextLine


        ###################################################################
        #getHelp block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def getHelp(self,verbosity):"+nextLine
        fileContents += self.indent2+"# Return help content"+nextLine
        fileContents += self.indent2+"if verbosity.lower() == \"long\":"+nextLine
        fileContents += self.indent3+"return (self.interpreterHelpLine,self.interpreterLongHelp)"+nextLine
        fileContents += self.indent2+"else:"+nextLine
        fileContents += self.indent3+"return (self.interpreterHelpLine,self.interpreterQuickHelp)"+nextLine

        ###################################################################
        #GenerateXML block
        ###################################################################
        fileContents += nextLine+nextLine+nextLine
        fileContents += self.indent+"def generateXML(self,tokenizer,line):"+nextLine
        fileContents += self.indent2+"# Tokenize the line"+nextLine
#        fileContents += self.indent2+"lineTokens = line.split()"+nextLine
        fileContents += self.indent2+"lineTokens = tokenizer.tokenize(line)"+nextLine
        fileContents += self.indent2+"self.xmlNewRequiredLines[:] = []"+nextLine
        fileContents += self.indent2+"for xL in self.xmlRequiredLines:"+nextLine
        fileContents += self.indent3+"self.xmlNewRequiredLines.append(xL)"+nextLine

        fileContents += self.indent2+"for ipar in range(len(self.xmlRequiredArgument)):"+nextLine
        fileContents += self.indent3+"line.replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])"+nextLine
        #required xml content
        fileContents += self.indent3+"for iRLine in range(len(self.xmlRequiredLines)):"+nextLine
        fileContents += self.indent4+"self.xmlNewRequiredLines[iRLine]=self.xmlNewRequiredLines[iRLine].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])"+nextLine
        fileContents += self.indent2+"for index,xmlLine in enumerate(self.xmlNewRequiredLines):"+nextLine
        fileContents += self.indent3+"self.xmlReturned.append(xmlLine)"+nextLine
        fileContents += self.indent3+"self.xmlPriorityCode.append(self.xmlRequiredLinePriority[index]) #required lines have priority code 2"+nextLine

        fileContents += self.indent2+"# Look over input line to see if any options are called out."+nextLine
        fileContents += self.indent2+"optCounter = 0"+nextLine
        fileContents += self.indent2+"optIndex = 0"+nextLine
        fileContents += self.indent2+"for optKey in self.parsingKeyOptional:"+nextLine
        fileContents += self.indent3+"# Tokenize the opt keys"+nextLine
        fileContents += self.indent3+"foundOptionalKey = False"+nextLine
        fileContents += self.indent3+"optKeyTokens = optKey.split()"+nextLine
        fileContents += self.indent3+"for iLT in range(len(lineTokens)):"+nextLine
        fileContents += self.indent4+"if lineTokens[iLT].lower() == optKeyTokens[0]:"+nextLine
        fileContents += self.indent5+"if len(optKeyTokens) == 1:"+nextLine
        fileContents += self.indent6+"optIndex = iLT"+nextLine
        fileContents += self.indent6+"foundOptionalKey = True"+nextLine
        fileContents += self.indent5+"else:"+nextLine

        fileContents += self.indent6+"for iPK in range(len(optKeyTokens)-1):"+nextLine
        fileContents += self.indent7+"optIndex = iLT"+nextLine
        fileContents += self.indent7+"if iLT+iPK+1 > len(lineTokens)-1:"+nextLine
        fileContents += self.indent8+"continue"+nextLine
        fileContents += self.indent7+"if optKeyTokens[iPK+1] == lineTokens[iLT+iPK+1].lower():"+nextLine
        fileContents += self.indent8+"if iPK+2 == len(optKeyTokens):"+nextLine
        fileContents += self.indent9+"foundOptionalKey = True"+nextLine
        fileContents += self.indent8+"else:"+nextLine
        fileContents += self.indent9+"continue"+nextLine


        fileContents += self.indent3+"#Found the key, now create the xml line"+nextLine
        fileContents += self.indent3+"if foundOptionalKey == True:"+nextLine
        fileContents += self.indent4+"self.Returned=copy.deepcopy(self.xmlOptionalLines[optCounter])"+nextLine

        fileContents += self.indent4+"for iopt in range(len(self.xmlOptionalLines[optCounter])):"+nextLine
        fileContents += self.indent5+"for ipar in range(len(self.xmlOptionalArgument[optCounter])):"+nextLine
        fileContents += self.indent6+"self.Returned[iopt] = self.Returned[iopt].replace(self.xmlOptionalArgument[optCounter][ipar],lineTokens[optIndex+int(self.xmlOptionalArgumentIndexes[optCounter][ipar])])"+nextLine
        fileContents += self.indent5+"for ipar in range(len(self.xmlRequiredArgument)):"+nextLine
        fileContents += self.indent6+"self.Returned[iopt] = self.Returned[iopt].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])"+nextLine
        fileContents += self.indent5+"self.xmlReturned.append(self.Returned[iopt])"+nextLine
        fileContents += self.indent5+"self.xmlPriorityCode.append(2) #optional lines have priority code 2"+nextLine
        fileContents += self.indent3+"optCounter += 1"+nextLine

        fileContents += self.indent2+"for xmlLine in self.xmlDefaultLines:"+nextLine
        fileContents += self.indent3+"self.xmlReturned.append(xmlLine)"+nextLine
        fileContents += self.indent3+"self.xmlPriorityCode.append(1) #optional lines have priority code 1"+nextLine

        fileContents += nextLine+self.indent2+"return (self.xmlReturned,self.xmlPriorityCode)"+nextLine

        parserFile.write(fileContents)
