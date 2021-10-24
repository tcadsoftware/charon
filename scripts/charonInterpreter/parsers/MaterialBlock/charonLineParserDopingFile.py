
from __future__ import print_function
import copy


class charonLineParserDopingFile:
    "DopingFile parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "DopingFile"
        self.parsingKey = "read doping from"
        self.parsingKeyOptional = []
        self.parsingKeyOptional.append("acceptor")
        self.parsingKeyOptional.append("donor")
        self.parsingKeyOptional.append("in 2d")
        self.parsingKeyOptional.append("in 3d")
        self.parsingKeyOptional.append("with buffer")
        self.interpreterHelpLine = "read doping from {filename} for [acceptor [donor [ in 2d [ in 3d [ with buffer {bufferValue}]]]]] "
        self.interpreterQuickHelp = "Read the doping from a tabulated file of doping concentrations"
        self.interpreterLongHelp = "Read the doping from a tabulated file of doping concentrations <> {filename} is the name of the file to be read <> Specify acceptor or donor <> Optionally add a buffer of {bufferValue} around the max and min spatial extents of the doping specifications."

        # Register the xml required lines
        self.xmlRequiredLines = []
        self.xmlRequiredLinePriority = []
        self.xmlRequiredLines.append("Charon->Closure Models->{MaterialBlockName}->Doping,Value,string,Function")
        self.xmlRequiredLines.append("Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},File Name,string,{filename}")
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlNewRequiredLines = []

        # Register the xml required arguments and their indexes
        self.xmlRequiredArgument = []
        self.xmlRequiredArgument.append("{filename}")
        self.xmlRequiredArgumentIndexes = []
        self.xmlRequiredArgumentIndexes.append("3")

        # Register the xml optional lines
        self.xmlOptionalLines = [[]]
        self.xmlOptionalLinePriority = [[]]
        self.xmlOptionalLines[0].append(" Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},Doping Type,string,Acceptor")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[1].append(" Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},Doping Type,string,Donor")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[2].append(" Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},Function Type,string,File2D")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[3].append(" Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},Function Type,string,File3D")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[4].append(" Charon->Closure Models->{MaterialBlockName}->Doping->Function {filename},Buffer,double,{bufferValue}")
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[3].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[4].append(2)

        # Register the xml optional arguments and their indexes
        self.xmlOptionalArgument = [[], [], [], [], ['{bufferValue}']]
        self.xmlOptionalArgumentIndexes = [[], [], [], [], [2]]

        # Register the xml default lines
        self.xmlDefaultLines = []
        self.xmlDefaultLinePriority = []

        self.xmlReturned = []
        self.xmlPriorityCode = []



    def isThisMe(self,tokenizer,line):
        # Tokenize the line
        lineTokens = tokenizer.tokenize(line)
        # Tokenize the parsing key
        parsingTokens = self.parsingKey.split()
        returnType = True
        for itoken in range(len(parsingTokens)):
            if itoken+1 > len(lineTokens):
                return False
            if lineTokens[itoken].lower() != parsingTokens[itoken].lower():
                returnType = False
        return returnType



    def getName(self):
        # Return parser name
         return self.parserName



    def getHelp(self,verbosity):
        # Return help content
        if verbosity.lower() == "long":
            return (self.interpreterHelpLine,self.interpreterLongHelp)
        else:
            return (self.interpreterHelpLine,self.interpreterQuickHelp)



    def generateXML(self,tokenizer,line):
        # Tokenize the line
        lineTokens = tokenizer.tokenize(line)
        self.xmlNewRequiredLines[:] = []
        for xL in self.xmlRequiredLines:
            self.xmlNewRequiredLines.append(xL)
        for ipar in range(len(self.xmlRequiredArgument)):
            line.replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
            for iRLine in range(len(self.xmlRequiredLines)):
                self.xmlNewRequiredLines[iRLine]=self.xmlNewRequiredLines[iRLine].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
        for index,xmlLine in enumerate(self.xmlNewRequiredLines):
            self.xmlReturned.append(xmlLine)
            self.xmlPriorityCode.append(self.xmlRequiredLinePriority[index]) #required lines have priority code 2
        # Look over input line to see if any options are called out.
        optCounter = 0
        optIndex = 0
        for optKey in self.parsingKeyOptional:
            # Tokenize the opt keys
            foundOptionalKey = False
            optKeyTokens = optKey.split()
            for iLT in range(len(lineTokens)):
                if lineTokens[iLT].lower() == optKeyTokens[0]:
                    if len(optKeyTokens) == 1:
                        optIndex = iLT
                        foundOptionalKey = True
                    else:
                        for iPK in range(len(optKeyTokens)-1):
                            optIndex = iLT
                            if iLT+iPK+1 > len(lineTokens)-1:
                                continue
                            if optKeyTokens[iPK+1] == lineTokens[iLT+iPK+1].lower():
                                if iPK+2 == len(optKeyTokens):
                                    foundOptionalKey = True
                                else:
                                    continue
            #Found the key, now create the xml line
            if foundOptionalKey == True:
                self.Returned=copy.deepcopy(self.xmlOptionalLines[optCounter])
                for iopt in range(len(self.xmlOptionalLines[optCounter])):
                    for ipar in range(len(self.xmlOptionalArgument[optCounter])):
                        self.Returned[iopt] = self.Returned[iopt].replace(self.xmlOptionalArgument[optCounter][ipar],lineTokens[optIndex+int(self.xmlOptionalArgumentIndexes[optCounter][ipar])])
                    for ipar in range(len(self.xmlRequiredArgument)):
                        self.Returned[iopt] = self.Returned[iopt].replace(self.xmlRequiredArgument[ipar],lineTokens[int(self.xmlRequiredArgumentIndexes[ipar])])
                    self.xmlReturned.append(self.Returned[iopt])
                    self.xmlPriorityCode.append(2) #optional lines have priority code 2
            optCounter += 1
        for xmlLine in self.xmlDefaultLines:
            self.xmlReturned.append(xmlLine)
            self.xmlPriorityCode.append(1) #optional lines have priority code 1

        return (self.xmlReturned,self.xmlPriorityCode)
