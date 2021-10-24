
from __future__ import print_function
import copy


class charonLineParserConstantCurrentContact:
    "ConstantCurrentContact parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "ConstantCurrentContact"
        self.parsingKey = "bc is current for"
        self.parsingKeyOptional = []
        self.parsingKeyOptional.append("with area")
        self.parsingKeyOptional.append("and length")
        self.parsingKeyOptional.append("with base doping type")
        self.interpreterHelpLine = "BC is current for {sidesetID} on {geometryBlock} fixed at {current} and initial voltage {voltage} [with area {cntarea} [ and length {cntlength} [ with base doping type {baseType} ]]] "
        self.interpreterQuickHelp = "Specify a current boundary condition"
        self.interpreterLongHelp = "Specify a current boundary condition  <> sidesetID is the contact name/type <> geometryBlock is the geometry name the contact is attached to <> fixed at current in A <> with initial voltage in volts"

        # Register the xml required lines
        self.xmlRequiredLines = []
        self.xmlRequiredLinePriority = []
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Type,string,Dirichlet")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Sideset ID,string,{sidesetID}")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Element Block ID,string,{geometryBlock}")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Equation Set Name,string,ALL_DOFS")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Strategy,string,Constant Current")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Current Value,double,{current}")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Initial Voltage,double,{voltage}")
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlRequiredLinePriority.append(2)
        self.xmlNewRequiredLines = []

        # Register the xml required arguments and their indexes
        self.xmlRequiredArgument = []
        self.xmlRequiredArgument.append("{sidesetID}")
        self.xmlRequiredArgument.append("{geometryBlock}")
        self.xmlRequiredArgument.append("{current}")
        self.xmlRequiredArgument.append("{voltage}")
        self.xmlRequiredArgumentIndexes = []
        self.xmlRequiredArgumentIndexes.append("4")
        self.xmlRequiredArgumentIndexes.append("6")
        self.xmlRequiredArgumentIndexes.append("9")
        self.xmlRequiredArgumentIndexes.append("13")

        # Register the xml optional lines
        self.xmlOptionalLines = [[]]
        self.xmlOptionalLinePriority = [[]]
        self.xmlOptionalLines[0].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Device Contact Area,double,{cntarea}")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Simulation Contact Length,double,{cntlength}")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[2].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,BJT1D Base Doping Type,string,{baseType}")
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[2].append(2)

        # Register the xml optional arguments and their indexes
        self.xmlOptionalArgument = [['{cntarea}'], ['{cntlength}'], ['{baseType}']]
        self.xmlOptionalArgumentIndexes = [[2], [2], [4]]

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
