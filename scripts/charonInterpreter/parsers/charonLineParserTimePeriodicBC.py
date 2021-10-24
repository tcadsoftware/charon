
from __future__ import print_function
import copy


class charonLineParserTimePeriodicBC:
    "TimePeriodicBC parser"

    def __init__(self):
        # Register the parsing keys
        self.parserName = "TimePeriodicBC"
        self.parsingKey = "bc is time periodic for"
        self.parsingKeyOptional = []
        self.parsingKeyOptional.append("offset by")
        self.parsingKeyOptional.append("with")
        self.parsingKeyOptional.append("plus")
        self.parsingKeyOptional.append("shifted by")
        self.interpreterHelpLine = "BC is time periodic for {sidesetID} on {geometryBlock} [[with {amplitude1} sin {frequency1} [plus {amplitude2} sin {frequency2} [shifted by {phaseShift}]]] offset by {DCoffset}] "
        self.interpreterQuickHelp = "Specify a temporally periodic potential applied on a contact."
        self.interpreterLongHelp = "Specify a temporally periodic potential applied on a contact  which can be specified as a superposition of two sinusoids and a DC voltage. <> sidesetID is the contact name/type <> geometryBlock is the geometry name the contact is attached to <> amplitude1 and amplitude2 are the amplitudes of the applied sinusoids  in Volts. <> frequency1 and frequency2 are the frequencies of the applied sinusoids  in Hz. <> phaseShift is the phase shift of the second applied sinusoid relative to the first sinusoid  as a percentage of 2\pi. <> DCoffset is the DC voltage applied."

        # Register the xml required lines
        self.xmlRequiredLines = []
        self.xmlRequiredLinePriority = []
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Type,string,Dirichlet")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Sideset ID,string,{sidesetID}")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Element Block ID,string,{geometryBlock}")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Equation Set Name,string,ALL_DOFS")
        self.xmlRequiredLines.append("Charon->Boundary Conditions->{sidesetID}ANONYMOUS,Strategy,string,Sinusoid")
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
        self.xmlRequiredArgumentIndexes = []
        self.xmlRequiredArgumentIndexes.append("5")
        self.xmlRequiredArgumentIndexes.append("7")

        # Register the xml optional lines
        self.xmlOptionalLines = [[]]
        self.xmlOptionalLinePriority = [[]]
        self.xmlOptionalLines[0].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,DC Offset,double,{phaseShift}")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Amplitude 1,double,{amplitude1}")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Frequency 1,double,{frequency1}")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Phase Shift 1,double,0.0")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Amplitude 2,double,0.0}")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Frequency 2,double,0.0")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Phase Shift 2,double,0.0")
        self.xmlOptionalLines[1].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Frequency 2,double,0.0")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[2].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Amplitude 2,double,{amplitude2}")
        self.xmlOptionalLines[2].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Frequency 2,double,{frequency2}")
        self.xmlOptionalLines.append([])
        self.xmlOptionalLines[3].append(" Charon->Boundary Conditions->{sidesetID}ANONYMOUS->Data,Phase Shift 2,double,{phaseShift}")
        self.xmlOptionalLinePriority[0].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority[1].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority[2].append(2)
        self.xmlOptionalLinePriority.append([])
        self.xmlOptionalLinePriority[3].append(2)

        # Register the xml optional arguments and their indexes
        self.xmlOptionalArgument = [['{DCoffset}'], ['{amplitude1}', '{frequency1}'], ['{amplitude2}', '{frequency2}'], ['{phaseShift}']]
        self.xmlOptionalArgumentIndexes = [[2], [1, 3], [1, 3], [2]]

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
