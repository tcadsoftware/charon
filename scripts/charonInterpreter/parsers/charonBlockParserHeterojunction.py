
class charonBlockParserHeterojunction:
    "Heterojunction parser"

    def __init__(self):
        # Register the parsing keys
        self.parsingBlockKey = "start heterojunction"
        self.parserBlockName = "Heterojunction"

        # Register the block arguments
        self.blockArgument = []
        self.blockArgument.append("{sidesetName}")
        self.blockArgumentIndexes = []
        self.blockArgumentIndexes.append(3)
        self.interpreterBlockHelpLine = "start heterojunction on {sidesetName} "

        # Register the xml default lines
        self.xmlDefaultLines = []
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole1ANONYMOUS,Type,string,Interface")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole1ANONYMOUS,Strategy,string,Interface Heterojunction")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole1ANONYMOUS,Sideset ID,string,{sidesetName}")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole1ANONYMOUS,Equation Set Name,string,HOLE_DENSITY1")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole1ANONYMOUS,Equation Set Name2,string,HOLE_DENSITY2")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole2ANONYMOUS,Type,string,Interface")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole2ANONYMOUS,Strategy,string,Interface Heterojunction")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole2ANONYMOUS,Sideset ID,string,{sidesetName}")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole2ANONYMOUS,Equation Set Name,string,HOLE_DENSITY2")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->hole2ANONYMOUS,Equation Set Name2,string,HOLE_DENSITY1")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec1ANONYMOUS,Type,string,Interface")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec1ANONYMOUS,Strategy,string,Interface Heterojunction")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec1ANONYMOUS,Sideset ID,string,{sidesetName}")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec1ANONYMOUS,Equation Set Name,string,ELECTRON_DENSITY1")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec1ANONYMOUS,Equation Set Name2,string,ELECTRON_DENSITY2")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec2ANONYMOUS,Type,string,Interface")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec2ANONYMOUS,Strategy,string,Interface Heterojunction")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec2ANONYMOUS,Sideset ID,string,{sidesetName}")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec2ANONYMOUS,Equation Set Name,string,ELECTRON_DENSITY2")
        self.xmlDefaultLines.append("Charon->Boundary Conditions->elec2ANONYMOUS,Equation Set Name2,string,ELECTRON_DENSITY1")

        self.xmlReturned = []
        self.xmlPriorityCode = []



    def isThisMe(self,tokenizer,line):
        # Tokenize the line
        lineTokens = tokenizer.tokenize(line)
        # Tokenize the parsing key
        parsingTokens = self.parsingBlockKey.split()
        returnType = True
        for itoken in range(len(parsingTokens)):
            if itoken+1 > len(lineTokens):
                return False
            if lineTokens[itoken].lower() != parsingTokens[itoken].lower():
                returnType = False
        return returnType



    def getName(self):
        # Return block parser name
         return self.parserBlockName



    def getHelpLine(self):
        return self.interpreterBlockHelpLine



    def generateXML(self,line):
        for xmlLine in self.xmlDefaultLines:
            self.xmlReturned.append(xmlLine)
            self.xmlPriorityCode.append(1) #optional lines have priority code 1

        return (self.xmlReturned,self.xmlPriorityCode)



    def generateBulkReplacements(self,tokenizer,line):
        lineTokens = tokenizer.tokenize(line)
        self.ArgReturned = []
        self.ArgReturnedValue = []
        for xmlLine in range(len(self.blockArgument)):
            self.ArgReturned.append(self.blockArgument[xmlLine])
            self.ArgReturnedValue.append(lineTokens[int(self.blockArgumentIndexes[xmlLine])])

        return (self.ArgReturned,self.ArgReturnedValue)
