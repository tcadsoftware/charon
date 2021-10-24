
class charonBlockParserCrowellSze:
    "CrowellSze parser"

    def __init__(self):
        # Register the parsing keys
        self.parsingBlockKey = "start crowell-sze parameters"
        self.parserBlockName = "CrowellSze"

        # Register the block arguments
        self.blockArgument = []
        self.blockArgumentIndexes = []
        self.interpreterBlockHelpLine = "start crowell-sze parameters "

        # Register the xml default lines
        self.xmlDefaultLines = []

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
