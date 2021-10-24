

class xmlLCMConverter:
    "This is an xml to LCM file converter"

    def __init__(self,filename):
        self.filename = filename
        self.inputLines = list(open(filename))
        self.nestingList = []
        self.LCMParameters = []

####################################################################################
#  Larryfy the xml file 
####################################################################################


    def convertFile(self):
        inCommentMode = False
        for iL in self.inputLines:
            lineTokens = iL.split()

            if len(lineTokens) == 0:
                continue

            if "<!--" in iL and "-->" in iL:  # one line comment
                stripThis = iL.partition('<!--')[-1].rpartition('-->')[0]
                iL = iL.replace("<!--"+stripThis+"-->","")
                lineTokens = iL.split()
                if len(lineTokens) == 0:
                    continue

            if "<!--" in iL:
                inCommentMode = True
                continue

            if "-->" in iL:
                inCommentMode = False
                continue

            if inCommentMode == True:
                continue

            ########################################################################
            #Add a new nest to the list.
            ########################################################################
            if lineTokens[0][:14] == "<ParameterList":
                if lineTokens[-1][-2:] != "/>":
                    c = '"'
                    quoteLocations = [pos for pos, char in enumerate(iL) if char == c]
                    if len(quoteLocations) == 0:
                        nestName = "ANONYMOUS"
                    else:
                        nestName = iL[quoteLocations[0]+1:quoteLocations[1]]
                    self.nestingList.append(nestName)

            ########################################################################
            #Remove a nest
            ########################################################################
#            if lineTokens[0] == "</ParameterList>":
            if "</ParameterList" in lineTokens[0]:
                if len(self.nestingList) > 0:
                    self.nestingList.pop()
                continue

            ########################################################################
            # extract a parameter
            ########################################################################
            if lineTokens[0] == "<Parameter":
                c = '"'
                parameterTokens = iL.split("=")
                getName = False
                getType = False
                getValue = False
                for pT in parameterTokens:
                    if getName == True:
                        quoteLocations = [pos for pos, char in enumerate(pT) if char == c]
                        parameterName = pT[quoteLocations[0]+1:quoteLocations[1]]
                        getName = False
                    if getType == True:
                        quoteLocations = [pos for pos, char in enumerate(pT) if char == c]
                        typeName = pT[quoteLocations[0]+1:quoteLocations[1]]
                        getType = False
                    if getValue == True:
                        quoteLocations = [pos for pos, char in enumerate(pT) if char == c]
                        try:
                            valueName = pT[quoteLocations[0]+1:quoteLocations[1]]
                        except Exception as e:
                            print ("exception occurred with line ",pT)
                        getValue = False
                    if pT[-4:] == "name":
                        getName = True
                    if pT[-4:] == "type":
                        getType = True
                    if pT[-5:] == "value":
                        getValue = True

                self.addNewParameter(parameterName,typeName,valueName)

####################################################################################
#  Larryfy the nested crap
####################################################################################
    def printNest(self):
        LCMLine = ""
        for nM in self.nestingList:
            LCMLine += nM+"->"

        print (LCMLine)

####################################################################################
#  Add a new parameter to the list
####################################################################################
    def addNewParameter(self,parameterName,typeName,valueName):
        LCMLine = ""
        for nM in self.nestingList:
            LCMLine += nM+"->"

        parameterString = LCMLine[:-2]+","+parameterName+","+typeName+","+valueName

        self.LCMParameters.append(parameterString)


####################################################################################
#  Print out the Larryfied parameter list
####################################################################################
    def printLCMParameters(self):
        for LCMP in self.LCMParameters:
            print (LCMP)

####################################################################################
#  write the Larryfied parameter list to a file
####################################################################################
    def writeLCMParameters(self):
        LCMParameterFile = open(self.filename+".LCMified","w+")
        for LCMP in self.LCMParameters:
            LCMParameterFile.write(LCMP+"\n")


####################################################################################
#  write the Larryfied parameter list to a file
####################################################################################
    def getLCMParameters(self):
        return self.LCMParameters
