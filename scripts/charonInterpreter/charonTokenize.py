
from itertools import combinations
from copy import deepcopy
import sys
import difflib


class charonTokenize:
    "This is a tokenizer for the charon Interpreter script"


###############################################################
# This is stateless
###############################################################

    def __init__(self):
        #Don't even do anything here
        self.ImAlive = True


###############################################################
# Tokenize a line
###############################################################

    def tokenize(self,line):

        ##########################################################3
        ## This mostly tokenizes on white space--handle cases 
        ## where an = has no space next to it
        ##########################################################3

        line = line.replace("="," = ")

        ##########################################################3
        ## Handle quoted sections
        ##########################################################3

        #Look for quoted sections
        bracketLeft = "\""
        bracketRight = "\""
        (quotedSection,qsTokenLocation) = self.getEnclosedString(bracketLeft,bracketRight,line)

        #Temporarily remove the quoted section from the line
        line = line.replace(bracketLeft+quotedSection+bracketRight,"")
        line = line.replace(","," ")

        #Tokenize the remainder on spaces
        lineTokens = line.split()

        #Finally, insert the quoted section into the appropriate token location
        if qsTokenLocation >-1:
            lineTokens.insert(qsTokenLocation,quotedSection)



        #Return the result
        return lineTokens

###############################################################
# get a single token from enclosed string
###############################################################

    def getEnclosedString(self,bracketLeft,bracketRight,line):
        #If there are quoted parts in the line, separate them out.
        quotedSection = line.lstrip().partition(bracketLeft)[-1].rpartition(bracketRight)[0]
        if len(quotedSection) == 0:
            return quotedSection,-1

        preQuote = line.split("\"")
        preQuoteTokens = preQuote[0].split()
        tokenLocation = len(preQuoteTokens)

        return quotedSection,tokenLocation


###############################################################
# get a set of all the possible option combinations
###############################################################


    def getOptionCombinations(self,possibleOptions):
        #Return all possible combinations of identified possible options

        # It's cheapest just to hardwire the possibilities for a small number of possibilities.
        numPossible = len(possibleOptions)
        possibleCombinations = []
        for index in range(numPossible):
            localComb = list(combinations(range(numPossible),index+1))
            for comb in localComb:
                possibleCombinations.append(list(comb))

        return possibleCombinations

###############################################################
# get a set of all the possible options
###############################################################


    def getOptions(self,line,options):

 
        #Loop over options 
        numPossible = 0
        amIPossible = [False]*len(options)
        for index,opt in enumerate(options):
            optTokens=opt.lower().split()
             #join into string 
            optString = ""
            optStringTokens = []
            for oT in optTokens:
                if "{" not in oT:
                    optStringTokens.append(oT)
                if "{" in oT:
                    break
            optString += ' '.join(optStringTokens)
            if optString in line:
                amIPossible[index] = True
                numPossible += 1


        possibleOptions = []
        possibleOptionsIndex = []
        for index,poss in enumerate(amIPossible):
            if poss == True:
                possibleOptions.append(options[index])
                possibleOptionsIndex.append(index)

        possibleCombinations = self.getOptionCombinations(possibleOptions)

        return (possibleOptions,possibleOptionsIndex,possibleCombinations)




    ###############################################################
    # Count the number if single instance options in a set
    ###############################################################


    def numSizeOneCount(self,optionVecBool):

        numCount = 0
        for oVB in optionVecBool:
            if oVB.count(False) == 1:
                numCount += 1

        return numCount



    ###############################################################
    # Count the number unvalidated options in optionsValidated list
    ###############################################################


    def optionValidatedFalseCount(self,optionVecBool):

        numCount = 0
        for oVB in optionVecBool:
            numCount += oVB.count(False)

        return numCount


    ###############################################################
    # set option validated
    ###############################################################

    def setOptionValidated(self,lineTokensBool,lineIndexes,optionValidated):
        #sufficient to check only the first line token for the option
       
        for lIVecIndex,lIVec in enumerate(lineIndexes):
            for lIIndex,lI in enumerate(lIVec):
                if lineTokensBool[lI] == True:
                    optionValidated[lIVecIndex][lIIndex] = True


    ###############################################################
    # validate line tokens
    ###############################################################

    def validateLineTokens(self,lineTokensBool,start,length):
        for index in range(length):
            lineTokensBool[start+index] = True


    ###############################################################
    # validate line tokens
    ###############################################################

    def getOptionVariables(self,options,indexes,lineTokens,startIndex,optionVariables,optionBool):
        for oIIndex,optionIndex in enumerate(indexes):
            optionBool[optionIndex] = True
            optionTokens = options[optionIndex].split()
            tempVec = []
            tempDict = {}
            counter = -1
            for oT in optionTokens:
                counter += 1
                if "{" in oT:
                    tempVec.append(lineTokens[startIndex[oIIndex]+counter])
                    tempDict[oT] = lineTokens[startIndex[oIIndex]+counter]
            optionVariables.append(tempDict)
            #optionVariables[oIIndex] = tempDict


    ###############################################################
    # get a set of all the possible option combinations
    ###############################################################


    ###############################################################
    # findOptions
    ###############################################################


    def findOptions(self,options,lineTokens,lineTokensLower):
        optionVariables = [{}]*len(options)
        optionBools = [False]*len(options)
        processLine = ' '.join(lineTokensLower)
        line = ' '.join(lineTokensLower)


        (possibleOptions,possibleOptionsIndex,possibleCombinations) = self.getOptions(line,options)

        condition = "not validated"
        # First do a simple check to see if there is only one possible.
        #line = ' '.join(lineTokensLower)

        allComboLineStrings = []
        allComboIndexes = []
        allComboLineTokenIndexes = []
        allComboTokenLengths = []
        lineTokensBool = [False]*len(lineTokensLower)

        for index,combo in enumerate(possibleCombinations):
            comboLineStrings = []
            comboIndexes = []
            comboLineTokenIndexes = []
            comboTokenLengths = []
            for comCheck in combo:

                (lineStrings,lineTokenIndexes) = self.retrieveLineStrings(possibleOptions[comCheck],lineTokensLower,lineTokensBool)

                comboLineStrings.append(lineStrings)
                comboIndexes.append(possibleOptionsIndex[comCheck])
                comboLineTokenIndexes.append(lineTokenIndexes)
                comboTokenLengths.append(len(possibleOptions[comCheck].split()))

            allComboLineStrings.append(comboLineStrings)
            allComboIndexes.append(comboIndexes)
            allComboLineTokenIndexes.append(comboLineTokenIndexes)
            allComboTokenLengths.append(comboTokenLengths)

        validatedCombinations = [False]*len(allComboLineStrings)
        validatedOptionStartIndex = []
        for aCLSIndex,aCLS in enumerate(allComboLineStrings):
            validatedOptionStartIndex.append([-1]*len(aCLS))


        for aCLSIndex,aCLS in enumerate(allComboLineStrings):

            lineTokensBool = [False]*len(lineTokensLower)
            optionValidated = []
            for cLSIndex,cLS in enumerate(aCLS):
                optionValidated.append([False]*len(cLS))

            continueProcessing = True
            safetyCount = 0
 
            failThisCombo = False
            optionsValidated = [False]*len(aCLS)

            while continueProcessing and safetyCount < 3:
                safetyCount += 1

                for cLSIndex,cLS in enumerate(aCLS):

                    numSizeOne = self.numSizeOneCount(optionValidated)

                    #Validated the option
                    (localLineStrings,lineTokenIndexes) = self.retrieveLineStrings(options[allComboIndexes[aCLSIndex][cLSIndex]],lineTokensLower,lineTokensBool)
                    if len(localLineStrings) == 1:
                        optionsValidated[cLSIndex] = True
                        self.validateLineTokens(lineTokensBool,lineTokenIndexes[0],len(options[allComboIndexes[aCLSIndex][cLSIndex]].split()))
                        validatedOptionStartIndex[aCLSIndex][cLSIndex] = lineTokenIndexes[0]

                    if lineTokensBool.count(False) == 0 and optionsValidated.count(False) > 0:  #This combo cannot be validated
                        failThisCombo = True
                        continueProcessing = False

                    if lineTokensBool.count(False) == 0 and optionsValidated.count(False) == 0:
                        continueProcessing = False
                        validatedCombinations[aCLSIndex] = True

                    if lineTokensBool.count(False) == 0 and optionsValidated.count(False) > 0: #used up the line but not the options-Fail
                        continueProcessing = False

                    if lineTokensBool.count(False) >0 and optionsValidated.count(False) == 0: #used up the options but not the line-Fail
                        continueProcessing = False


        if validatedCombinations.count(True) == 0:
            condition = "not validated"

        if validatedCombinations.count(True) > 1:
            condition = "ambiguous"

        retOptionVariables = []

        if validatedCombinations.count(True) == 1: # get option variables
            condition = "validated"
            for vCIndex,vC in enumerate(validatedCombinations):
                if vC == True:
                    self.getOptionVariables(options,allComboIndexes[vCIndex],lineTokens,validatedOptionStartIndex[vCIndex],retOptionVariables,optionBools)

        optCounter = 0
        for oBIndex, oB in enumerate(optionBools):
            if oB == True:
                optionVariables[oBIndex] = retOptionVariables[optCounter]
                optCounter += 1

        return (condition,optionBools,optionVariables)

    ################################################################################
    ##  Retrieve line strings
    ################################################################################


    def retrieveLineStrings(self,option,lineTokens,lineTokensBool):
        lineStrings = []
        lineStringsOptionsIndex = []
        lineTokensIndexes = []
        localTokens = deepcopy(lineTokens)
        #First, get the option up to the first variable
        subOption = option.lower()
        #Hamhanded, but works
        subOptionTokens = option.lower().split()
        sanitySubTokens = []
        addToSanity = True
        for sOIndex,sO in enumerate(subOptionTokens):
            if "{" in sO:
                subOptionTokens[sOIndex] = "wildcard"
        lineCheckString = ' '.join(localTokens).lower()
        maxCounter = len(localTokens)+1
        counter = 0
        remainingOptions = True 
        while remainingOptions == True and len(localTokens)>0 and counter<maxCounter:
            counter += 1

            for index,lT in enumerate(localTokens):
                thisLineString = ""
                thisLineIndex = -1
                #create wildcarded line string tokens
                if index+len(subOptionTokens) > len(localTokens): #can do no more with this line
                    remainingOptions = False
                    break
                compTokens = deepcopy(localTokens[index:index+len(subOptionTokens)])
                indexTokens = deepcopy(localTokens[index:index+len(subOptionTokens)])
                for lTIndex,lT in enumerate(compTokens):
                    if subOptionTokens[lTIndex] == "wildcard":
                        compTokens[lTIndex] = "wildcard"

                if subOptionTokens == compTokens:
                    if len(lineTokensIndexes) == 0:
                        startIndex = 0
                    else:
                        startIndex = lineTokensIndexes[-1]+1
                    thisLineIndex = self.getListTokensStartIndex(lineTokens,indexTokens,startIndex)
                    if thisLineIndex == None:  #cou'd not find matching line
                        continue
                    if lineTokensBool[thisLineIndex:thisLineIndex+len(option.split())].count(True) > 0:
                        continue
                    thisLineString = ' '.join(localTokens[index:index+len(option.split())]).strip()
                    del localTokens[index:index+len(option.split())]
                    break
            if thisLineString is not "":
                lineStrings.append(thisLineString)
                lineTokensIndexes.append(thisLineIndex)

        return (lineStrings,lineTokensIndexes)


    ################################################################################
    ##  Remove lower tokens
    ################################################################################

    def removeLowerTokens(self,rString,options,lowerTokens):
        rTokens = rString

    ################################################################################
    ##  Find the index of sub tokens
    ################################################################################

    def getListTokensStartIndex(self,lineTokens,lineSubTokens,startIndex=0):
        returnIndex = None
        for index,lT in enumerate(lineTokens[startIndex:]):
            if lineSubTokens == lineTokens[index+startIndex:index+startIndex+len(lineSubTokens)]:
                returnIndex = index+startIndex
                break

        return returnIndex




 
    ################################################################################
    ##  IsItMe
    ################################################################################

    def isItMe(self,line,primary,options):
        #Process out things we don't need and buffer things like =
        primary = primary.replace("="," = ")
        primaryVariables = {}
        optionVariables = []*len(options)
        suggestedSyntax = ""

        foundOptions = []
        for index,opt in enumerate(options):
            foundOptions.append(False)
            # remove () for this test
            options[index] = options[index].replace("="," = ")
            
        #First step.  Buffer lien "=" and tokenize with both lower case and case preserved
        if "=" in line:
            line = line.replace("="," = ")
        lineTokens = self.tokenize(line)
        lineTokensLower = self.tokenize(line.lower())

        #check the line for primary
        primaryTokens = primary.split()


        affirmPrimary = True
        localVariables = []

        for index,pT in enumerate(primaryTokens):
            if "{" in pT:  #This is a variable and amounts to a wildcard for this check
                localVariables.append(lineTokens[index]) #Need to preserve case for variables
                primaryVariables[pT] = lineTokens[index]
                continue
            if pT != lineTokensLower[index]:
                affirmPrimary = False

                #Assess for similarity
                for pTIndex,pT in enumerate(primaryTokens):
                    if pTIndex == len(lineTokensLower):
                        break
                    if "{" not in pT:
                        primaryExpected = pT.lower()+" "
                        primaryCheck = lineTokensLower[pTIndex]+" "

                primaryExpected = primaryExpected.strip()
                primaryCheck = primaryCheck.strip()
                match = difflib.SequenceMatcher(a=primaryExpected, b=primaryCheck).ratio()
                if match > 0.85:
                    suggestedSyntax = primary+"\n"+"with options:"
                    for opt in options:
                        suggestedSyntax += "\n       "+opt

                return (affirmPrimary,foundOptions,primaryVariables,optionVariables,suggestedSyntax)

        #Check for options
        if len(lineTokens) == len(primaryTokens):
            lineTokens = []
        if len(lineTokens) > len(primaryTokens): # sanity check
            lineTokens = lineTokens[len(primaryTokens):]
            lineTokensLower = lineTokensLower[len(primaryTokens):]

            #(condition,optionVariables) = self.findOptions(lineTokens,lineTokensLower)
            (condition,foundOptions,optionVariables) = self.findOptions(options,lineTokens,lineTokensLower)

            if condition == "not validated" or condition == "ambiguous":
                print ()
                print ("Error!  This line validated for the primary, but chirp cannot resolve the options")
                print ("\n",line.strip(),"\n")
                print ("The validated primary is:")
                print ("    ",primary,"\n")
                print ("The valid options are:")
                for opt in options:
                    print("     ",opt)
                sys.exit(1)


        return (affirmPrimary,foundOptions,primaryVariables,optionVariables,suggestedSyntax)

