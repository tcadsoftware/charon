

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


    def getEnclosedString(self,bracketLeft,bracketRight,line):
        #If there are quoted parts in the line, separate them out.
        quotedSection = line.lstrip().partition(bracketLeft)[-1].rpartition(bracketRight)[0]
        if len(quotedSection) == 0:
            return quotedSection,-1

        preQuote = line.split("\"")
        preQuoteTokens = preQuote[0].split()
        tokenLocation = len(preQuoteTokens)

        return quotedSection,tokenLocation


