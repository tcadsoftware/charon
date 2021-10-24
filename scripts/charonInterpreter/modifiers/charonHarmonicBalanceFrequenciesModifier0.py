
class charonHarmonicBalanceFrequenciesModifier0:
    "class for modifying the charonHarmonicBalanceFrequenciesModifier0 parameterList"


    def __init__(self):
        self.modifierName = "charonHarmonicBalanceFrequenciesModifier0"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        for lineNumber, line in enumerate(pLList):
            if line.find("FREQUENCIES") >= 0:
                # the line to replace will look like this:
                # 'FREQUENCIES1e1   3e1   2e2FREQUENCIES'
                # and it should become
                # '1e1,3e1,2e2'
                [begin, replacement, end] = pLList[lineNumber].split('FREQUENCIES')
                replacement = ','.join(replacement.split())
                pLList[lineNumber] = begin+replacement+end
        return(pLList)
