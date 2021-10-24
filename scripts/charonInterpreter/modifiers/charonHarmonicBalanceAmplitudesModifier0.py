
class charonHarmonicBalanceAmplitudesModifier0:
    "class for modifying the charonHarmonicBalanceAmplitudesModifier0 parameterList"


    def __init__(self):
        self.modifierName = "charonHarmonicBalanceAmplitudesModifier0"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        for lineNumber, line in enumerate(pLList):
            if line.find("AMPLITUDES") >= 0:
                # the line to replace will look like this:
                # 'AMPLITUDES1e1   3e1   2e2AMPLITUDES'
                # and it should become
                # '1e1,3e1,2e2'
                [begin, replacement, end] = pLList[lineNumber].split('AMPLITUDES')
                replacement = ','.join(replacement.split())
                pLList[lineNumber] = begin+replacement+end
        return(pLList)
