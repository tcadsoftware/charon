
class charonOhmicBCModifier1:
    "class for modifying the charonOhmicBCModifier1 parameterList"


    def __init__(self):
        self.modifierName = "charonOhmicBCModifier1"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        for lineNumber, line in pLList:
            print "This is just a test for having multiple modifiers in a single parser",lineNumber, line
