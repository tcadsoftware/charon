
class charon1DBJTBaseModifier1:
    "class for modifying the charon1DBJTBaseModifier1 parameterList"


    def __init__(self):
        self.modifierName = "charon1DBJTBaseModifier1"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        print "This is just a test for having multiple modifiers in a single parser"
        return pLList
