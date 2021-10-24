
class charonInitialConditionsModifier1:
    "class for modifying the charonInitialConditionsModifier1 parameterList"


    def __init__(self):
        self.modifierName = "charonInitialConditionsModifier1"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        pLList[:] = [l.replace("LATTICE_TEMPERATURE,Value,string","Lattice Temperature,Value,string") for l in pLList]
        return(pLList)
