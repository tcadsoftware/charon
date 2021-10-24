
class charonInitialConditionsModifier0:
    "class for modifying the charonInitialConditionsModifier0 parameterList"


    def __init__(self):
        self.modifierName = "charonInitialConditionsModifier0"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        pLList[:] = [l.replace("LATTICE_TEMPERATURE,Value,double","Lattice Temperature,Value,double") for l in pLList]
        return(pLList)
