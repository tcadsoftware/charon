
class charonHarmonicBalanceAnalysisModifier0:
    "class for modifying the charonHarmonicBalanceAnalysisModifier0 parameterList"


    def __init__(self):
        self.modifierName = "charonHarmonicBalanceAnalysisModifier0"


    def getName(self):
        return self.modifierName


    def testForModification(self,pLList):
        # if a time domain equation set is specified here, move it to the Frequency Domain Options within the appropriate physics block
        for lineNumber, line in enumerate(pLList):
            if(line.find("Charon->Physics Blocks->") >= 0):
                # grab the physics block name
                physicsBlockName = line.split("->")[2]
                pLList[lineNumber] = line.replace("Charon->Physics Blocks->"+physicsBlockName+"->ANONYMOUS,Type,string,","Charon->Physics Blocks->"+physicsBlockName+"->ANONYMOUS->Options,Time Domain Equation Set,string,")
        # set the equation set type to Frequency Domain within the appropriate physics block
        for lineNumber, line in enumerate(pLList):
            if(line.find("Charon->Physics Blocks->") >= 0):
                # grab the physics block name
                physicsBlockName = line.split("->")[2]
                pLList[lineNumber] = line.replace("Charon->Physics Blocks->"+physicsBlockName+"->ANONYMOUS,HBType","Charon->Physics Blocks->"+physicsBlockName+"->ANONYMOUS,Type")
        return pLList
