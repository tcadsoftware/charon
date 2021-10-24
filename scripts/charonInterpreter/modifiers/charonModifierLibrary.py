
from .charonContactOnInsulatorModifier0 import *
from .charon1DBJTBaseModifier0 import *
from .charonInitialConditionsModifier0 import *
from .charonInitialConditionsModifier1 import *
from .charonOhmicBCModifier0 import *
from .charonFixedChargeModifier0 import *
from .charonChargeDensityModifier0 import *
from .charonVoltageModifier0 import *
from .charonHarmonicBalanceAnalysisModifier0 import *
from .charonHarmonicBalancePhaseShiftsModifier0 import *
from .charonHarmonicBalanceFrequenciesModifier0 import *
from .charonHarmonicBalanceAmplitudesModifier0 import *


class modifierLibrary:
    "modifier library"


    def __init__(self):
        # create the modifier objects
        self.modifierList = []
        self.modifierList.append(charonContactOnInsulatorModifier0())
        self.modifierList.append(charon1DBJTBaseModifier0())
        self.modifierList.append(charonInitialConditionsModifier0())
        self.modifierList.append(charonInitialConditionsModifier1())
        self.modifierList.append(charonOhmicBCModifier0())
        self.modifierList.append(charonFixedChargeModifier0())
        self.modifierList.append(charonChargeDensityModifier0())
        self.modifierList.append(charonVoltageModifier0())
        self.modifierList.append(charonHarmonicBalanceAnalysisModifier0())
        self.modifierList.append(charonHarmonicBalancePhaseShiftsModifier0())
        self.modifierList.append(charonHarmonicBalanceFrequenciesModifier0())
        self.modifierList.append(charonHarmonicBalanceAmplitudesModifier0())


    def executeModifiers(self,useModList,pLList):
        # Loop over the requested modifiers
        for uML in useModList:
            for mL in self.modifierList:
                if uML == mL.getName():
                    pLList = mL.testForModification(pLList)

        return pLList
