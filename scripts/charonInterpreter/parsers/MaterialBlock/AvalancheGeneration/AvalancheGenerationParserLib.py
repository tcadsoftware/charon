
try:
    import coloramaDISABLED as colors
except ImportError:
    class stubColors:
        "subs for colors when colors doesn't exist on system"
    
        def __init__(self):
            self.Fore = colorClass()
            self.Back = colorClass()
            self.Style = styleClass()
    
    class colorClass():
        "stubbed color class"
    
        def __init__(self):
            self.BLACK = ""
            self.BLUE = ""
            self.WHITE = ""
            self.RED = ""
            self.GREEN = ""
    
    class styleClass():
        "stubbed style class"
    
        def __init__(self):
            self.RESET_ALL = ""
    colors = stubColors()

import sys
from .charonLineParserElectron_ah import *
from .charonLineParserElectron_al import *
from .charonLineParserHole_bh import *
from .charonLineParserHole_hbaromega import *
from .charonLineParserMinField import *
from .charonLineParserElectron_bh import *
from .charonLineParserHole_bl import *
from .charonLineParserElectron_bl import *
from .charonLineParserThreshBehavior import *
from .charonLineParserElectron_hbaromega import *
from .charonLineParserDrivingForce import *
from .charonLineParserholeDFRefDens import *
from .charonLineParserElectron_E0 import *
from .charonLineParserHole_al import *
from .charonLineParserHole_ah import *
from .charonLineParserElectronDFRefDens import *
from .charonLineParserHole_E0 import *
from .charonBlockParserVanOverstraeten import *
from .charonBlockParserCrowellSze import *
from .charonBlockParserMedici import *
from .VanOverstraeten.VanOverstraetenParserLib import *
from .CrowellSze.CrowellSzeParserLib import *
from .Medici.MediciParserLib import *



class AvalancheGenerationParserLib:
    "This is the  AvalancheGenerationParserLib parser library "


    def __init__(self):
        # set the parser library name 
        self.parserLibName = "AvalancheGenerationParserLib"
        # create the linparser objects 
        self.lineParsers = []
        self.lineParsers.append(charonLineParserElectron_ah())
        self.lineParsers.append(charonLineParserElectron_al())
        self.lineParsers.append(charonLineParserHole_bh())
        self.lineParsers.append(charonLineParserHole_hbaromega())
        self.lineParsers.append(charonLineParserMinField())
        self.lineParsers.append(charonLineParserElectron_bh())
        self.lineParsers.append(charonLineParserHole_bl())
        self.lineParsers.append(charonLineParserElectron_bl())
        self.lineParsers.append(charonLineParserThreshBehavior())
        self.lineParsers.append(charonLineParserElectron_hbaromega())
        self.lineParsers.append(charonLineParserDrivingForce())
        self.lineParsers.append(charonLineParserholeDFRefDens())
        self.lineParsers.append(charonLineParserElectron_E0())
        self.lineParsers.append(charonLineParserHole_al())
        self.lineParsers.append(charonLineParserHole_ah())
        self.lineParsers.append(charonLineParserElectronDFRefDens())
        self.lineParsers.append(charonLineParserHole_E0())
        # create the blockparser objects 
        self.blockParsers = []
        self.blockParsers.append([charonBlockParserVanOverstraeten(),VanOverstraetenParserLib()])
        self.blockParsers.append([charonBlockParserCrowellSze(),CrowellSzeParserLib()])
        self.blockParsers.append([charonBlockParserMedici(),MediciParserLib()])
        # create the parserLibrary objects 
        parserLibraries = []
        parserLibraries.append(VanOverstraetenParserLib())
        parserLibraries.append(CrowellSzeParserLib())
        parserLibraries.append(MediciParserLib())


    def isThisMyLine(self,tokenizer,line):
        for lP in self.lineParsers:
            self.isThisMe = lP.isThisMe(tokenizer,line)
            if self.isThisMe == True:
                return (True,lP)
        return (False,None)


    def isThisMyBlock(self,tokenizer,line):
        for bP in self.blockParsers:
            self.isThisMe = bP[0].isThisMe(tokenizer,line)
            if self.isThisMe == True:
                return (True,bP[0],bP[1])
        return (False,None,None)


    def generateHelp(self,genHelp,indent):
        self.addIndent = "     "
        cRStyle = ""
        for lP in self.lineParsers:
            (self.helpLine,self.helpContent) = lP.getHelp(genHelp)
            self.helpContentList = self.helpContent.split("<>")
            print (cRStyle+indent+colors.Fore.RED+colors.Back.WHITE+self.helpLine)
            cRStyle = "\n"
            for hCL in self.helpContentList:
                print ("\t"+indent+colors.Fore.BLUE+colors.Back.WHITE+hCL.lstrip())
        for bP in range(len(self.blockParsers)):
            print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().lstrip())
            self.blockParsers[bP][1].generateHelp(genHelp,indent+self.addIndent)
            print (indent+colors.Fore.GREEN+colors.Back.WHITE+self.blockParsers[bP][0].getHelpLine().replace("start","end").lstrip())
            print (indent+colors.Style.RESET_ALL)


    def getName(self):
        return self.parserLibName